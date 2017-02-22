#include "surface.h"
#include "STL_reader.h"


Surface::Surface(std::string filename, double pumpingFactor, double temperature,
    std::string label, bool flipNormals, double scale, double edgeLength)
    : pumpingFactor_(pumpingFactor), temperature_(temperature),
      pumpedParticles_(0), collisionCounter_(0), label_(label) {
    loadFromSTL(filename, scale, edgeLength);
    std::cout << "Building AABB tree...\n";
    buildAABBTree();
    std::cout << "Computing area cumulative distribution function...\n";
    computeAreaCDF();
    std::cout << "Computing face normals...\n";
    computeFaceNormals();
    computeFaceRotations(flipNormals);
}

void Surface::setAccommodationCoefficients(
    std::unordered_map<Element, double, element_hash> coefficients) {
    accommodationCoefficients_ = coefficients;
}

double Surface::getAccommodationCoefficient(Element element) const {
    std::unordered_map<Element, double, element_hash>::const_iterator coeff_it =
        accommodationCoefficients_.find(element);
    if (coeff_it == accommodationCoefficients_.end()) {
        return DEFAULT_ACCOMMODATION_COEFFICIENT;
    }

    return coeff_it->second;
}

void Surface::computeFirstIntersection(const Ray& r,
    Ray_intersection& isect, bool useSkipFunctor, const Skip &skip) const {
    if (useSkipFunctor) {
        isect = tree_.first_intersection(r, skip);
    } else {
        isect = tree_.first_intersection(r);
    }
}

void Surface::computeAllIntersections(const Ray& r,
    std::back_insert_iterator<std::vector<Ray_intersection>> it) const {
    tree_.all_intersections(r, it);
}

bool Surface::isLoaded() const {
    return mesh_.is_valid() && !mesh_.is_empty();
}

double Surface::getTemperature() const {
    return temperature_;
}

unsigned long Surface::getPumpedParticles() const {
    return pumpedParticles_;
}

bool Surface::checkIfPumped(Rng& rng) const {
    return uni01(rng) < pumpingFactor_;
}

void Surface::addPumpedParticle() {
    pumpedParticles_ += 1;
}

void Surface::buildAABBTree() {
    if (!isLoaded()) return;

    tree_.rebuild(faces(mesh_).first, faces(mesh_).second, mesh_);
    tree_.build();
    // This is not probably needed as we do only intersection queries ATM
    // tree_.accelerate_distance_queries();
    bbox_ = tree_.bbox();
}

void Surface::computeFaceNormals() {
    if (!isLoaded()) return;

    faceNormals_ = mesh_.add_property_map<face_descriptor, Vector>
        ("f:normals", CGAL::NULL_VECTOR).first;
    CGAL::Polygon_mesh_processing::compute_face_normals(mesh_,
        faceNormals_);
}

void Surface::computeFaceRotations(bool flipNormals) {
    if (!isLoaded()) return;

    K::Aff_transformation_3 ident;
    faceRotations_ = mesh_.add_property_map<face_descriptor, K::Aff_transformation_3>
        ("f:rotations", ident).first;

    for (face_descriptor fd : mesh_.faces()) {
        halfedge_descriptor i = mesh_.halfedge(fd);
        Point v1 = mesh_.point(mesh_.source(i)),
              v2 = mesh_.point(mesh_.source(mesh_.next(i)));

        Vector t1 = v2 - v1;
        t1 = t1 / std::sqrt(t1.squared_length());
        Vector n = faceNormals_[fd];
        if (flipNormals) n = n * -1.0;
        Vector t2 = CGAL::cross_product(n, t1);
        t2 = t2 / std::sqrt(t2.squared_length());

        faceRotations_[fd] = K::Aff_transformation_3(
            t1.x(), t2.x(), n.x(),
            t1.y(), t2.y(), n.y(),
            t1.z(), t2.z(), n.z()
        );
    }
}

void Surface::computeAreaCDF() {
    if (!isLoaded()) return;

    areaCDF_.resize(mesh_.number_of_faces());

    double cdfValue = 0.0;
    for (face_descriptor fd : mesh_.faces()) {
        double faceArea = CGAL::Polygon_mesh_processing::face_area(fd, mesh_,
            CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh_.points()));
        cdfValue += faceArea;
        areaCDF_[fd] = cdfValue;
    }

    for (face_descriptor fd : mesh_.faces()) {
        areaCDF_[fd] /= cdfValue;
    }
}

bool Surface::loadFromSTL(std::string filename, double scale,
    double edgeLength) {
    std::ifstream ifs;
    Surface_mesh mesh;

    ifs.open(filename, std::ifstream::in);

    if (!ifs.is_open()) {
        throw std::invalid_argument("Couldn't open STL file");
        return false;
    }

    std::cout << "Loading surface " << filename << "...\n";
    std::vector<std::array<double, 3>> points;
    std::vector<std::array<int, 3>> facets;

    bool ret = CGAL::read_STL(ifs, points, facets, true);
    std::cout << "File has " << points.size() << " points and " <<
        facets.size() << " facets\n";

    if (scale != 1.0) {
        for (size_t i = 0; i < points.size(); ++i) {
            points[i][0] *= scale;
            points[i][1] *= scale;
            points[i][2] *= scale;
        }
    }

    // Remove degenerate triangles
    size_t nRemovedTriangles = 0;
    for (std::vector<std::array<int, 3>>::iterator it = facets.begin();
        it < facets.end(); ) {
        size_t i0 = (*it)[0], i1 = (*it)[1], i2 = (*it)[2];
        Point p0(points[i0][0], points[i0][1], points[i0][2]);
        Point p1(points[i1][0], points[i1][1], points[i1][2]);
        Point p2(points[i2][0], points[i2][1], points[i2][2]);
        Triangle t(p0, p1, p2);
        if (t.is_degenerate()) {
            facets.erase(it);
            nRemovedTriangles += 1;
        } else {
            ++it;
        }
    }
    std::cout << nRemovedTriangles << " degenerate triangles removed\n";

    if (ret) {
        std::cout << "orient_polygon_soup...\n";
        ret =
            CGAL::Polygon_mesh_processing::orient_polygon_soup(points, facets);
    }
    if (ret) {
        if (CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(facets)) {
            std::cout << "The input is a polygon mesh\n";
        } else {
            throw std::invalid_argument("The input is not a polygon mesh");
        }
    }

    if (ret) {
        std::cout << "polygon_soup_to_polygon_mesh...\n";
        CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points,
            facets, mesh);
        int nRemoved = CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
        std::cout << nRemoved << " isolated vertices removed\n";
        ret = mesh.is_valid(true) && !mesh.is_empty();
    }

    if (ret) {
        std::cout << "Mesh assembly ok\n";
    } else {
        throw std::invalid_argument("Mesh assembly failed");
    }

    if (ret) {
        if (edgeLength > 0.0) {
            std::cout << "Splitting long edges in the mesh...\n";
            CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh),
                edgeLength, mesh);
        }
        std::cout << "Loaded " << filename << " with " <<
            mesh.number_of_faces() << " faces and " <<
            mesh.number_of_vertices() << " vertices." << std::endl;

        mesh_ = mesh;
    }

    ifs.close();

    if (!ret) {
        throw std::invalid_argument("Couldn't decode STL file");
    }

    return ret;
}

std::tuple<Point, face_descriptor> Surface::getRandomPoint(Rng& rng) const {
    double rnd = uni01(rng);

#ifdef DEBUG
    size_t nFaces = areaCDF_.size();
    std::cout << "Number of faces: " << nFaces << std::endl <<
        "Random number given: " << rnd << std::endl;
#endif

    size_t faceIndex = std::upper_bound(areaCDF_.begin(), areaCDF_.end(),
        rnd) - areaCDF_.begin();

    face_descriptor fd(faceIndex);
#ifdef DEBUG
    std::cout << "Found index: " << fd << std::endl <<
        "Corresponding areaCDF value: " << areaCDF_[fd] << std::endl;
#endif

    halfedge_descriptor i = mesh_.halfedge(fd), i2 = mesh_.next(i);
    Point v1 = mesh_.point(mesh_.source(i)),
          v2 = mesh_.point(mesh_.source(i2)),
          v3 = mesh_.point(mesh_.source(mesh_.next(i2)));

    double a1 = uni01(rng), a2 = uni01(rng);
    if (a1 + a2 > 1.0) {
        a1 = 1.0 - a1;
        a2 = 1.0 - a2;
    }

    double b = 1.0 - a1 - a2;
    Point p(b*v1.x() + a1*v2.x() + a2*v3.x(),
            b*v1.y() + a1*v2.y() + a2*v3.y(),
            b*v1.z() + a1*v2.z() + a2*v3.z());

#ifdef DEBUG
    std::cout << "Found point: (" << p.x() << ", " << p.y() << ", " <<
        p.z() << ")" << std::endl << std::endl;
#endif

    return std::make_tuple(p, fd);
}

Direction Surface::generateCosineLawDirection(face_descriptor fd, Rng& rng) const {
    double phi = 2.0*M_PI * uni01(rng);
    double theta = std::asin(std::sqrt(uni01(rng)));
    double st = std::sin(theta);

    Direction d(st * std::cos(phi),
                st * std::sin(phi),
                std::cos(theta));

    return faceRotations_[fd].transform(d);
}

Bbox Surface::bbox() const {
    return CGAL::Polygon_mesh_processing::bbox_3(mesh_);
}

std::string Surface::getLabel() const {
    return label_;
}

void Surface::incrementCollisionCounter() {
    collisionCounter_ += 1;
}

unsigned long Surface::getCollisionCount() const {
    return collisionCounter_;
}

