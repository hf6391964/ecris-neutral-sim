#include "surface.h"

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
        faceNormals_,
        CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh_.points()).
        geom_traits(K()));
}

void Surface::computeFaceRotations() {
    if (!isLoaded()) return;

    K::Aff_transformation_3 ident;
    faceRotations_ = mesh_.add_property_map<face_descriptor, K::Aff_transformation_3>
        ("f:rotations", ident).first;

    for (face_descriptor fd : mesh_.faces()) {
        halfedge_descriptor i = mesh_.halfedge(fd), i2 = mesh_.next(i);
        Point v1 = mesh_.point(mesh_.source(i)),
              v2 = mesh_.point(mesh_.source(i2)),
              v3 = mesh_.point(mesh_.source(mesh_.next(i2)));

        Vector t1 = v2 - v1;
        Vector t2 = v3 - v1;
        t1 = t1 / std::sqrt(t1.squared_length());
        t2 = t2 / std::sqrt(t2.squared_length());
        Vector n = CGAL::cross_product(t2, t1);
        n = n / -std::sqrt(n.squared_length());

        faceRotations_[fd] = K::Aff_transformation_3(
            t1.x(), t2.x(), n.x(),
            t1.y(), t2.y(), n.y(),
            t1.z(), t2.z(), n.z()
        );
    }
}

void Surface::computeFaceMidpoints() {
    if (!isLoaded()) return;

    Surface_mesh::Property_map<face_descriptor, Point> faceMidpoints =
        mesh_.add_property_map<face_descriptor, Point>("f:midpoints").first;

    // TODO if needed
}

void Surface::computeAreaCDF() {
    if (!isLoaded()) return;

    double meshArea = CGAL::Polygon_mesh_processing::area(mesh_,
        CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh_.points()).
        geom_traits(K()));

    areaCDF_.resize(mesh_.number_of_faces());

    double cdfValue = 0.0;
    for (face_descriptor fd : mesh_.faces()) {
        double faceArea = CGAL::Polygon_mesh_processing::face_area(fd, mesh_,
            CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh_.points()));
        cdfValue += faceArea / meshArea;
        areaCDF_[fd] = cdfValue;
    }
}

Point Surface::getRandomPoint(Rng &rng) {
    size_t nFaces = areaCDF_.size();

    double rnd = rng.uniform(0.0, 1.0);

    std::cout << "Number of faces: " << nFaces << std::endl <<
        "Random number given: " << rnd << std::endl;

    size_t faceIndex = std::upper_bound(areaCDF_.begin(), areaCDF_.end(),
        rnd) - areaCDF_.begin();

    face_descriptor fd(faceIndex);
    std::cout << "Found index: " << (int)fd << std::endl <<
        "Corresponding areaCDF value: " << areaCDF_[fd] << std::endl;

    halfedge_descriptor i = mesh_.halfedge(fd), i2 = mesh_.next(i);
    Point v1 = mesh_.point(mesh_.source(i)),
          v2 = mesh_.point(mesh_.source(i2)),
          v3 = mesh_.point(mesh_.source(mesh_.next(i2)));

    // TODO validate this
    double a1 = rng.uniform(0.0, 1.0), a2 = rng.uniform(0.0, 1.0);
    if (a1 + a2 > 1.0) {
        a1 = 1.0 - a1;
        a2 = 1.0 - a2;
    }

    double b = 1.0 - a1 - a2;
    Point p(b*v1.x() + a1*v2.x() + a2*v3.x(),
            b*v1.y() + a1*v2.y() + a2*v3.y(),
            b*v1.z() + a1*v2.z() + a2*v3.z());

    std::cout << "Found point: (" << p.x() << ", " << p.y() << ", " <<
        p.z() << ")" << std::endl << std::endl;

    return p;
}

bool Surface::loadFromSTL(std::string filename, double avgTriangleArea) {
    std::ifstream ifs;
    Surface surf;
    Surface_mesh mesh;

    ifs.open(filename, std::ifstream::in);

    if (!ifs.is_open()) return false;

    std::vector<std::array<double, 3>> points;
    std::vector<std::array<int, 3>> facets;

    bool ret = CGAL::read_STL(ifs, points, facets, true);
    if (ret) {
        ret =
            CGAL::Polygon_mesh_processing::orient_polygon_soup(points, facets);
    }
    if (ret) {
        CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points,
            facets, mesh);
        ret = mesh.is_valid() && !mesh.is_empty();
    }
    if (ret) {
        if (avgTriangleArea > 0.0) {
            double totalArea = CGAL::Polygon_mesh_processing::area(mesh,
                CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh.points()).
                geom_traits(K()));
            double avgArea = totalArea / mesh.number_of_faces();
            double densityControlFactor = avgArea / avgTriangleArea;

            std::vector<Surface_mesh::Face_index> new_facets;
            std::vector<Surface_mesh::Vertex_index> new_vertex;

            CGAL::Polygon_mesh_processing::refine(mesh, faces(mesh),
                std::back_inserter(new_facets),
                std::back_inserter(new_vertex),
                CGAL::Polygon_mesh_processing::parameters::density_control_factor(densityControlFactor).
                vertex_point_map(mesh.points()));
        }
        std::cout << "Loaded " << filename << " with " <<
            mesh.number_of_faces() << " faces and " <<
            mesh.number_of_vertices() << " vertices." << std::endl;

        mesh_ = mesh;
    }

    ifs.close();
    return false;
}

bool Surface::isLoaded() {
    return mesh_.is_valid() && !mesh_.is_empty();
}

bool Surface::computeIntersection(Ray r) {
    Point source = r.source();

    std::list<Ray_intersection> intersections;
    tree_.all_intersections(r, std::back_inserter(intersections));

    bool found = false;
    double nearestDistance = 0.0;
    Point *nearestPoint;
    face_descriptor nearestFaceId;
    Point *p = NULL;
    for (Ray_intersection intersection : intersections) {
        if (intersection && (p = boost::get<Point>(&(intersection->first)))) {
            double distance = CGAL::squared_distance(source, *p);
            if (!found || distance < nearestDistance) {
                found = true;
                nearestPoint = p;
                nearestDistance = distance;
                nearestFaceId = intersection->second;
            }
        }
    }

    if (found) {
        std::cout << "Found intersection at point (" << nearestPoint->x() <<
            ", " << nearestPoint->y() << ", " << nearestPoint->z() << ")" <<
            std::endl;
        std::cout << "Face id: " << nearestFaceId << std::endl;
        Vector n = faceNormals_[nearestFaceId];
        std::cout << "Face normal: (" << n.x() << ", " << n.y() <<
            ", " << n.z() << ")" << std::endl << std::endl;
    }

    return found;
}

