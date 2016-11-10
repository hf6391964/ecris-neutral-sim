#include "surface.h"
#include "STL_reader.h"

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

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
            double totalArea = 0.0;
            for (face_descriptor fd : mesh.faces()) {
                totalArea += CGAL::Polygon_mesh_processing::face_area(fd, mesh,
                    CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh.points()));
            }
            double avgArea = totalArea / mesh.number_of_faces();
            double densityControlFactor = avgArea / avgTriangleArea;

            std::vector<Surface_mesh::Face_index> new_facets;
            std::vector<Surface_mesh::Vertex_index> new_vertex;

            CGAL::Polygon_mesh_processing::refine(mesh, faces(mesh),
                std::back_inserter(new_facets),
                std::back_inserter(new_vertex),
                CGAL::Polygon_mesh_processing::parameters::density_control_factor(densityControlFactor));
        }
        std::cout << "Loaded " << filename << " with " <<
            mesh.number_of_faces() << " faces and " <<
            mesh.number_of_vertices() << " vertices." << std::endl;

        mesh_ = mesh;
    }

    ifs.close();
    return false;
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
    std::cout << "Found index: " << (int)fd << std::endl <<
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

