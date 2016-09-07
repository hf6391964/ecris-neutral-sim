#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <random>
#include <cmath>

#include "cgal_and_typedefs.h"

bool stlToPolyhedron(std::string filename, Polyhedron &mesh,
    double density_control_factor=-1.0) {
    std::ifstream ifs;

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
        ret = mesh.is_valid() && !mesh.empty();
    }
    if (ret) {
        if (density_control_factor > 0.0) {
            std::vector<Polyhedron::Facet_handle> new_facets;
            std::vector<Polyhedron::Vertex_handle> new_vertex;
            CGAL::Polygon_mesh_processing::refine(mesh, faces(mesh),
                std::back_inserter(new_facets),
                std::back_inserter(new_vertex),
                CGAL::Polygon_mesh_processing::parameters::density_control_factor(density_control_factor));
        }
        std::cout << "Loaded " << filename << " with " <<
            mesh.size_of_facets() << " facets and " <<
            mesh.size_of_vertices() << " vertices." << std::endl;
    }

    ifs.close();
    return ret;
}

static std::uniform_real_distribution<double> theta_dist(0, M_PI);
static std::uniform_real_distribution<double> phi_dist(0, 2.0*M_PI);
static std::default_random_engine re;
Ray randomRay() {
    double theta = theta_dist(re);
    double phi = phi_dist(re);

    return Ray(Point(0, 0, 0), Direction(std::sin(theta)*std::cos(phi),
                                         std::sin(theta)*std::sin(phi),
                                         std::cos(theta)));
}

int main() {
    Polyhedron mesh_cap1;
    Polyhedron mesh_cap2;
    Polyhedron mesh_walls;
    Polyhedron mesh_all;

    stlToPolyhedron("models/cylinder/cylinder_cap1.stl", mesh_cap1);
    stlToPolyhedron("models/cylinder/cylinder_cap2.stl", mesh_cap2);
    stlToPolyhedron("models/cylinder/cylinder_walls.stl", mesh_walls, 20.0);
    stlToPolyhedron("models/cylinder/cylinder_all.stl", mesh_all);
    std::cout << std::endl;

    const int N_RAYS = 1000000;
    Tree tree(faces(mesh_walls).first, faces(mesh_walls).second, mesh_walls);

    for (int i = 0; i < N_RAYS; i++) {
        Ray r = randomRay();
        Direction d = r.direction();
        /* std::cout << "Ray " << i << ": (" << d.dx() << ", " << d.dy() << */
        /*     ", " << d.dz() << ")" << std::endl; */
        /* std::cout << "do_intersect: " << tree.do_intersect(r) << std::endl; */
        Ray_intersection isect = tree.any_intersection(r);
        Point *p;
        if (isect && (p = boost::get<Point>(&(isect->first)))) {
            /* std::cout << "Intersection detected at " << "(" << p->x() << */
            /*     ", " << p->y() << ", " << p->z() << ")" << std::endl; */
            /* std::cout << "xy radius: " << */
                /* std::sqrt(p->x()*p->x() + p->y()*p->y()) << std::endl; */
        } else {
            /* std::cout << "No intersection" << std::endl; */
        }
        /* std::cout << std::endl; */
    }

    return 0;
}

