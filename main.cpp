#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/STL_reader.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;

bool stlToPolyhedron(std::string filename, Polyhedron &mesh) {
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

    ifs.close();
    return ret;
}

int main() {
    Polyhedron cap1;

    bool success = stlToPolyhedron("models/cylinder/cylinder_cap1.stl", cap1);

    if (success) {
        std::cout << "Great!" << std::endl;
    }

    return 0;
}

