#include "surface.h"

Surface Surface::loadFromSTL(std::string filename, double avgTriangleArea) {
    std::ifstream ifs;
    Surface surf;
    Polyhedron mesh;

    ifs.open(filename, std::ifstream::in);

    if (!ifs.is_open()) return surf;

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
        if (avgTriangleArea > 0.0) {
            double totalArea = CGAL::Polygon_mesh_processing::area(mesh);
            double avgArea = totalArea / mesh.size_of_facets();
            double densityControlFactor = avgArea / avgTriangleArea;

            std::vector<Polyhedron::Facet_handle> new_facets;
            std::vector<Polyhedron::Vertex_handle> new_vertex;

            CGAL::Polygon_mesh_processing::refine(mesh, faces(mesh),
                std::back_inserter(new_facets),
                std::back_inserter(new_vertex),
                CGAL::Polygon_mesh_processing::parameters::density_control_factor(densityControlFactor));
        }
        std::cout << "Loaded " << filename << " with " <<
            mesh.size_of_facets() << " facets and " <<
            mesh.size_of_vertices() << " vertices." << std::endl;

        surf.mesh_ = mesh;
    }

    ifs.close();
    return surf;
}

bool Surface::isLoaded() {
    return mesh_.is_valid() && !mesh_.is_empty();
}
