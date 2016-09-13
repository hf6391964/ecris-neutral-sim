#ifndef SURFACE_H
#define SURFACE_H

#include <string>
#include <algorithm>

#include "cgal_and_typedefs.h"

class Surface;

struct IntersectionPoint {
    Surface& surface;
    face_descriptor faceId;
    Point point;
};

class Surface {
    Tree tree_;
    Bbox bbox_;
    std::vector<double> areaCDF_;
    Surface_mesh mesh_;
    Surface_mesh::Property_map<face_descriptor, Vector> faceNormals_;
    Surface_mesh::Property_map<face_descriptor, K::Aff_transformation_3>
        faceRotations_;
    double stickingFactor_;
    double temperature_;
    double emissionRate_;

    public:
        Surface() {};
        Surface(std::string filename, double stickingFactor, double temperature,
            double emissionRate, double avgTriangleArea = -1.0)
            : stickingFactor_(stickingFactor), temperature_(temperature),
              emissionRate_(emissionRate) {
            loadFromSTL(filename, avgTriangleArea);
            buildAABBTree();
            computeAreaCDF();
            computeFaceNormals();
            computeFaceRotations();
        }

        bool loadFromSTL(std::string filename, double avgTriangleArea = -1.0);
        void buildAABBTree();
        void computeAreaCDF();
        void computeFaceNormals();
        void computeFaceRotations();

        std::tuple<Point, face_descriptor> getRandomPoint(Rng& rng) const;
        bool computeIntersection(const Ray& r) const;
        Direction generateCosineLawDirection(face_descriptor fd, Rng& rng) const;

        bool isLoaded() const {
            return mesh_.is_valid() && !mesh_.is_empty();
        };

        bool isEmissive() const {
            return emissionRate_ > 0.0;
        }

        double getEmissionRate() const {
            return isEmissive() ? emissionRate_ : 0.0;
        }

        Bbox bbox() const {
            return CGAL::Polygon_mesh_processing::bbox_3(mesh_,
                CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh_.points()));
        }
};

#endif
