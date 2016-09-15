#ifndef SURFACE_H
#define SURFACE_H

#include <string>
#include <algorithm>

#include "cgal_and_typedefs.h"

class Surface;

struct IntersectionPoint {
    Surface* pSurface = NULL;
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
    double pumpingFactor_;
    double temperature_;
    double emissionRate_;
    unsigned long pumpedParticles_ = 0;

    public:
        Surface() {};
        Surface(std::string filename, double pumpingFactor, double temperature,
            double emissionRate, double avgTriangleArea = -1.0)
            : pumpingFactor_(pumpingFactor), temperature_(temperature),
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

        template<typename OutputIterator>
        void computeIntersections(const Ray& r, OutputIterator out) const {
            tree_.all_intersections(r, out);
        }

        Direction generateCosineLawDirection(face_descriptor fd, Rng& rng) const;

        bool isLoaded() const {
            return mesh_.is_valid() && !mesh_.is_empty();
        }

        bool isEmissive() const {
            return emissionRate_ > 0.0;
        }

        double getEmissionRate() const {
            return isEmissive() ? emissionRate_ : 0.0;
        }

        double getTemperature() const {
            return temperature_;
        }

        unsigned long getPumpedParticles() const {
            return pumpedParticles_;
        }

        Bbox bbox() const {
            return CGAL::Polygon_mesh_processing::bbox_3(mesh_,
                CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh_.points()));
        }

        bool checkIfPumped(Rng& rng) const {
            return uni01(rng) < pumpingFactor_;
        }

        void addPumpedParticle() {
            pumpedParticles_ += 1;
        }
};

#endif
