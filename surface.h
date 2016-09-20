#ifndef SURFACE_H
#define SURFACE_H

#include <string>
#include <algorithm>
#include <atomic>

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
    std::atomic_ulong pumpedParticles_;

    public:
        Surface() {
            pumpedParticles_ = 0;
        }
        Surface(std::string filename, double pumpingFactor, double temperature,
            double emissionRate, double avgTriangleArea = -1.0)
            : pumpingFactor_(pumpingFactor), temperature_(temperature),
              emissionRate_(emissionRate) {
            pumpedParticles_ = 0;
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

        void computeFirstIntersection(const Ray& r,
            Ray_intersection& isect) const {
            isect = tree_.first_intersection(r);
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
