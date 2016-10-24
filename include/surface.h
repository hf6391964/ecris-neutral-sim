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
    std::atomic_ulong pumpedParticles_;

    public:
        Surface() {
            pumpedParticles_ = 0;
        }
        Surface(std::string filename, double pumpingFactor, double temperature,
            bool flipNormals = false, double avgTriangleArea = -1.0)
            : pumpingFactor_(pumpingFactor), temperature_(temperature) {
            pumpedParticles_ = 0;
            loadFromSTL(filename, avgTriangleArea);
            buildAABBTree();
            computeAreaCDF();
            computeFaceNormals();
            computeFaceRotations(flipNormals);
        }

        bool loadFromSTL(std::string filename, double avgTriangleArea = -1.0);
        void buildAABBTree();
        void computeAreaCDF();
        void computeFaceNormals();
        void computeFaceRotations(bool flipNormals);

        std::tuple<Point, face_descriptor> getRandomPoint(Rng& rng) const;

        void computeFirstIntersection(const Ray& r,
            Ray_intersection& isect) const {
            isect = tree_.first_intersection(r);
        }

        void computeAllIntersections(const Ray& r,
            std::back_insert_iterator<std::vector<Ray_intersection>> it) const {
            tree_.all_intersections(r, it);
        }

        Direction generateCosineLawDirection(face_descriptor fd, Rng& rng) const;

        bool isLoaded() const {
            return mesh_.is_valid() && !mesh_.is_empty();
        }

        double getTemperature() const {
            return temperature_;
        }

        unsigned long getPumpedParticles() const {
            return pumpedParticles_;
        }

        bool checkIfPumped(Rng& rng) const {
            return uni01(rng) < pumpingFactor_;
        }

        void addPumpedParticle() {
            pumpedParticles_ += 1;
        }

        Bbox bbox() const;
};

#endif