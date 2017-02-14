#pragma once

#include <string>
#include <algorithm>
#include <atomic>
#include <unordered_map>

#include "cgal_and_typedefs.h"
#include "element_data.h"

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
    std::atomic_ulong collisionCounter_;
    std::unordered_map<Element, double, element_hash> accommodationCoefficients_;
    const double DEFAULT_ACCOMMODATION_COEFFICIENT = 1.0;
    std::string label_;

    public:
        Surface(std::string filename, double pumpingFactor, double temperature,
            std::string label, bool flipNormals = false,
            double scale = 1.0,
            double edgeLength = 0.0);

        void setAccommodationCoefficients(
            std::unordered_map<Element, double, element_hash> coefficients);

        double getAccommodationCoefficient(Element element) const;

        bool loadFromSTL(std::string filename, double scale = 1.0,
            double edgeLength = 0.0);
        void buildAABBTree();
        void computeAreaCDF();
        void computeFaceNormals();
        void computeFaceRotations(bool flipNormals);

        std::tuple<Point, face_descriptor> getRandomPoint(Rng& rng) const;

        void computeFirstIntersection(const Ray& r,
            Ray_intersection& isect) const;

        void computeAllIntersections(const Ray& r,
            std::back_insert_iterator<std::vector<Ray_intersection>> it) const;

        Direction generateCosineLawDirection(face_descriptor fd, Rng& rng) const;

        bool isLoaded() const;

        double getTemperature() const;

        unsigned long getPumpedParticles() const;

        bool checkIfPumped(Rng& rng) const;

        void addPumpedParticle();

        void incrementCollisionCounter();

        unsigned long getCollisionCount() const;

        Bbox bbox() const;

        std::string getLabel() const;
};

