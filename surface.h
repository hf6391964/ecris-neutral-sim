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
        bool loadFromSTL(std::string filename, double avgTriangleArea = -1.0);
        bool isLoaded() const;
        void buildAABBTree();
        void computeAreaCDF();
        void computeFaceNormals();
        void computeFaceRotations();
        Point getRandomPoint(Rng& rng) const;
        bool computeIntersection(const Ray& r) const;
        Direction generateCosineLawDirection(face_descriptor fd, Rng& rng) const;

        Surface() {};
        Surface(std::string filename, double avgTriangleArea = -1.0) {
            loadFromSTL(filename, avgTriangleArea);
            buildAABBTree();
            computeAreaCDF();
            computeFaceNormals();
            computeFaceRotations();
        }
};

#endif
