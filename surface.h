#ifndef SURFACE_H
#define SURFACE_H

#include <string>
#include <algorithm>

#include "cgal_and_typedefs.h"

class Surface {
    Surface_mesh mesh_;
    Surface_mesh::Property_map<face_descriptor, Vector> faceNormals_;
    double stickingFactor_;
    double sojournTime_;
    double temperature_;
    double emissionRate_;
    Tree tree_;
    Bbox bbox_;
    std::vector<double> areaCDF_;

    public:
        bool loadFromSTL(std::string filename, double avgTriangleArea = -1.0);
        bool isLoaded();
        void buildAABBTree();
        void computeAreaCDF();
        void computeFaceNormals();
        void computeFaceMidpoints();
        Point getRandomPoint(Rng &rng);
        bool computeIntersection(Ray r);

        Surface() {};
        Surface(std::string filename, double avgTriangleArea = -1.0) {
            loadFromSTL(filename, avgTriangleArea);
            buildAABBTree();
            computeAreaCDF();
            computeFaceNormals();
        }
};

#endif
