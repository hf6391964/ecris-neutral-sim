#include "wallneutralization.h"
#include "electronmodel.h"


WallNeutralization::WallNeutralization(double confinementTime,
    const std::string &wallFilename,
    const std::string &end1Filename,
    const std::string &end2Filename,
    const SurfaceCollection &surfaces)
    : NeutralizationChannel(confinementTime) {
    std::vector<Point> wallEndpoints =
        ElectronModel::parseElectronEndpoints(wallFilename);
    std::vector<Point> end1Endpoints =
        ElectronModel::parseElectronEndpoints(end1Filename);
    std::vector<Point> end2Endpoints =
        ElectronModel::parseElectronEndpoints(end2Filename);

    // Find intersections with radial wall
    for (Point p : wallEndpoints) {
        // Construct a ray originating from the z axis to the wall intersection
        // point
        Ray r(Point(0.0, 0.0, p.z()), p);
        IntersectionPoint ip;
        if (surfaces.findClosestIntersection(r, ip)) {
            wallPoints_.push_back(ip);
        }
    }
    // Find intersections with end 1
    for (Point p : end1Endpoints) {
        // Construct a ray starting from origin, going through the point p
        Ray r(Point(0.0, 0.0, 0.0), p);
        IntersectionPoint ip;
        if (surfaces.findClosestIntersection(r, ip)) {
            end1Points_.push_back(ip);
        }
    }
    // Find intersections with end 2
    for (Point p : end2Endpoints) {
        // Construct a ray starting from origin, going through the point p
        Ray r(Point(0.0, 0.0, 0.0), p);
        IntersectionPoint ip;
        if (surfaces.findClosestIntersection(r, ip)) {
            end2Points_.push_back(ip);
        }
    }
}

