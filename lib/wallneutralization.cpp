#include <algorithm>

#include "wallneutralization.h"
#include "electronmodel.h"


WallNeutralization::WallNeutralization(
    std::shared_ptr<ParticlePopulation> population,
    double confinementTime,
    const std::string &wallFilename,
    const std::string &end1Filename,
    const std::string &end2Filename,
    const SurfaceCollection &surfaces)
    : NeutralizationChannel(confinementTime), population_(population) {
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
    pointSets_.push_back(&wallPoints_);

    // Find intersections with end 1
    for (Point p : end1Endpoints) {
        // Construct a ray starting from origin, going through the point p
        Ray r(Point(0.0, 0.0, 0.0), p);
        IntersectionPoint ip;
        if (surfaces.findClosestIntersection(r, ip)) {
            end1Points_.push_back(ip);
        }
    }
    pointSets_.push_back(&end1Points_);

    // Find intersections with end 2
    for (Point p : end2Endpoints) {
        // Construct a ray starting from origin, going through the point p
        Ray r(Point(0.0, 0.0, 0.0), p);
        IntersectionPoint ip;
        if (surfaces.findClosestIntersection(r, ip)) {
            end2Points_.push_back(ip);
        }
    }
    pointSets_.push_back(&end2Points_);

    unsigned long count = 0;
    std::vector<unsigned long> countVector(pointSets_.size());
    for (std::vector<IntersectionPoint> *points : pointSets_) {
        count += points->size();
        countVector.push_back(count);
    }

    for (unsigned long c : countVector) {
        cumulativeNormalizedEndpointCount_.push_back((double)c / count);
    }
}

Particle WallNeutralization::sampleNeutralProduct(Rng &rng,
    const Particle &sourceParticle, double decayTime) const {
    Particle result(sourceParticle.getElement(),
        sourceParticle.getTime() + decayTime);

    // 1. Sample particle position from the wall distribution population
    double x = uni01(rng);
    size_t i = std::lower_bound(cumulativeNormalizedEndpointCount_.begin(),
        cumulativeNormalizedEndpointCount_.end(), x) -
        cumulativeNormalizedEndpointCount_.begin();
    std::vector<IntersectionPoint> *points = pointSets_[i];
    size_t j = uni01(rng) * points->size();
    IntersectionPoint ip = points->at(j);
    result.setNextIntersection(ip);
    result.setPosition(ip.point);

    // 2. Sample particle velocity from the plasma 1+ distribution
    Vector vel = population_->getRandomParticleVelocity(rng);
    result.setVelocity(vel);

    // 3. Thermal accommodation (handled by Surface)
    result.goToIntersection(rng);

    // 4. Assume particle is neutralized on the wall and return it as the
    // product

    return result;
}

