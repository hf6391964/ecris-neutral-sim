#include <algorithm>

#include "wallneutralization.h"
#include "electronmodel.h"


WallNeutralization::WallNeutralization(
    std::shared_ptr<ParticlePopulation> population,
    double confinementTime,
    const std::string &wallFilename,
    const std::string &end1Filename,
    const std::string &end2Filename,
    const SurfaceCollection &surfaces,
    double spatialResolution, double angularResolution)
    : NeutralizationChannel(confinementTime, "wall neutralization"),
      population_(population),
      spatialResolution_(spatialResolution), angularResolution_(angularResolution),
      surfaces_(surfaces) {

    std::vector<Point> wallEndpoints =
        ElectronModel::parseElectronEndpoints(wallFilename);
    std::vector<Point> end1Endpoints =
        ElectronModel::parseElectronEndpoints(end1Filename);
    std::vector<Point> end2Endpoints =
        ElectronModel::parseElectronEndpoints(end2Filename);

    for (Point p : wallEndpoints) {
        double theta = std::min(M_PI, std::max(-M_PI, std::atan2(p.y(), p.x())));
        double thetaBinned = angularResolution_ * Util::fastFloor(theta / angularResolution_);
        double zBinned = spatialResolution_ * Util::fastFloor(p.z() / spatialResolution_);
        double r = std::sqrt(p.x()*p.x() + p.y()*p.y());
        Endpoint ep = { r, thetaBinned, zBinned };
        wallPoints_.push_back(ep);
    }
    pointSets_.push_back(std::move(wallPoints_));

    for (Point p : end1Endpoints) {
        double xBinned = spatialResolution_ * Util::fastFloor(p.x() / spatialResolution);
        double yBinned = spatialResolution_ * Util::fastFloor(p.y() / spatialResolution);
        Endpoint ep = { xBinned, yBinned, p.z() };
        end1Points_.push_back(ep);
    }
    pointSets_.push_back(std::move(end1Points_));

    for (Point p : end2Endpoints) {
        double xBinned = spatialResolution_ * Util::fastFloor(p.x() / spatialResolution);
        double yBinned = spatialResolution_ * Util::fastFloor(p.y() / spatialResolution);
        Endpoint ep = { xBinned, yBinned, p.z() };
        end2Points_.push_back(ep);
    }
    pointSets_.push_back(std::move(end2Points_));

    unsigned long count = 0;
    std::vector<unsigned long> countVector;
    for (EndpointSetVector::const_iterator itPoints = pointSets_.begin();
        itPoints != pointSets_.end(); ++itPoints) {
        count += (*itPoints).size();
        countVector.push_back(count);
    }

    for (unsigned long c : countVector) {
        cumulativeNormalizedEndpointCount_.push_back((double)c / count);
    }
}

Particle WallNeutralization::sampleNeutralProduct(Rng &rng,
    const Particle &sourceParticle, double decayTime) const {
    Particle result(sourceParticle.getElement(), 0.0);
    // Time will be set later
    double resultTime = sourceParticle.getTime() + decayTime;

    // 1. Sample the position of the wall intersection
    size_t i = std::lower_bound(cumulativeNormalizedEndpointCount_.begin(),
        cumulativeNormalizedEndpointCount_.end(), uni01(rng)) -
        cumulativeNormalizedEndpointCount_.begin();
    size_t j = uni01(rng) * pointSets_.at(i).size();
    Endpoint ep = pointSets_.at(i).at(j);
    Ray ray;
    if (i == 0) {
        // Endpoint is on the radial wall
        double r = ep.coord1,
               theta = ep.coord2 + uni01(rng) * angularResolution_,
               z = ep.coord3 + uni01(rng) * spatialResolution_;

        ray = Ray(Point(0.0, 0.0, z), Point(
            r * std::cos(theta),
            r * std::sin(theta),
            z
        ));
    } else {
        // Endpoint is on either of the ends
        ray = Ray(Point(0.0, 0.0, 0.0), Point(
            ep.coord1 + uni01(rng) * spatialResolution_,
            ep.coord2 + uni01(rng) * spatialResolution_,
            ep.coord3
        ));
    }

    // 2. Sample particle speed from the plasma 1+ distribution
    double speed = population_->getRandomParticleSpeed(rng);

    // Ray from the particle start point to the wall endpoint
    IntersectionPoint ip;
    if (!surfaces_.findClosestIntersection(ray, ip)) {
        std::cout << i << std::endl;
        std::cout << "ray = " << ray << "\n";
        throw std::out_of_range("Could not find wall intersection point ");
    }

    result.setVelocity(speed, ray.direction());
    result.setNextIntersection(ip);

    // 3. Thermal accommodation (handled by Surface)
    result.goToIntersection(rng);
    result.setTime(resultTime);

    // 4. Assume particle is neutralized on the wall and return it as the
    // product

    return result;
}

