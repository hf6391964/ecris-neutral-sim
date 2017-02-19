#pragma once

#include <string>

#include "neutralizationchannel.h"
#include "surfacecollection.h"
#include "particlepopulation.h"
#include "particle.h"

struct Endpoint {
    double coord1, coord2, coord3;
    // Cylinder points: r, theta, z
    // Cartesian points: x, y, z
};
typedef std::vector<std::vector<Endpoint>> EndpointSetVector;

class WallNeutralization : public NeutralizationChannel {
    std::vector<Endpoint> wallPoints_, end1Points_, end2Points_;
    EndpointSetVector pointSets_;
    std::vector<double> cumulativeNormalizedEndpointCount_;
    const ParticlePopulation &population_;
    double confinementTime_;
    double spatialResolution_;
    double angularResolution_;
    const SurfaceCollection &surfaces_;

    public:
        WallNeutralization(
            const ParticlePopulation &population_,
            double confinementTime,
            const std::string &wallFilename,
            const std::string &end1Filename,
            const std::string &end2Filename,
            const SurfaceCollection &surfaces,
            double spatialResolution = 0.001,
            double angularResolution = 0.018);

        Particle sampleNeutralProduct(Rng &rng,
            const Particle &sourceParticle, double decayTime) const;
};

