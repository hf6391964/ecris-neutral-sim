#pragma once

#include <string>

#include "neutralizationchannel.h"
#include "surfacecollection.h"
#include "particlepopulation.h"
#include "particle.h"

class WallNeutralization : NeutralizationChannel {
    std::vector<IntersectionPoint> wallPoints_, end1Points_, end2Points_;
    std::vector<std::vector<IntersectionPoint> *> pointSets_;
    std::vector<double> cumulativeNormalizedEndpointCount_;
    std::shared_ptr<ParticlePopulation> population_;
    double confinementTime_;

    public:
        WallNeutralization(
            std::shared_ptr<ParticlePopulation> population_,
            double confinementTime,
            const std::string &wallFilename,
            const std::string &end1Filename,
            const std::string &end2Filename,
            const SurfaceCollection &surfaces);

        Particle sampleNeutralProduct(Rng &rng,
            const Particle &sourceParticle, double decayTime) const;
};

