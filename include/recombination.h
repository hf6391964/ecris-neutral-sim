#pragma once

#include <string>

#include "neutralizationchannel.h"
#include "particlepopulation.h"
#include "particle.h"

class Recombination : public NeutralizationChannel {
    std::shared_ptr<ParticlePopulation> population_;

    public:
        Recombination(std::shared_ptr<ParticlePopulation> population,
            double reactionRate);

        Particle sampleNeutralProduct(Rng &rng,
            const Particle &sourceParticle, double decayTime) const;
};

