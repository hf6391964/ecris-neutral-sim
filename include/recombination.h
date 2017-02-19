#pragma once

#include <string>

#include "neutralizationchannel.h"
#include "particlepopulation.h"
#include "particle.h"

class Recombination : public NeutralizationChannel {
    const ParticlePopulation &population_;

    public:
        Recombination(const ParticlePopulation &population,
            double reactionRate);

        Particle sampleNeutralProduct(Rng &rng,
            const Particle &sourceParticle, double decayTime) const;
};

