#pragma once

#include "particle.h"

class NeutralizationChannel {
    double decayTime_ = 0.0;

    public:
        NeutralizationChannel(double decayTime) : decayTime_(decayTime) {}

        double timeToReaction(Rng &rng) const;

        virtual Particle sampleNeutralProduct(Rng &rng,
            const Particle &sourceParticle, double decayTime) const = 0;
};

