#pragma once

#include "particle.h"

class NeutralizationChannel {
    double decayTime_ = 0.0;
    std::string label_;

    public:
        NeutralizationChannel(double decayTime, const std::string& label)
            : decayTime_(decayTime), label_(label) {}

        double timeToReaction(Rng &rng) const;

        std::string getLabel() const;

        virtual Particle sampleNeutralProduct(Rng &rng,
            const Particle &sourceParticle, double decayTime) const = 0;
};

