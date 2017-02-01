#pragma once

#include "particle.h"

class NeutralizationChannel {
    double decayTime_ = 0.0;
    std::string label_;
    std::atomic_ulong reactionCounter_;

    public:
        NeutralizationChannel(double decayTime, const std::string& label)
            : decayTime_(decayTime), label_(label), reactionCounter_(0) {}

        double timeToReaction(Rng &rng) const;

        std::string getLabel() const;

        virtual Particle sampleNeutralProduct(Rng &rng,
            const Particle &sourceParticle, double decayTime) const = 0;

        void incrementReactionCounter();

        unsigned long getReactionCount() const;
};

