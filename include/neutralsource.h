#pragma once

#include "particle.h"

class NeutralSource {
    std::string label_;
    public:
        NeutralSource(std::string label) : label_(label) {}
        std::string getLabel() const { return label_; }
        virtual Particle generateParticle(Rng& rng, double time) const = 0;
        // Return emission rate in particles/s
        virtual double getEmissionRate() const = 0;
};

