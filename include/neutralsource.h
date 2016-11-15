#pragma once

#include "particle.h"

class NeutralSource {
    public:
        virtual void generateParticle(Particle& particle, Rng& rng) const = 0;
};

