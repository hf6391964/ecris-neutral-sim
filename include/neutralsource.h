#pragma once

#include "particle.h"

class NeutralSource {
    public:
        virtual Particle generateParticle(Rng& rng, double time) const = 0;
};

