#pragma once

#include "particle.h"

class NeutralSource {
    public:
        virtual Particle generateParticle(Rng& rng) const = 0;
};

