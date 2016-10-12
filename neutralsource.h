#ifndef NEUTRALSOURCE_H
#define NEUTRALSOURCE_H

#include "particle.h"

class NeutralSource {
    public:
        virtual void generateParticle(Particle& particle, Rng& rng) const = 0;
};

#endif

