#pragma once

#include "neutralsource.h"
#include "surface.h"
#include "particle.h"

class SurfaceEmission : public NeutralSource {
    Surface* pSurface_;
    double emissionRate_;
    public:
        SurfaceEmission(Surface* pSurface, double emissionRate_);

        void generateParticle(Particle& particle, Rng& rng) const;
};

