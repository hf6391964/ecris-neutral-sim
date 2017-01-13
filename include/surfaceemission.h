#pragma once

#include "neutralsource.h"
#include "surface.h"
#include "particle.h"
#include "element_data.h"

class SurfaceEmission : public NeutralSource {
    Surface* pSurface_;
    double emissionRate_;
    Element element_;
    const ElementData *elementData_;

    public:
        SurfaceEmission(Surface* pSurface, double emissionRate_,
            Element element_);

        Particle generateParticle(Rng& rng) const;
};

