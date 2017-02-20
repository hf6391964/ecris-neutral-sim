#pragma once

#include "neutralsource.h"
#include "surface.h"
#include "particle.h"
#include "element_data.h"

class SurfaceEmission : public NeutralSource {
    std::shared_ptr<Surface> pSurface_;
    double emissionRate_;
    Element element_;
    const ElementData *elementData_;

    public:
        SurfaceEmission(std::shared_ptr<Surface> pSurface, double emissionRate_,
            Element element_, std::string label);

        Particle generateParticle(Rng& rng, double time) const;
};

