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
        // The emission rate is in [ Pa m^3 / s ]
        SurfaceEmission(std::shared_ptr<Surface> pSurface, double emissionRate,
            Element element_, std::string label);

        Particle generateParticle(Rng& rng, double time) const;

        double getEmissionRate() const {
            return emissionRate_ /
                (pSurface_->getTemperature() * EV_TO_JOULE);
        }
};

