#pragma once

#include "cgal_and_typedefs.h"
#include "spatialdistribution.h"
#include "constants.h"
#include "element_data.h"

class ParticlePopulation {
    protected:
        Element element_ = ELEMENT_NONE;
        double particleMass_eV_;
        int chargeState_;
        DensityDistribution densityDistribution_;

    public:
        ParticlePopulation(Element element, int chargeState,
            DensityDistribution densityDistribution)
            : element_(element), chargeState_(chargeState),
              densityDistribution_(densityDistribution) {
            particleMass_eV_ = ELEMENT_DATA.at(element_)->mass *
                ATOMIC_MASS_TO_EV;
        }

        ParticlePopulation(double particleMass_eV, int chargeState,
            DensityDistribution densityDistribution)
            : particleMass_eV_(particleMass_eV), chargeState_(chargeState),
              densityDistribution_(densityDistribution) {}

        virtual ~ParticlePopulation() {};

        Element getElement() const {
            return element_;
        }

        virtual Vector getRandomParticleVelocity(Rng &rng) const = 0;
        Point getRandomParticlePosition(Rng &rng) const;

        double getParticleMass_eV() const { return particleMass_eV_; }
        int getChargeState() const { return chargeState_; }

        double getDensityAt(const Point &p) const;

        virtual double getRelativeSpeed(double particleSpeed, monte_state *ms,
            gsl_rng *rng) const = 0;

        virtual double calculateRateCoefficient(double particleSpeed,
            monte_state *ms, gsl_rng *rng,
            double (*fCrossSection)(double, void *), void *fArgs) const = 0;

        void setCoordinateTransformation(const Aff_transformation
            &transformation);
        void removeCoordinateTransformation();
};

