#pragma once

#include "cgal_and_typedefs.h"
#include "spatialdistribution.h"

class ParticlePopulation {
    protected:
        double particleMass_eV_;
        int chargeState_;
        DensityDistribution densityDistribution_;

    public:
        ParticlePopulation(double particleMass_eV, int chargeState,
            DensityDistribution densityDistribution)
            : particleMass_eV_(particleMass_eV), chargeState_(chargeState),
              densityDistribution_(densityDistribution) {}
        virtual ~ParticlePopulation() {};

        virtual Vector getRandomParticleVelocity(Rng &rng) const = 0;
        Point getRandomParticlePosition(Rng &rng) const;

        double getParticleMass_eV() const { return particleMass_eV_; }
        int getChargeState() const { return chargeState_; }

        double getDensityAt(const Point &p) const;

        virtual double getRelativeSpeed(double particleSpeed, monte_state *ms,
            gsl_rng *rng) const = 0;

        void setCoordinateTransformation(const Aff_transformation
            &transformation);
        void removeCoordinateTransformation();
};

