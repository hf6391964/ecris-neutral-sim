#pragma once

#include "cgal_and_typedefs.h"
#include "densitydistribution.h"

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

        virtual Vector getRandomParticleVelocity(Rng &rng) const = 0;
        Point getRandomParticlePosition(Rng &rng) const;

        double getDensityAt(const Point &p) const;

        virtual double getRelativeSpeed(double particleSpeed, monte_state *ms,
            gsl_rng *rng) const = 0;

        void setCoordinateTransformation(const Aff_transformation
            &transformation);
        void removeCoordinateTransformation();
};

