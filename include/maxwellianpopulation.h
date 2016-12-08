#pragma once

#include "cgal_and_typedefs.h"
#include "particlepopulation.h"

class MaxwellianPopulation : ParticlePopulation {
    double temperature_eV_;

    public:
        MaxwellianPopulation(double particleMass_eV, int chargeState, double
            temperature_eV, DensityDistribution distr)
            : ParticlePopulation(particleMass_eV, chargeState, distr),
              temperature_eV_(temperature_eV) {}

        Vector getRandomParticleVelocity(Rng &rng) const;
        double getRelativeSpeed(double particleSpeed, monte_state *ms,
            gsl_rng *rng) const;
};
