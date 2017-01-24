#pragma once

#include "cgal_and_typedefs.h"
#include "particlepopulation.h"

class MaxwellianPopulation : public ParticlePopulation {
    double temperature_eV_;

    public:
        MaxwellianPopulation(Element element, int chargeState,
            double temperature_eV, DensityDistribution distr)
            : ParticlePopulation(element, chargeState, distr),
              temperature_eV_(temperature_eV) {}
        MaxwellianPopulation(double particleMass_eV, int chargeState,
            double temperature_eV, DensityDistribution distr)
            : ParticlePopulation(particleMass_eV, chargeState, distr),
              temperature_eV_(temperature_eV) {}

        Vector getRandomParticleVelocity(Rng &rng) const;
        double getRelativeSpeed(double particleSpeed, monte_state *ms,
            gsl_rng *rng) const;
        double calculateRateCoefficient(double particleSpeed,
            monte_state *ms, gsl_rng *rng,
            double (*fCrossSection)(double, void *), void *fArgs) const;
        double getTemperature() const;
};

