#include "maxwellianpopulation.h"

Vector MaxwellianPopulation::getRandomParticleVelocity(Rng &rng) const {
    return Util::getMBVelocity(rng, temperature_eV_, particleMass_eV_);
}

double MaxwellianPopulation::getRelativeSpeed(double particleSpeed,
    monte_state *ms, gsl_rng *rng) const {
    return Util::calculateMBRelativeSpeed(particleSpeed, temperature_eV_,
        particleMass_eV_, rng, ms);
}

double MaxwellianPopulation::calculateRateCoefficient(double particleSpeed,
    monte_state *ms, gsl_rng *rng,
    double (*fCrossSection)(double, void *), void *fArgs) const {
    return Util::calculateMBRelativeRateCoeff(particleSpeed, temperature_eV_,
        particleMass_eV_, rng, ms, fCrossSection, fArgs);
}

