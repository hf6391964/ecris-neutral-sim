#include "particlepopulation.h"

double ParticlePopulation::getDensityAt(const Point &p) const {
    return densityDistribution_->getValueAt(p);
}

double ParticlePopulation::getRandomParticleSpeed(Rng &rng) const {
    return std::sqrt(getRandomParticleVelocity(rng).squared_length());
}

Point ParticlePopulation::getRandomParticlePosition(Rng &rng) const {
    return densityDistribution_->getRandomPosition(rng);
}

