#include "particlepopulation.h"

double ParticlePopulation::getDensityAt(const Point &p) const {
    return densityDistribution_.getValueAt(p);
}

Point ParticlePopulation::getRandomParticlePosition(Rng &rng) const {
    // TODO implementation
    return Point(0.0, 0.0, 0.0);
}

void ParticlePopulation::setCoordinateTransformation(const Aff_transformation
    &transformation) {
    densityDistribution_.setCoordinateTransformation(transformation);
}
void ParticlePopulation::removeCoordinateTransformation() {
    densityDistribution_.removeCoordinateTransformation();
}

