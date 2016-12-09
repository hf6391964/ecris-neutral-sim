#include "simpleplasmamodel.h"

SimplePlasmaModel::SimplePlasmaModel(std::string electronDensityFilename,
    double electronWeight, std::vector<double> ionRelativeDensities,
    double electronTemperature, std::vector<double> ionTemperatures,
    double ionMassEv) : ionMassEv_(ionMassEv) {
    particlePopulations_.resize(ionRelativeDensities.size() + 1);

    // Calculate the density distributions
    // TODO insert check for validity of the distribution
    // i.e. if the reading succeeded
    DensityDistribution electronDensity(electronDensityFilename,
        electronWeight);

    // Populate the ion densities
    int i = 1;
    for (double relDensity : ionRelativeDensities) {
        DensityDistribution ionDensityDistribution(electronDensity,
            relDensity);
        particlePopulations_[i] = new MaxwellianPopulation(ionMassEv, i,
            ionTemperatures[i], ionDensityDistribution);
        i += 1;
    }
    particlePopulations_[0] =
        new MaxwellianPopulation(ELECTRON_MASS_EV, -1, electronTemperature,
            electronDensity);

    maxChargeState_ = particlePopulations_.size() - 1;
}

SimplePlasmaModel::~SimplePlasmaModel() {
    for (ParticlePopulation *population : particlePopulations_) {
        if (population != NULL) {
            delete population;
        }
    }
}

void SimplePlasmaModel::setCoordinateTransformation(
    const Aff_transformation &tf) {
    for (ParticlePopulation *population : particlePopulations_) {
        population->setCoordinateTransformation(tf);
    }
}

double SimplePlasmaModel::getIonDensityAt(const Point &p,
    unsigned int chargeState) const {
    if (chargeState > maxChargeState_) return 0.0;
    return particlePopulations_[chargeState]->getDensityAt(p);
}

double SimplePlasmaModel::getElectronDensityAt(const Point &p) const {
    return getIonDensityAt(p, 0);
}

Vector SimplePlasmaModel::getIsotropicIonVelocity(Rng &rng,
    unsigned int chargeState) const {
    if (chargeState > maxChargeState_) return Vector(0.0, 0.0, 0.0);
    return particlePopulations_[chargeState]->getRandomParticleVelocity(rng);
}

Vector SimplePlasmaModel::getIsotropicElectronVelocity(Rng &rng) const {
    return getIsotropicIonVelocity(rng, 0);
}

