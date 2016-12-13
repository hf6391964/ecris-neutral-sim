#include "simpleplasmamodel.h"

#include "maxwellianpopulation.h"
#include "spatialdistribution.h"
#include "chargeexchangereaction.h"
#include "electronionizationreaction.h"

SimplePlasmaModel::SimplePlasmaModel(std::string electronDensityFilename,
    double electronWeight, std::vector<double> ionRelativeDensities,
    double electronTemperature, std::vector<double> ionTemperatures,
    const ElementData &elementData) : elementData_(elementData) {
    particlePopulations_.resize(ionRelativeDensities.size() + 1);

    // Calculate the density distributions
    // TODO insert check for validity of the distribution
    // i.e. if the reading succeeded
    DensityDistribution electronDensity(electronDensityFilename,
        electronWeight);

    ionMassEv_ = elementData_.mass * ATOMIC_MASS_TO_EV;
    // Populate the ion densities
    int i = 1;
    maxIonTemp_ = 0.0;
    for (double relDensity : ionRelativeDensities) {
        DensityDistribution ionDensityDistribution(electronDensity,
            relDensity);
        particlePopulations_[i] = std::shared_ptr<ParticlePopulation>(
            new MaxwellianPopulation(ionMassEv_, i,
            ionTemperatures[i], ionDensityDistribution));
        i += 1;
        maxIonTemp_ = std::max(maxIonTemp_, ionTemperatures[i]);
    }
    particlePopulations_[0] = std::shared_ptr<ParticlePopulation>(
        new MaxwellianPopulation(ELECTRON_MASS_EV, -1, electronTemperature,
            electronDensity));

    maxChargeState_ = particlePopulations_.size() - 1;
}

SimplePlasmaModel::~SimplePlasmaModel() {
}

void SimplePlasmaModel::populateCollisionReactions(
    CollisionGenerator &generator, simthreadresources *thread_res) const {
    std::cout << "Populating collision reactions...\n";
    // Populate the electron ionization reaction
    generator.addCollisionReaction(new ElectronIonizationReaction(
        particlePopulations_[0], elementData_.ionizationParameters));

    // Populate charge exchange reactions
    for (unsigned int q = 1; q <= maxChargeState_; ++q) {
        generator.addCollisionReaction(new ChargeExchangeReaction(
            particlePopulations_[q], elementData_.ionizationParameters.I_le));
    }

    double mbave = Util::getMBAverage(maxIonTemp_, ionMassEv_);

    generator.precomputeReactionRates(8.0*mbave, mbave/1000.0, thread_res);
}

void SimplePlasmaModel::setCoordinateTransformation(
    const Aff_transformation &tf) {
    for (auto iPopulation = particlePopulations_.begin();
         iPopulation != particlePopulations_.end(); ++iPopulation) {
        (*iPopulation)->setCoordinateTransformation(tf);
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

