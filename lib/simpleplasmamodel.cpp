#include "simpleplasmamodel.h"

#include "maxwellianpopulation.h"
#include "spatialdistribution.h"
#include "chargeexchangereaction.h"
#include "electronionizationreaction.h"
#include "flychkparser.h"
#include "wallneutralization.h"
#include "recombination.h"

SimplePlasmaModel::SimplePlasmaModel(std::string electronDensityFilename,
    double electronWeight, std::vector<double> ionRelativeDensities,
    double electronTemperature, std::vector<double> ionTemperatures,
    Element element) : element_(element) {
    elementData_ = ELEMENT_DATA.at(element_);

    particlePopulations_.resize(ionRelativeDensities.size() + 1);

    // Calculate the density distributions
    // TODO insert check for validity of the distribution
    // i.e. if the reading succeeded
    DensityDistribution electronDensity(electronDensityFilename,
        electronWeight);

    ionMassEv_ = elementData_->mass;
    // Populate the ion densities
    int q = 1;
    maxIonTemp_ = 0.0;
    for (double relDensity : ionRelativeDensities) {
        DensityDistribution ionDensityDistribution(electronDensity,
            relDensity);
        particlePopulations_[q] = std::shared_ptr<MaxwellianPopulation>(
            new MaxwellianPopulation(element_, q,
            ionTemperatures[q - 1], ionDensityDistribution));
        maxIonTemp_ = std::max(maxIonTemp_, ionTemperatures[q - 1]);
        q += 1;
    }
    particlePopulations_[0] = std::shared_ptr<MaxwellianPopulation>(
        new MaxwellianPopulation(ELECTRON_MASS_EV, -1, electronTemperature,
            electronDensity));

    maxChargeState_ = ionRelativeDensities.size();
}

SimplePlasmaModel::~SimplePlasmaModel() {
}

void SimplePlasmaModel::populateCollisionReactions(
    CollisionGenerator &generator, simthreadresources *thread_res) const {
    // Populate the electron ionization reaction
    generator.addCollisionReaction(new ElectronIonizationReaction(
        particlePopulations_[0], elementData_->ionizationParameters));

    // Populate charge exchange reactions
    for (unsigned int q = 1; q <= maxChargeState_; ++q) {
        generator.addCollisionReaction(new ChargeExchangeReaction(
            particlePopulations_[q], elementData_->ionizationParameters.I_le));
    }

    double mbave = Util::getMBAverage(maxIonTemp_, ionMassEv_);

    generator.precomputeReactionRates(8.0*mbave, mbave/10.0, thread_res);
}

void SimplePlasmaModel::populateNeutralizationReactions(
    NeutralizationGenerator &generator,
    double confinementTime,
    const std::string &wallFilename,
    const std::string &end1Filename,
    const std::string &end2Filename,
    const std::string &flychkFilename,
    const SurfaceCollection &surfaces,
    double averageElectronDensity) const {
    // Wall neutralization
    generator.addNeutralizationChannel(new WallNeutralization(
        particlePopulations_[1], confinementTime, wallFilename,
        end1Filename, end2Filename, surfaces));

    // Recombination: sample recombination rate from FLYCHK data
    FlychkParser parser(flychkFilename);
    double recombinationRateCoefficient;
    if (!parser.getTotalRateCoefficient(
        particlePopulations_[1]->getTemperature(), 1,
        recombinationRateCoefficient)) {
        throw std::invalid_argument(
            "Couldn't find recombination rate coefficient from FLYCHK");
    }
    double recombinationRate =
        averageElectronDensity * recombinationRateCoefficient;
    generator.addNeutralizationChannel(
        new Recombination(particlePopulations_[1], recombinationRate));
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

