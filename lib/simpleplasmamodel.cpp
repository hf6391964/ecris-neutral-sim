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

    ionMassEv_ = elementData_->mass * ATOMIC_MASS_TO_EV;
    maxIonTemp_ = *std::max_element(ionTemperatures.begin(),
        ionTemperatures.end());
    particlePopulations_.at(0) = std::make_shared<MaxwellianPopulation>(
        ELECTRON_MASS_EV, -1, electronTemperature,
        std::make_shared<DensityDistribution>(electronDensityFilename,
        electronWeight));
    // Populate the ion densities
    int q = 1;
    for (double relDensity : ionRelativeDensities) {
        particlePopulations_.at(q) = std::make_shared<MaxwellianPopulation>(
            element_, q, ionTemperatures.at(q - 1),
            std::make_shared<DensityDistribution>(
                particlePopulations_.at(0)->getDensityDistribution(), relDensity));
        q += 1;
    }

    maxChargeState_ = ionRelativeDensities.size();
}

void SimplePlasmaModel::populateCollisionReactions(
    CollisionGenerator &generator, uint_least32_t seed, double
    speedStepSize) const {
    // Populate the electron ionization reaction
    generator.addCollisionReaction(std::make_unique<ElectronIonizationReaction>(
        particlePopulations_.at(0), elementData_->ionizationParameters));

    // Populate charge exchange reactions
    for (unsigned int q = 1; q <= maxChargeState_; ++q) {
        generator.addCollisionReaction(std::make_unique<ChargeExchangeReaction>(
            particlePopulations_.at(q), elementData_->ionizationParameters.I_le));
    }

    double mbave = Util::getMBAverage(maxIonTemp_, ionMassEv_);

    generator.precomputeReactionRates(8.0*mbave, speedStepSize * mbave,
        seed);
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
    generator.addNeutralizationChannel(std::make_unique<WallNeutralization>(
        particlePopulations_.at(1), confinementTime, wallFilename,
        end1Filename, end2Filename, surfaces));

    // Recombination: sample recombination rate from FLYCHK data
    FlychkParser parser(flychkFilename);
    double recombinationRateCoefficient;
    if (!parser.getTotalRateCoefficient(
        particlePopulations_.at(0)->getTemperature(), 1,
        recombinationRateCoefficient)) {
        throw std::invalid_argument(
            "Couldn't find recombination rate coefficient from FLYCHK");
    }
    double recombinationRate =
        averageElectronDensity * recombinationRateCoefficient;
    generator.addNeutralizationChannel(std::make_unique<Recombination>(
        particlePopulations_.at(1), recombinationRate));
}

