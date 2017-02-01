#include "neutralizationgenerator.h"
#include "logger.h"


void NeutralizationGenerator::addNeutralizationChannel(
    NeutralizationChannel *channel) {
    channels_.push_back(std::unique_ptr<NeutralizationChannel>(channel));
}

Particle NeutralizationGenerator::sampleNeutralizationReaction(Rng &rng,
    const Particle &sourceParticle) const {
    if (channels_.size() == 0) {
        throw std::length_error(
            "Cannot sample neutralization: no reactions added");
    }

    NeutralizationChannel *chosenChannel = channels_[0].get();
    double minTime = channels_[0]->timeToReaction(rng);
    double time;
    for (auto i = channels_.begin() + 1; i != channels_.end(); ++i) {
        time = (*i)->timeToReaction(rng);
        logger << "Neutralization channel " << (*i)->getLabel() <<
            ", decay time = " << time;
        if (time < minTime) {
            minTime = time;
            chosenChannel = (*i).get();
        }
    }

    logger << "Chosen neutralization channel " << chosenChannel->getLabel() <<
        ", decay time = " << minTime;

    chosenChannel->incrementReactionCounter();

    return chosenChannel->sampleNeutralProduct(rng, sourceParticle, minTime);
}

