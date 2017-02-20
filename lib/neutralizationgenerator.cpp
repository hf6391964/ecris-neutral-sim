#include "neutralizationgenerator.h"
#include "logger.h"


void NeutralizationGenerator::addNeutralizationChannel(
    std::unique_ptr<NeutralizationChannel> channel) {
    channels_.push_back(std::move(channel));
}

Particle NeutralizationGenerator::sampleNeutralizationReaction(Rng &rng,
    const Particle &sourceParticle) const {
    if (channels_.size() == 0) {
        throw std::length_error(
            "Cannot sample neutralization: no reactions added");
    }

    NeutralizationChannel *chosenChannel = NULL;
    double maxRate = 0.0;
    double time;
    for (auto i = channels_.begin(); i != channels_.end(); ++i) {
        time = (*i)->timeToReaction(rng);
        logger << "Neutralization channel " << (*i)->getLabel() <<
            ", decay time = " << time << '\n';
        double rate = 1.0 / time;
        if (1.0 / time > maxRate) {
            maxRate = rate;
            chosenChannel = (*i).get();
        }
    }

    if (chosenChannel == NULL) {
        throw std::range_error("No neutralization reactions specified, can't continue");
    }

    double minTime = 1.0 / maxRate;
    logger << "Chosen neutralization channel " << chosenChannel->getLabel() <<
        ", decay time = " << minTime << '\n';

    chosenChannel->incrementReactionCounter();

    return chosenChannel->sampleNeutralProduct(rng, sourceParticle, minTime);
}

void NeutralizationGenerator::writeStatistics(Logger& log) const {
    log << "Neutralization channel statistics:\n";
    for (NeutralizationChannelVector::const_iterator it = channels_.begin();
        it != channels_.end(); ++it) {
        log << (*it)->getLabel() << ": " << (*it)->getReactionCount() <<
            " neutralization reactions\n";
    }
    log << '\n';
}

