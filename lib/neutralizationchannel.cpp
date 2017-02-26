#include "neutralizationchannel.h"


double NeutralizationChannel::timeToReaction(Rng &rng) const {
    return -decayTime_ * std::log(1.0 - uni01(rng));
}

std::string NeutralizationChannel::getLabel() const {
    return label_;
}

void NeutralizationChannel::incrementReactionCounter() {
    reactionCounter_ += 1;
}

unsigned long NeutralizationChannel::getReactionCount() const {
    return reactionCounter_;
}

