#include "neutralizationchannel.h"


double NeutralizationChannel::timeToReaction(Rng &rng) const {
    double x = uni01(rng);

    return -decayTime_ * std::log(x * decayTime_);
}

std::string NeutralizationChannel::getLabel() const {
    return label_;
}

