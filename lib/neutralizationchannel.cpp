#include "neutralizationchannel.h"


double NeutralizationChannel::timeToReaction(Rng &rng) const {
    double x = uni01(rng);

    return -decayTime_ * std::log(x * decayTime_);
}

