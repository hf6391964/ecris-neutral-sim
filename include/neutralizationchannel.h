#pragma once

#include "cgal_and_typedefs.h"

class NeutralizationChannel {
    double decayTime_ = 0.0;

    public:
        NeutralizationChannel(double decayTime) : decayTime_(decayTime) {}

        double timeToReaction(Rng &rng) const;
};

