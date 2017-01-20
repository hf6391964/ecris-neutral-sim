#pragma once

#include <string>

#include "neutralizationchannel.h"
#include "particlepopulation.h"

class Recombination : NeutralizationChannel {
    std::shared_ptr<ParticlePopulation> population_;

    public:
        Recombination(std::shared_ptr<ParticlePopulation> population,
            double reactionRate);
};

