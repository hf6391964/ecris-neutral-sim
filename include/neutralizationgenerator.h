#pragma once

#include "neutralizationchannel.h"

class NeutralizationGenerator {
    std::vector<std::unique_ptr<NeutralizationChannel>> channels_;

    public:
        void addNeutralizationChannel(NeutralizationChannel *channel);

        Particle sampleNeutralizationReaction(Rng &rng,
            const Particle &sourceParticle) const;
};

