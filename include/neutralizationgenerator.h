#pragma once

#include "neutralizationchannel.h"
#include "logger.h"

typedef std::vector<std::unique_ptr<NeutralizationChannel>>
    NeutralizationChannelVector;

class NeutralizationGenerator {
    NeutralizationChannelVector channels_;

    public:
        void addNeutralizationChannel(
            std::unique_ptr<NeutralizationChannel> channel);

        Particle sampleNeutralizationReaction(Rng &rng,
            const Particle &sourceParticle) const;

        void writeStatistics(Logger& log) const;
};

