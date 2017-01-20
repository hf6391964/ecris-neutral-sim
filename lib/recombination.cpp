#include "recombination.h"
#include "flychkparser.h"

Recombination::Recombination(std::shared_ptr<ParticlePopulation> population,
    double reactionRate) :
    NeutralizationChannel(1.0 / reactionRate), population_(population) {}

