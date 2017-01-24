#include "recombination.h"
#include "flychkparser.h"

Recombination::Recombination(std::shared_ptr<ParticlePopulation> population,
    double reactionRate) :
    NeutralizationChannel(1.0 / reactionRate), population_(population) {}

Particle Recombination::sampleNeutralProduct(Rng &rng,
    const Particle &sourceParticle, double decayTime) const {
    Particle result(sourceParticle.getElement(),
        sourceParticle.getTime() + decayTime);
    Point pos = population_->getRandomParticlePosition(rng);
    Vector vel = population_->getRandomParticleVelocity(rng);
    result.setVelocity(vel);
    result.setPosition(pos);
    return result;
}

