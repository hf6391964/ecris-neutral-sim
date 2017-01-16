#include "collisionreaction.h"

CollisionReaction::CollisionReaction(
    std::shared_ptr<ParticlePopulation> population) : population_(population) {}

std::shared_ptr<ParticlePopulation> CollisionReaction::getPopulation() const {
    return population_;
}
