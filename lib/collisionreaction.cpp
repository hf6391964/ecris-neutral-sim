#include "collisionreaction.h"

CollisionReaction::CollisionReaction(
    const ParticlePopulation &population) : reactionCounter_(0),
    population_(population) {}

const ParticlePopulation &CollisionReaction::getPopulation() const {
    return population_;
}

std::string CollisionReaction::getLabel() const {
    return label_ + " (" + population_.getLabel() + ')';
}

void CollisionReaction::incrementReactionCounter() {
    reactionCounter_ += 1;
}

unsigned long CollisionReaction::getReactionCount() const {
    return reactionCounter_;
}

