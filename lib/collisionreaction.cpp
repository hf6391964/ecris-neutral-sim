#include "collisionreaction.h"

CollisionReaction::CollisionReaction(
    std::shared_ptr<ParticlePopulation> population) : population_(population) {
    reactionCounter_ = 0;
}

std::shared_ptr<ParticlePopulation> CollisionReaction::getPopulation() const {
    return population_;
}

std::string CollisionReaction::getLabel() const {
    return label_ + " " + population_->getLabel();
}

void CollisionReaction::incrementReactionCounter() {
    reactionCounter_ += 1;
}

unsigned long CollisionReaction::getReactionCount() const {
    return reactionCounter_;
}

