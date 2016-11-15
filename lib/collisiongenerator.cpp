#include "collisiongenerator.h"

void CollisionGenerator::addCollisionReaction(CollisionReaction *reaction) {
    collisionReactions_.push_back(reaction);
}

