#pragma once

#include <vector>

#include "collisionreaction.h"

class CollisionGenerator {
    std::vector<CollisionReaction *> collisionReactions_;
    Grid grid_;
    double *totalReactionRate_ = NULL;
    double *cumulativeProbability_ = NULL;
    double majorantReactionRate_ = 0.0;

    public:
        CollisionGenerator(const Grid &grid);
        ~CollisionGenerator();

        void addCollisionReaction(CollisionReaction *);
        void precomputeReactions();
};

