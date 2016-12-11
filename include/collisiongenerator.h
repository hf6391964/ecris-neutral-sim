#pragma once

#include <vector>

#include "collisionreaction.h"
#include "grid.h"

class CollisionGenerator {
    std::vector<CollisionReaction *> collisionReactions_;
    Grid grid_;
    double *totalReactionRate_ = NULL;
    double *cumulativeProbability_ = NULL;
    double *majorantReactionRate_ = NULL;
    size_t nSpeedSteps_ = 0;
    double speedStepSize_ = 0.0;

    void _cleanup();

    public:
        CollisionGenerator(const Grid &grid);
        ~CollisionGenerator();

        void addCollisionReaction(CollisionReaction *);
        void precomputeReactionRates(double maxSpeed,
            double speedStepSize, simthreadresources &thread_res);
};

