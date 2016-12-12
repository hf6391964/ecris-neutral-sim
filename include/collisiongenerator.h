#pragma once

#include <vector>

#include "collisionreaction.h"
#include "grid.h"

class CollisionGenerator {
    std::vector<std::unique_ptr<CollisionReaction>> collisionReactions_;
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

        // NOTE: ownership of the pointed reaction will be transferred to the
        // instance of CollisionGenerator!
        void addCollisionReaction(CollisionReaction *reaction);
        void precomputeReactionRates(double maxSpeed,
            double speedStepSize, simthreadresources &thread_res);
};

