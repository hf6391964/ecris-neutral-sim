#pragma once

#include <vector>

#include "collisionreaction.h"
#include "grid.h"
#include "logger.h"

typedef std::vector<std::unique_ptr<CollisionReaction>> CollisionReactionVector;

class CollisionGenerator {
    CollisionReactionVector collisionReactions_;
    Grid grid_;
    double *rateCoefficients_ = NULL;
    double *majorantReactionRate_ = NULL;
    size_t nSpeedSteps_ = 0;
    double speedStepSize_ = 0.0;
    size_t nReactions_ = 0;

    void _cleanup();

    public:
        CollisionGenerator(const Grid &grid);
        ~CollisionGenerator();

        // NOTE: ownership of the pointed reaction will be transferred to the
        // instance of CollisionGenerator!
        void addCollisionReaction(CollisionReaction *reaction);
        void precomputeReactionRates(double maxSpeed,
            double speedStepSize, simthreadresources *thread_res);

        // Gives the mean free time for a particle with given speed
        double getMeanFreeTime(double particleSpeed) const;

        CollisionReaction *sampleCollision(Rng &rng, const Point &p,
            double particleSpeed, double dt) const;

        void writeStatistics(Logger &log) const;
};

