#include "collisiongenerator.h"

CollisionGenerator::CollisionGenerator(const Grid &grid)
    : grid_(grid) {
}

CollisionGenerator::~CollisionGenerator() {
    if (totalReactionRate_ != NULL) {
        delete[] totalReactionRate_;
    }
    if (cumulativeProbability_ != NULL) {
        delete[] cumulativeProbability_;
    }
}

void CollisionGenerator::addCollisionReaction(CollisionReaction *reaction) {
    collisionReactions_.push_back(reaction);
}

void CollisionGenerator::precomputeReactions() {
    size_t nReactions = collisionReactions_.size();
    size_t gridSize = grid_.arraySize();

    // Deallocate previously allocated resources
    if (totalReactionRate_ != NULL) {
        delete[] totalReactionRate_;
    }
    if (cumulativeProbability_ != NULL) {
        delete[] cumulativeProbability_;
    }

    // Allocate necessary resources
    totalReactionRate_ = new double[gridSize];
    cumulativeProbability_ = new double[gridSize * nReactions];
    majorantReactionRate_ = 0.0;

    // First populate cumulative probability array with reaction rates.
    // It will be normalized later
    for (size_t i = 0; i < gridSize; i++) {
        Point p;
        if (!grid_.getCellMidpoint(i, p)) continue;

        int j = 0;
        double totalRate = 0.0;
        double majorant = 0.0;
        for (CollisionReaction *reaction : collisionReactions_) {
            double rate = reaction->getMeanReactionRate(p, 0.0);
            totalRate += rate;
            cumulativeProbability_[nReactions * i + j] = rate;
            majorant += reaction->getMajorantReactionRate(p, 0.0);
            j += 1;
        }

        totalReactionRate_[i] = totalRate;
        majorantReactionRate_ = std::max(majorantReactionRate_, majorant);

        for (int k = 0; k < j; k++) {
            cumulativeProbability_[nReactions*i + k] /= totalRate;
        }
    }
}

