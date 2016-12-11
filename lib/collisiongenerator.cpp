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
    if (majorantReactionRate_ != NULL) {
        delete[] majorantReactionRate_;
    }
}

void CollisionGenerator::addCollisionReaction(CollisionReaction *reaction) {
    collisionReactions_.push_back(reaction);
}

void CollisionGenerator::precomputeReactionRates(double maxSpeed,
    double speedStepSize, simthreadresources &thread_res) {
    size_t nReactions = collisionReactions_.size();
    size_t gridSize = grid_.arraySize();

    speedStepSize_ = speedStepSize;
    nSpeedSteps_ = std::ceil(maxSpeed / speedStepSize);

    // Deallocate previously allocated resources
    if (totalReactionRate_ != NULL) {
        delete[] totalReactionRate_;
    }
    if (cumulativeProbability_ != NULL) {
        delete[] cumulativeProbability_;
    }
    if (majorantReactionRate_ != NULL) {
        delete[] majorantReactionRate_;
    }

    // Allocate necessary resources
    totalReactionRate_ = new double[nSpeedSteps_ * gridSize];
    cumulativeProbability_ =
        new double[nSpeedSteps_ * gridSize * nReactions];
    majorantReactionRate_ = new double[nSpeedSteps_];

    for (size_t iv = 0; iv < nSpeedSteps_; ++iv) {
        double *pTotalReactionRate = &totalReactionRate_[iv * gridSize];
        double *pCumulativeProbability =
            &cumulativeProbability_[iv * gridSize * nReactions];

        double particleSpeed = speedStepSize_ * iv;

        for (size_t i = 0; i < gridSize; ++i) {
            Point p;
            if (!grid_.getCellMidpoint(i, p)) continue;

            int j = 0;
            double totalRate = 0.0;
            double majorant = 0.0;
            for (CollisionReaction *reaction : collisionReactions_) {
                double relativeSpeed =
                    reaction->getPopulation().getRelativeSpeed(particleSpeed,
                    thread_res.ms, thread_res.gslrng);
                double rate = reaction->getMeanReactionRate(p, relativeSpeed);
                totalRate += rate;
                pCumulativeProbability[nReactions * i + j] = rate;
                majorant += reaction->getMajorantReactionRate(p, relativeSpeed);
                j += 1;
            }

            pTotalReactionRate[i] = totalRate;
            majorantReactionRate_[iv] = std::max(majorantReactionRate_[iv], majorant);

            for (int k = 0; k < j; ++k) {
                pCumulativeProbability[nReactions*i + k] /= totalRate;
            }
        }
    }
}

