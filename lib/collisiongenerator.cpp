#include "collisiongenerator.h"

CollisionGenerator::CollisionGenerator(const Grid &grid)
    : grid_(grid) {
}

CollisionGenerator::~CollisionGenerator() {
    _cleanup();
}

void CollisionGenerator::_cleanup() {
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
}

void CollisionGenerator::addCollisionReaction(CollisionReaction *reaction) {
    collisionReactions_.push_back(std::unique_ptr<CollisionReaction>(reaction));
}

void CollisionGenerator::precomputeReactionRates(double maxSpeed,
    double speedStepSize, simthreadresources &thread_res) {
    size_t nReactions = collisionReactions_.size();
    size_t gridSize = grid_.arraySize();

    speedStepSize_ = speedStepSize;
    nSpeedSteps_ = std::ceil(maxSpeed / speedStepSize);

    _cleanup();

    // Allocate necessary resources
    totalReactionRate_ = new double[nSpeedSteps_ * gridSize];
    cumulativeProbability_ =
        new double[nSpeedSteps_ * gridSize * nReactions];
    majorantReactionRate_ = new double[nSpeedSteps_];

    for (size_t iv = 0; iv < nSpeedSteps_; ++iv) {
        double *pTotalReactionRate = &totalReactionRate_[iv * gridSize];
        double *pCumulativeProbability =
            &cumulativeProbability_[iv * gridSize * nReactions];
        majorantReactionRate_[iv] = 0.0;
        double particleSpeed = speedStepSize_ * iv;

        for (size_t i = 0; i < gridSize; ++i) {
            Point p;
            if (!grid_.getCellMidpoint(i, p)) continue;

            int j = 0;
            double totalRate = 0.0;
            for (auto iReaction = collisionReactions_.begin();
                 iReaction != collisionReactions_.end(); ++i) {
                double rate = (*iReaction)->getReactionRate(p, particleSpeed,
                    thread_res);
                totalRate += rate;
                // Sum rates to get cumulative rate, will be normalized later
                // to yield probability
                pCumulativeProbability[nReactions * i + j] = totalRate;
                j += 1;
            }

            pTotalReactionRate[i] = totalRate;
            majorantReactionRate_[iv] = std::max(majorantReactionRate_[iv],
                totalRate);

            for (int k = 0; k < j; ++k) {
                // Normalize the previously stored cumulative rates into
                // probabilities
                pCumulativeProbability[nReactions*i + k] /= totalRate;
            }
        }
    }
}

double CollisionGenerator::getMeanFreeTime(double particleSpeed) const {
    size_t velIndex = particleSpeed / speedStepSize_;
    if (velIndex >= nSpeedSteps_) return -1.0;

    return 1.0 / majorantReactionRate_[velIndex];
}

