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
    double speedStepSize, simthreadresources *thread_res) {

    nReactions_ = collisionReactions_.size();
    gridSize_ = grid_.arraySize();

    speedStepSize_ = speedStepSize;
    nSpeedSteps_ = std::ceil(maxSpeed / speedStepSize);

    _cleanup();

    size_t nBytes = nSpeedSteps_ * (1 + gridSize_ * (1 + nReactions_));
    std::cout << "Precomputed values will take up " << nBytes << " bytes\n";
    // Allocate necessary resources
    totalReactionRate_ = new double[nSpeedSteps_ * gridSize_];
    cumulativeProbability_ =
        new double[nSpeedSteps_ * gridSize_ * nReactions_];
    majorantReactionRate_ = new double[nSpeedSteps_];

    std::cout << "Precomputing collision rates...\n";
    for (size_t iv = 0; iv < nSpeedSteps_; ++iv) {
        double *pTotalReactionRate = &totalReactionRate_[iv * gridSize_];
        double *pCumulativeProbability =
            &cumulativeProbability_[iv * gridSize_ * nReactions_];
        majorantReactionRate_[iv] = 0.0;
        double particleSpeed = speedStepSize_ * iv;

        for (size_t i = 0; i < gridSize_; ++i) {
            Point p;
            if (!grid_.getCellMidpoint(i, p)) continue;

            int j = 0;
            double totalRate = 0.0;
            for (auto iReaction = collisionReactions_.begin();
                 iReaction != collisionReactions_.end(); ++iReaction) {
                double rate = (*iReaction)->getReactionRate(p, particleSpeed,
                    thread_res);
                totalRate += rate;
                // Sum rates to get cumulative rate, will be normalized later
                // to yield probability
                pCumulativeProbability[nReactions_ * i + j] = totalRate;
                j += 1;
            }

            pTotalReactionRate[i] = totalRate;
            majorantReactionRate_[iv] = std::max(majorantReactionRate_[iv],
                totalRate);

            for (int k = 0; k < j; ++k) {
                // Normalize the previously stored cumulative rates into
                // probabilities
                pCumulativeProbability[nReactions_*i + k] /= totalRate;
            }
        }
    }

    std::cout << "Precomputation done.\n";
}

double CollisionGenerator::getMeanFreeTime(double particleSpeed) const {
    size_t velIndex = particleSpeed / speedStepSize_;
    if (velIndex >= nSpeedSteps_) return -1.0;

    return 1.0 / majorantReactionRate_[velIndex];
}

CollisionReaction *CollisionGenerator::sampleCollision(Rng &rng,
    const Point &p, double particleSpeed, double dt) const {
    size_t velIndex = particleSpeed / speedStepSize_;
    size_t spatialIndex;
    if (velIndex >= nSpeedSteps_ || grid_.arrayIndex(p, spatialIndex)) {
        return NULL;
    }

    size_t velSpatialIndex = velIndex * gridSize_ + spatialIndex;
    double totalRate = totalReactionRate_[velSpatialIndex];

    double x = uni01(rng);
    if (x > std::exp(-totalRate * dt)) {  // Reaction happens
        double y = uni01(rng);
        double *pCumulativeProbability =
            &cumulativeProbability_[nReactions_ * velSpatialIndex];
        size_t iReaction;
        for (iReaction = 0; iReaction < nReactions_; ++iReaction) {
            if (y < pCumulativeProbability[iReaction]) {
                break;
            }
        }
        return collisionReactions_[iReaction].get();
    }

    return NULL;
}

