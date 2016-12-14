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
    nSpeedSteps_ = std::ceil(maxSpeed / speedStepSize) + 1;

    _cleanup();

    size_t nBytes = nSpeedSteps_ * (1 + gridSize_ * (1 + nReactions_));
    std::cout << "Precomputed values will take up " << nBytes << " bytes\n";
    // Allocate necessary resources
    totalReactionRate_ = new double[nSpeedSteps_ * gridSize_];
    cumulativeProbability_ =
        new double[nSpeedSteps_ * gridSize_ * nReactions_];
    majorantReactionRate_ = new double[nSpeedSteps_];
    double *rateCoeffs = new double[nReactions_];

    std::vector<std::shared_ptr<ParticlePopulation>> reactionPopulations;
    for (size_t ir = 0; ir < nReactions_; ++ir) {
        reactionPopulations.push_back(collisionReactions_[ir]->getPopulation());
    }

    std::cout << "Precomputing collision rates...\n";
    for (size_t iv = 0; iv < nSpeedSteps_; ++iv) {
        double *pTotalReactionRate = &totalReactionRate_[iv * gridSize_];
        double *pCumulativeProbability =
            &cumulativeProbability_[iv * gridSize_ * nReactions_];
        majorantReactionRate_[iv] = 0.0;
        double particleSpeed = speedStepSize_ * iv;

        for (size_t ir = 0; ir < nReactions_; ++ir) {
            rateCoeffs[ir] = collisionReactions_[ir]->getRateCoefficient(
                particleSpeed, thread_res);
            /* std::cout << "rate coeff " << ir << ": " << rateCoeffs[ir] << std::endl; */
        }

        for (size_t i = 0; i < gridSize_; ++i) {
            Point p;
            if (!grid_.getCellMidpoint(i, p)) continue;

            double totalRate = 0.0;
            for (size_t ir = 0; ir < nReactions_; ++ir) {
                double rate = reactionPopulations[ir]->getDensityAt(p) *
                    rateCoeffs[ir];
                totalRate += rate;
                // Sum rates to get cumulative rate, will be normalized later
                // to yield probability
                pCumulativeProbability[nReactions_ * i + ir] = totalRate;
            }

            pTotalReactionRate[i] = totalRate;
            majorantReactionRate_[iv] = std::max(majorantReactionRate_[iv],
                totalRate);

            for (size_t k = 0; k < nReactions_; ++k) {
                // Normalize the previously stored cumulative rates into
                // probabilities
                pCumulativeProbability[nReactions_*i + k] /= totalRate;
            }
        }

        /* std::cout << "Majorant rate for v = " << particleSpeed << " is " << majorantReactionRate_[iv] << "\n"; */
    }

    delete[] rateCoeffs;

    std::cout << "Precomputation done.\n";
}

double CollisionGenerator::getMeanFreeTime(double particleSpeed) const {
    size_t velIndex = particleSpeed / speedStepSize_;
    if (velIndex >= nSpeedSteps_ - 1) return -1.0;

    double t = (particleSpeed - velIndex * speedStepSize_) / speedStepSize_;
    double interpolatedRate = (1.0 - t) * majorantReactionRate_[velIndex] +
        t * majorantReactionRate_[velIndex + 1];

    return 1.0 / interpolatedRate;
}

CollisionReaction *CollisionGenerator::sampleCollision(Rng &rng,
    const Point &p, double particleSpeed, double dt) const {
    size_t velIndex = particleSpeed / speedStepSize_;
    size_t spatialIndex;
    if (velIndex >= nSpeedSteps_ - 1 || !grid_.arrayIndex(p, spatialIndex)) {
        return NULL;
    }

    size_t velSpatialIndex = velIndex * gridSize_ + spatialIndex;
    size_t velSpatialIndex1 = (velIndex + 1) * gridSize_ + spatialIndex;
    double t = (particleSpeed - velIndex * speedStepSize_) / speedStepSize_;
    double totalRate = (1.0 - t) * totalReactionRate_[velSpatialIndex] +
        t * totalReactionRate_[velSpatialIndex1];

    double x = uni01(rng);
    if (x > std::exp(-totalRate * dt)) {  // Reaction happens
        double y = uni01(rng);
        double *pCumulativeProbability =
            &cumulativeProbability_[nReactions_ * velSpatialIndex];
        double *pCumulativeProbability1 =
            &cumulativeProbability_[nReactions_ * velSpatialIndex1];
        size_t iReaction;
        for (iReaction = 0; iReaction < nReactions_; ++iReaction) {
            double prob = (1.0 - t) * pCumulativeProbability[iReaction] +
                t * pCumulativeProbability1[iReaction];
            if (y < prob) {
                break;
            }
        }
        return collisionReactions_[iReaction].get();
    }

    return NULL;
}

