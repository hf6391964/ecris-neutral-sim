#include "collisiongenerator.h"

CollisionGenerator::CollisionGenerator(const Grid &grid)
    : grid_(grid) {
}

void CollisionGenerator::addCollisionReaction(std::unique_ptr<CollisionReaction> reaction) {
    collisionReactions_.push_back(std::move(reaction));
}

void CollisionGenerator::precomputeReactionRates(double maxSpeed,
    double speedStepSize, uint_least32_t seed) {

    nReactions_ = collisionReactions_.size();
    size_t gridSize = grid_.arraySize();

    speedStepSize_ = speedStepSize;
    nSpeedSteps_ = std::ceil(maxSpeed / speedStepSize) + 1;

    size_t nBytes = nSpeedSteps_ * (1 + nReactions_) * sizeof(double);
    std::cout << "Precomputed values will take up " << nBytes << " bytes\n";
    // Allocate necessary resources
    majorantReactionRate_ = std::vector<double>(nSpeedSteps_);
    rateCoefficients_ = RateCoefficientVector(nSpeedSteps_);

    std::vector<double> densityMaxima(nReactions_, 0.0);
    for (size_t i = 0; i < gridSize; ++i) {
        Point p;
        if (!grid_.getCellMidpoint(i, p)) continue;
        for (size_t ir = 0; ir < nReactions_; ++ir) {
            densityMaxima.at(ir) = std::max(densityMaxima.at(ir),
                collisionReactions_.at(ir)->getPopulation()->getDensityAt(p));
        }
    }

    std::cout << "Precomputing collision rates...\n";
    mc_integrate_resources mc_res(seed);
    for (size_t iv = 0; iv < nSpeedSteps_; ++iv) {
        double particleSpeed = speedStepSize_ * iv;
        rateCoefficients_.at(iv) = std::vector<double>(nReactions_);

        // Calculate rate coefficients separately for each reaction
        // Find the majorant reaction rate at the given velocity
        majorantReactionRate_.at(iv) = 0.0;
        for (size_t ir = 0; ir < nReactions_; ++ir) {
            rateCoefficients_.at(iv).at(ir) =
                collisionReactions_.at(ir)->getRateCoefficient(particleSpeed, mc_res);
            majorantReactionRate_.at(iv) += densityMaxima.at(ir) *
                rateCoefficients_.at(iv).at(ir);
        }
    }

    std::cout << "Precomputation done.\n";
}

double CollisionGenerator::getMeanFreeTime(double particleSpeed) const {
    size_t velIndex = particleSpeed / speedStepSize_;
    if (velIndex >= nSpeedSteps_ - 1) {
        std::string error = "Couldn't get mean free time for particle speed ";
        error += std::to_string(particleSpeed);
        throw std::out_of_range(error);
    }

    double t = (particleSpeed - velIndex * speedStepSize_) / speedStepSize_;
    double interpolatedRate = (1.0 - t) * majorantReactionRate_.at(velIndex) +
        t * majorantReactionRate_.at(velIndex + 1);

    return 1.0 / interpolatedRate;
}

CollisionReaction *CollisionGenerator::sampleCollision(Rng &rng,
    const Point &p, double particleSpeed, double dt) const {
    size_t velIndex = particleSpeed / speedStepSize_;
    size_t spatialIndex;
    if (velIndex >= nSpeedSteps_ - 1 || !grid_.arrayIndex(p, spatialIndex)) {
        return nullptr;
    }

    double totalReactionRate1 = 0.0, totalReactionRate2 = 0.0;
    std::vector<double> cumulativeReactionRate1(nReactions_),
        cumulativeReactionRate2(nReactions_);
    for (size_t ir = 0; ir < nReactions_; ++ir) {
        double density = collisionReactions_.at(ir)->getPopulation()->getDensityAt(p);
        double rate1 = density * rateCoefficients_.at(velIndex).at(ir),
            rate2 = density * rateCoefficients_.at(velIndex+1).at(ir);
        totalReactionRate1 += rate1;
        totalReactionRate2 += rate2;
        cumulativeReactionRate1.at(ir) = totalReactionRate1;
        cumulativeReactionRate2.at(ir) = totalReactionRate2;
    }

    double t = (particleSpeed - velIndex * speedStepSize_) / speedStepSize_;
    double totalRate = (1.0 - t) * totalReactionRate1 + t * totalReactionRate2;

    double x = uni01(rng);
    if (x > std::exp(-totalRate * dt)) {  // Reaction happens
        double y = uni01(rng) * totalRate;
        size_t iReaction;
        for (iReaction = 0; iReaction < nReactions_; ++iReaction) {
            double probability = (1.0 - t) * cumulativeReactionRate1.at(iReaction) +
                t * cumulativeReactionRate2.at(iReaction);
            if (y < probability) {
                break;
            }
        }
        collisionReactions_.at(iReaction)->incrementReactionCounter();
        return collisionReactions_.at(iReaction).get();
    }

    return nullptr;
}

void CollisionGenerator::writeStatistics(Logger &log) const {
    log << "Collision reaction statistics:\n";
    for (CollisionReactionVector::const_iterator it = collisionReactions_.begin();
        it != collisionReactions_.end(); ++it) {
        log << (*it)->getLabel() << ": " << (*it)->getReactionCount() <<
            " reactions\n";
    }
    log << '\n';
}

