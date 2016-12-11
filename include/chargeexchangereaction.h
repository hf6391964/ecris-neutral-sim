#pragma once

#include "collisionreaction.h"
#include "particlepopulation.h"

class ChargeExchangeReaction : public CollisionReaction {
    double ionMeanSpeed_ = 0.0;
    double ionMajorantSpeed_ = 0.0;
    unsigned int chargeState_ = 0;
    double ionizationPotentialEv_ = 0.0;
    double projectileMass_ = 0.0;
    double mullerSalzbornCrossSection_ = 0.0;

    public:
        ChargeExchangeReaction(const ParticlePopulation &population,
            double ionMeanSpeed, double ionMajorantSpeed,
            double ionizationPotentialEv);

        double getReactionRate(const Point &p, double particleSpeed,
            simthreadresources &thread_res) const;

        std::vector<Particle> computeReactionProducts(Rng &rng,
            const Point &, const Particle &target) const;

        double getCrossSection(double) const;

        static double crossSection(double, void *args);
};

