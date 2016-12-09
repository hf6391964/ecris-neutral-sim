#pragma once

#include "collisionreaction.h"
#include "particlepopulation.h"

class ChargeExchangeReaction : CollisionReaction {
    double ionMeanSpeed_ = 0.0;
    double ionMajorantSpeed_ = 0.0;
    ParticlePopulation &population_;
    unsigned int chargeState_ = 0;
    double ionizationPotentialEv_ = 0.0;
    double projectileMass_ = 0.0;
    double mullerSalzbornCrossSection_ = 0.0;

    public:
        ChargeExchangeReaction(ParticlePopulation &population,
            double ionMeanSpeed, double ionMajorantSpeed,
            double ionizationPotentialEv);

        double getMeanReactionRate(const Point &p, double relativeSpeed)
            const;
        double getMajorantReactionRate(const Point &p,
            double relativeSpeed) const;

        std::vector<Particle> computeReactionProducts(Rng &rng,
            const Point &, const Particle &target) const;

        double getCrossSection(double) const;
};

