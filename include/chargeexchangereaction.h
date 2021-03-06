#pragma once

#include "collisionreaction.h"
#include "particlepopulation.h"

class ChargeExchangeReaction : public CollisionReaction {
    unsigned int chargeState_ = 0;
    double ionizationPotentialEv_ = 0.0;
    double mullerSalzbornCrossSection_ = 0.0;

    public:
        ChargeExchangeReaction(
            const std::shared_ptr<ParticlePopulation> &population,
            double ionizationPotentialEv);

        double getReactionRate(const Point &p, double particleSpeed,
            mc_integrate_resources &mc_res) const;

        double getRateCoefficient(double particleSpeed,
            mc_integrate_resources &mc_res) const;

        CollisionProducts computeReactionProducts(Rng &rng,
            const Particle &target) const;

        double getCrossSection(double) const;

        static double crossSection(double, void *args);
};

