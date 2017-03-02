#pragma once

#include "collisionreaction.h"
#include "element_data.h"
#include "particlepopulation.h"

class ElectronIonizationReaction : public CollisionReaction {
    IonizationParameters ionizationParameters_;

    public:
        ElectronIonizationReaction(
            const std::shared_ptr<ParticlePopulation> &population,
            const IonizationParameters &ip);

        double getReactionRate(const Point &p, double particleSpeed,
            mc_integrate_resources &mc_res) const;

        double getRateCoefficient(double particleSpeed,
            mc_integrate_resources &mc_res) const;

        CollisionProducts computeReactionProducts(Rng &rng,
            const Point &p, const Particle &target) const;

        double getCrossSection(double speed) const;

        static double crossSection(double velocity, void *p);
};

