#pragma once

#include "collisionreaction.h"
#include "element_data.h"
#include "particlepopulation.h"

class ElectronIonizationReaction : public CollisionReaction {
    IonizationParameters ionizationParameters_;

    public:
        ElectronIonizationReaction(std::shared_ptr<ParticlePopulation> population,
            const IonizationParameters &ip);

        double getReactionRate(const Point &p, double particleSpeed,
            simthreadresources *thread_res) const;

        double getRateCoefficient(double particleSpeed,
            simthreadresources *thread_res) const;

        std::vector<Particle> computeReactionProducts(Rng &rng,
            const Point &p, const Particle &target) const;

        double getCrossSection(double speed) const;

        static double crossSection(double velocity, void *p);
};

