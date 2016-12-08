#pragma once

#include "collisionreaction.h"
#include "element_data.h"

class ElectronIonizationReaction : CollisionReaction {
    double electronMeanSpeed_ = 0.0;
    double electronMajorantSpeed_ = 0.0;
    PlasmaDensities &plasmaDensities_;
    IonizationParameters ionizationParameters_;

    public:
        ElectronIonizationReaction(PlasmaDensities &densities,
            double electronMeanSpeed, double electronMajorantSpeed,
            const IonizationParameters &ip);

        double getMeanReactionRate(const Point &p, double relativeSpeed) const;
        double getMajorantReactionRate(const Point &p, double relativeSpeed)
            const;

        std::vector<Particle> computeReactionProducts(Rng &rng,
            const Point &p, const Particle &target) const;

        double getCrossSection(double speed) const;

        static double crossSection(double velocity, void *p);
};
