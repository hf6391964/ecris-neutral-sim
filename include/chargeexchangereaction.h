#pragma once

#include "collisionreaction.h"

class ChargeExchangeReaction : CollisionReaction {
    double ionMeanSpeed_ = 0.0;
    double ionMajorantSpeed_ = 0.0;
    PlasmaDensities &plasmaDensities_;
    unsigned int chargeState_ = 0;
    double ionizationPotentialEv_ = 0.0;
    double mullerSalzbornCrossSection_ = 0.0;

    public:
        ChargeExchangeReaction(PlasmaDensities &densities,
            double ionMeanSpeed, double ionMajorantSpeed,
            unsigned int chargeState, double ionizationPotentialEv);

        double getMeanReactionRate(const Point &p) const;
        double getMajorantReactionRate(const Point &p) const;

        double getCrossSection(const Point &, double) const;
};

