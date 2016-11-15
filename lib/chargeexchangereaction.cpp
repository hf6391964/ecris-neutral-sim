#include "chargeexchangereaction.h"

ChargeExchangeReaction::ChargeExchangeReaction(PlasmaDensities &densities,
    double ionMeanSpeed, double ionMajorantSpeed,
    unsigned int chargeState, double ionizationPotentialEv)
    : ionMeanSpeed_(ionMeanSpeed), ionMajorantSpeed_(ionMajorantSpeed),
      plasmaDensities_(densities),
      chargeState_(chargeState), ionizationPotentialEv_(ionizationPotentialEv) {
    mullerSalzbornCrossSection_ = 1.43e-12 *
        std::pow((double)chargeState_, 1.17) *
        std::pow(ionizationPotentialEv_, -2.76);
}

double ChargeExchangeReaction::getCrossSection(const Point &, double) const {
    return mullerSalzbornCrossSection_;
}

double ChargeExchangeReaction::getMeanReactionRate(const Point &p) const {
    return getCrossSection(p, ionMeanSpeed_) * ionMeanSpeed_ *
        plasmaDensities_.getIonDensityAt(p, chargeState_);
}

double ChargeExchangeReaction::getMajorantReactionRate(const Point &p) const {
    return getCrossSection(p, ionMajorantSpeed_) * ionMajorantSpeed_ *
        plasmaDensities_.getIonDensityAt(p, chargeState_);
}

