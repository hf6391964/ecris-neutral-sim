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

std::vector<Particle> ChargeExchangeReaction::computeReactionProducts(
    Rng &rng, const Point &, const Particle &target) const {
    // Elastic collision kinematics calculated in center of mass frame

    // Projectile particle velocity is isotropic Maxwell-Boltzmann:
    Vector ionVelocity = plasmaDensities_.getIsotropicIonVelocity(
        rng, chargeState_);
    Vector neutralVelocity = target.getVelocity();
    double ionMass = plasmaDensities_.getIonMass();
    double neutralMass = target.getMass_eV();

    Vector vcm = (ionMass * ionVelocity + neutralMass * neutralVelocity) /
        (ionMass + neutralMass);

    Vector neutralVelocityCm = neutralVelocity - vcm;

    double momentumNorm = std::sqrt(neutralVelocityCm.squared_length() /
        neutralMass);

    Direction dir = Util::getIsotropicSphereDirection(rng);
    Vector productVel = vcm + dir.vector() * momentumNorm / ionMass;

    Particle neutralProduct;
    neutralProduct.setMass_eV(ionMass);
    neutralProduct.setVelocity(productVel);

    std::vector<Particle> products;
    products.push_back(neutralProduct);

    return products;
}

