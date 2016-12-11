#include "chargeexchangereaction.h"

ChargeExchangeReaction::ChargeExchangeReaction(
    const ParticlePopulation &population,
    double ionMeanSpeed, double ionMajorantSpeed, double ionizationPotentialEv)
    : CollisionReaction(population),
      ionMeanSpeed_(ionMeanSpeed), ionMajorantSpeed_(ionMajorantSpeed),
      ionizationPotentialEv_(ionizationPotentialEv) {
    chargeState_ = population_.getChargeState();
    mullerSalzbornCrossSection_ = 1.43e-16 *
        std::pow((double)chargeState_, 1.17) *
        std::pow(ionizationPotentialEv_, -2.76);
}

double ChargeExchangeReaction::getCrossSection(double) const {
    return mullerSalzbornCrossSection_;
}

double ChargeExchangeReaction::getMeanReactionRate(const Point &p, double) const {
    return getCrossSection(ionMeanSpeed_) * ionMeanSpeed_ *
        population_.getDensityAt(p);
}

double ChargeExchangeReaction::getMajorantReactionRate(const Point &p,
    double relativeSpeed) const {
    return getCrossSection(ionMajorantSpeed_) * ionMajorantSpeed_ *
        population_.getDensityAt(p);
}

std::vector<Particle> ChargeExchangeReaction::computeReactionProducts(
    Rng &rng, const Point &, const Particle &target) const {
    // Elastic collision kinematics calculated in center of mass frame

    // Projectile particle velocity is isotropic Maxwell-Boltzmann:
    Vector ionVelocity = population_.getRandomParticleVelocity(rng);
    Vector neutralVelocity = target.getVelocity();
    double ionMass = population_.getParticleMass_eV();
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

