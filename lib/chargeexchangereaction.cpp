#include "chargeexchangereaction.h"

ChargeExchangeReaction::ChargeExchangeReaction(
    const std::shared_ptr<ParticlePopulation> &population,
    double ionizationPotentialEv)
    : CollisionReaction(population),
      ionizationPotentialEv_(ionizationPotentialEv) {
    chargeState_ = population_->getChargeState();
    label_ = "charge exchange";
    mullerSalzbornCrossSection_ = 1.43e-16 *
        std::pow(chargeState_, 1.17) *
        std::pow(ionizationPotentialEv_, -2.76);
}

double ChargeExchangeReaction::getCrossSection(double) const {
    return mullerSalzbornCrossSection_;
}

double ChargeExchangeReaction::getReactionRate(const Point &p,
    double particleSpeed, mc_integrate_resources &mc_res) const {
    return population_->getDensityAt(p) *
        getRateCoefficient(particleSpeed, mc_res);
}

double ChargeExchangeReaction::getRateCoefficient(double particleSpeed,
    mc_integrate_resources &mc_res) const {
    double tmpCrossSection = mullerSalzbornCrossSection_;
    return population_->calculateRateCoefficient(particleSpeed, mc_res.ms,
        mc_res.gslrng, crossSection, static_cast<void *>(&tmpCrossSection));
}

double ChargeExchangeReaction::crossSection(double, void *args) {
    return *(static_cast<double *>(args));
}

CollisionProducts ChargeExchangeReaction::computeReactionProducts(
    Rng &rng, const Point &, const Particle &target) const {
    std::vector<Particle> products;
    if (chargeState_ == 1) {
        // Elastic collision kinematics calculated in center of mass frame
        // Projectile particle velocity is isotropic Maxwell-Boltzmann:
        Vector ionVelocity = population_->getRandomParticleVelocity(rng);
        Vector neutralVelocity = target.getVelocity();
        double ionMass = population_->getParticleMass_eV();
        double neutralMass = target.getMass_eV();

        Vector vcm = (ionMass * ionVelocity + neutralMass * neutralVelocity) /
            (ionMass + neutralMass);

        Vector neutralVelocityCm = neutralVelocity - vcm;

        double momentumNorm = std::sqrt(neutralVelocityCm.squared_length() /
            neutralMass);

        Direction dir = Util::getIsotropicSphereDirection(rng);
        Vector productVel = vcm + dir.vector() * momentumNorm / ionMass;

        Particle neutralProduct(population_->getElement(), target.getTime());
        neutralProduct.setVelocity(productVel);

        products.push_back(neutralProduct);

        return std::make_pair(products, 0);
    }

    // If the bombarding population is greater than 1+, there are no neutral
    // products but the neutral is ionized and plasma charge state is offset by
    // one
    return std::make_pair(products, 1);
}

