#include "chargeexchangereaction.h"

ChargeExchangeReaction::ChargeExchangeReaction(
    std::shared_ptr<ParticlePopulation> population,
    double ionizationPotentialEv)
    : CollisionReaction(population),
      ionizationPotentialEv_(ionizationPotentialEv) {
    label_ = "charge exchange";
    chargeState_ = population_->getChargeState();
    mullerSalzbornCrossSection_ = 1.43e-16 *
        std::pow((double)chargeState_, 1.17) *
        std::pow(ionizationPotentialEv_, -2.76);
}

double ChargeExchangeReaction::getCrossSection(double) const {
    return mullerSalzbornCrossSection_;
}

double ChargeExchangeReaction::getReactionRate(const Point &p,
    double particleSpeed, simthreadresources *thread_res) const {
    return population_->getDensityAt(p) *
        getRateCoefficient(particleSpeed, thread_res);
}

double ChargeExchangeReaction::getRateCoefficient(double particleSpeed,
    simthreadresources *thread_res) const {
    return population_->calculateRateCoefficient(particleSpeed, thread_res->ms,
        thread_res->gslrng, crossSection, (void *)&mullerSalzbornCrossSection_);
}

double ChargeExchangeReaction::crossSection(double, void *args) {
    return *((double *)args);
}

CollisionProducts ChargeExchangeReaction::computeReactionProducts(
    Rng &rng, const Point &, const Particle &target) const {
    // Elastic collision kinematics calculated in center of mass frame
    std::vector<Particle> products;

    if (chargeState_ == 1) {
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
    }

    return std::make_pair(products, 0);
}

