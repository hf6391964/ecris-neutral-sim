#include "integration/gsl_integration.h"
#include "electronionizationreaction.h"

ElectronIonizationReaction::ElectronIonizationReaction(
    const ParticlePopulation &population,
    const IonizationParameters &ip)
    : CollisionReaction(population),
      ionizationParameters_(ip) {
    label_ = "ionization";
}

double ElectronIonizationReaction::getReactionRate(const Point &p,
    double particleSpeed, simthreadresources *thread_res) const {
    return population_.getDensityAt(p) * getRateCoefficient(particleSpeed,
        thread_res);
}

double ElectronIonizationReaction::getRateCoefficient(
    double particleSpeed, simthreadresources *thread_res) const {
    return population_.calculateRateCoefficient(particleSpeed, thread_res->ms,
        thread_res->gslrng, crossSection, (void *)&ionizationParameters_);
}

CollisionProducts ElectronIonizationReaction::computeReactionProducts(
    Rng &, const Point &, const Particle &) const {
    return std::make_pair(std::vector<Particle>(), 1);
}

double ElectronIonizationReaction::getCrossSection(double speed) const {
    return crossSection(speed, (void *)&ionizationParameters_);
}

double ElectronIonizationReaction::crossSection(double v, void *p) {
    // Calculate ionization energy in eV
    double E = 0.5 * ELECTRON_MASS_EV * v*v / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);

    double A, I;
    double *B;
    IonizationParameters *ip = (IonizationParameters *)p;
    if (E < ip->I_le) return 0.0;  // Can't ionize if impact energy lower than
                                   // ionization energy
    if (E < ip->ionizationEnergyLimit * ip->I_le) {
        A = ip->A_le;
        I = ip->I_le;
        B = ip->B_le;
    } else {
        A = ip->A_he;
        I = ip->I_he;
        B = ip->B_he;
    }

    double aterm = A * std::log(E / I);
    // Modified Horner scheme
    double delta = 1.0 - I / E;
    double bterms = 0.0;
    for (int i = 4; i >= 0; --i) {
        bterms = B[i] + bterms*delta;
    }
    bterms *= delta;
    return (aterm + bterms) / (I*E) * 1e-17;
}

