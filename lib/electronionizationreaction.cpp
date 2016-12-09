#include "integration/gsl_integration.h"
#include "electronionizationreaction.h"

ElectronIonizationReaction::ElectronIonizationReaction(
    ParticlePopulation &population, double electronMeanSpeed,
    double electronMajorantSpeed, const IonizationParameters &ip)
    : electronMeanSpeed_(electronMeanSpeed),
      electronMajorantSpeed_(electronMajorantSpeed),
      population_(population), ionizationParameters_(ip) {
}

double ElectronIonizationReaction::getMeanReactionRate(const Point &p,
    double relativeSpeed) const {
    // TODO return precomputed values from table
    /* return plasmaDensities_.getElectronDensityAt(p) * */
    /*     Util::calculateMBRateCoefficient( */
    /*         plasmaDensities_.getElectronTemperature(), */
    /*         ELECTRON_MASS_EV, crossSection, (void *)&ionizationParameters_, ws); */
    return 0.0;
}

double ElectronIonizationReaction::getMajorantReactionRate(const Point &p,
    double relativeSpeed) const {
    /* return plasmaDensities_.getElectronDensityAt(p) * */
    /*     getCrossSection(p, electronMajorantSpeed_) * electronMajorantSpeed_; */
    // TODO calculate based on highest density and precomputed rate coeffs
    return 0.0;
}

std::vector<Particle> ElectronIonizationReaction::computeReactionProducts(
    Rng &, const Point &, const Particle &) const {
    return std::vector<Particle>();
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
    return (aterm + bterms) / (I*E);
}

