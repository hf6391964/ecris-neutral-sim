#include "integration/gsl_integration.h"
#include "electronionizationreaction.h"

ElectronIonizationReaction::ElectronIonizationReaction(
    PlasmaDensities &densities, double electronMeanSpeed,
    double electronMajorantSpeed, const IonizationParameters &ip)
    : electronMeanSpeed_(electronMeanSpeed),
      electronMajorantSpeed_(electronMajorantSpeed),
      plasmaDensities_(densities), ionizationParameters_(ip) {
}

double ElectronIonizationReaction::getMeanReactionRate(const Point &p) const {
    // TODO rate coefficient calculation
    return plasmaDensities_.getElectronDensityAt(p) *
        Util::calculateMBRateCoefficient(
            plasmaDensities_.getElectronTemperature(),
            ELECTRON_MASS_EV, crossSection, (void *)&ionizationParameters_);
}

double ElectronIonizationReaction::getMajorantReactionRate(const Point &p)
    const {
    return plasmaDensities_.getElectronDensityAt(p) *
        getCrossSection(p, electronMajorantSpeed_) * electronMajorantSpeed_;
}

std::vector<Particle> ElectronIonizationReaction::computeReactionProducts(
    Rng &, const Point &, const Particle &) const {
    return std::vector<Particle>();
}

double ElectronIonizationReaction::getCrossSection(const Point &,
    double speed) const {
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

