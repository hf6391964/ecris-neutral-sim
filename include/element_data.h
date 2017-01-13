#pragma once

#include <unordered_map>

// Ionization paramaters as per
// https://www.nist.gov/sites/default/files/documents/srd/jpcrd347.pdf
struct IonizationParameters {
    double I_le;  // Ionization potential, units eV
    // Units chosen so that the coefficients can be read directly from Table 1
    // in above source.
    double A_le;  // Parameter a in units 1e-17 eV^2 m^2 = 1e-13 eV^2 cm^2
    double B_le[5];  // Parameters Bi in units 1e-17 eV^2 m^2

    // If this is 0, low energy parameters are always used
    double ionizationEnergyLimit;

    // High energy (E > ionizationEnergyLimit) parameters
    double I_he;
    double A_he;
    double B_he[5];
};

struct ElementData {
    unsigned int Z;  // proton number
    double mass;  // in atomic mass units

    IonizationParameters ionizationParameters;
};

enum Element { ARGON };

extern const ElementData ARGON_DATA;

extern const std::unordered_map<Element, const ElementData *> ELEMENT_DATA;

