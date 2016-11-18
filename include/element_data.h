#pragma once

// Ionization paramaters as per
// https://www.nist.gov/sites/default/files/documents/srd/jpcrd347.pdf
struct IonizationParameters {
    double I;
    double A;
    double B[5];
};

struct ElementData {
    unsigned int Z;  // proton number
    double mass;  // in atomic mass units

    // If this is 0, low energy parameters are always used
    double ionizationHighEnergyLimit;
    IonizationParameters ionizationParametersLowEnergy;
    // Parameters for impact energies larger than the high energy limit
    IonizationParameters ionizationParametersHighEnergy;
};

