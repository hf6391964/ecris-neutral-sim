#pragma once

#define MAX_SHELLS 20

struct ElementData {
    unsigned int Z;  // proton number
    double mass;  // in atomic mass units
    // Shell ordering: TODO fill
    // 1s, 2s, 2p, 3s, 3p, ...
    double ionizationEnergies[MAX_SHELLS];
    unsigned int electronPopulation[MAX_SHELLS];
};

