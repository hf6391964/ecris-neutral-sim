#ifndef ELEMENT_DATA_H
#define ELEMENT_DATA_H

#define MAX_Z 60

struct ElementData {
    unsigned int Z;  // proton number
    double mass;  // in atomic mass units
    double ionizationEnergies[MAX_Z];
};

#endif

