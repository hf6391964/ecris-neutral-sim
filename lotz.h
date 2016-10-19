#ifndef LOTZ_H
#define LOTZ_H

#include "element_data.h"

struct LotzParameters {
    double ai[MAX_SHELLS];
    double bi[MAX_SHELLS];
    double ci[MAX_SHELLS];
};

double ionizationCrossSection(unsigned int N, double electronEnergy,
    LotzParameters parameters, double ionizationEnergies[],
    unsigned int electronPopulation[]);

#endif

