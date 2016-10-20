#include <cmath>

#include "lotz.h"

double ionizationCrossSection(unsigned int N, double electronEnergy,
    LotzParameters parameters, double ionizationEnergies[],
    unsigned int electronPopulation[]) {
    double crossSection = 0.0;

    for (unsigned int i = 0; i < N; i++) {
        double EoverPi = electronEnergy / ionizationEnergies[i];
        crossSection += parameters.ai[i] * (double)electronPopulation[i] *
            std::log(EoverPi) / (electronEnergy * ionizationEnergies[i]) *
            (1.0 - parameters.bi[i] * std::exp(-parameters.ci[i] *
                                               (EoverPi - 1.0)));
    }

    return crossSection;
}
