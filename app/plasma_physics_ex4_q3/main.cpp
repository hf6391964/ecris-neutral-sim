#include <cmath>
#include <fstream>
#include <iostream>

#include "util.h"

double sigma_exc(double speed, void *) {
    const double E = ELECTRON_MASS_EV * speed * speed / (2.0 * SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    const double DELTA_E = 12.754,
                 A1 = 3.651, A2 = -0.8405, A3 = 1.2365, A4 = 2.5236;
    if (E < DELTA_E) return 0.0;
    const double x = E / DELTA_E;
    const double x_recip = 1.0 / x;

    double val = 5.984 / (DELTA_E * x) * std::pow(1.0 - x_recip, 0.55) *
        (A1 + x_recip * (A2 + A3 * x_recip) + A4 * std::log(x));
    return val;
}

double sigma_ion(double speed, void *) {
    const double E = ELECTRON_MASS_EV * speed * speed / (2.0 * SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    const double DELTA_E = 15.42;
    if (E < DELTA_E) return 0.0;
    const double C1 = 2.05 * DELTA_E;
    const double x = E / DELTA_E;

    double val = 1.828 / x * std::pow(1.0 - std::pow(x, -0.92), 2.19) *
        std::log(C1 * x);
    return val;
}

int main() {
    const double EXP_MIN = 0.0, EXP_MAX = 3.0;
    const int N_INTERVALS = 100;
    std::ofstream fexc("exc.csv");
    std::ofstream fion("ion.csv");
    for (int i = 0; i <= N_INTERVALS; i++) {
        double exponent = EXP_MIN +
            (EXP_MAX - EXP_MIN) * (double)i / N_INTERVALS;
        double T = std::pow(10.0, exponent);
        double rc_exc = Util::calculateMBRateCoefficient(T, ELECTRON_MASS_EV,
            sigma_exc, NULL);
        double rc_ion = Util::calculateMBRateCoefficient(T, ELECTRON_MASS_EV,
            sigma_ion, NULL);
        fexc << T << ',' << rc_exc << '\n';
        fion << T << ',' << rc_ion << '\n';
    }
    fexc.close();
    fion.close();
}

