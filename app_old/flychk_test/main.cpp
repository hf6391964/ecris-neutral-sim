#include <iostream>
#include <fstream>
#include <cmath>

#include "flychkparser.h"

int main() {
    FlychkParser parser("../../data/rt.018.dat");
    double minexp = 0.0, maxexp = 6.0;
    size_t nSteps = 300;

    std::ofstream datafile;

    datafile.open("rr.csv");
    for (size_t i = 0; i < nSteps; ++i) {
        double exponent = minexp + i * (maxexp - minexp) / nSteps;
        double T = std::pow(10, exponent);

        double rcoeff_lin, rcoeff_log;
        if (parser.getRateCoefficient(T, 1, "rr", rcoeff_lin, false) &&
            parser.getRateCoefficient(T, 1, "rr", rcoeff_log, true)) {
            datafile << T << ' ' << rcoeff_lin << ' ' << rcoeff_log << std::endl;
        }
    }
    datafile.close();

    datafile.open("dr.csv");
    for (size_t i = 0; i < nSteps; ++i) {
        double exponent = minexp + i * (maxexp - minexp) / nSteps;
        double T = std::pow(10, exponent);

        double rcoeff_lin, rcoeff_log;
        if (parser.getRateCoefficient(T, 1, "dr", rcoeff_lin, false) &&
            parser.getRateCoefficient(T, 1, "dr", rcoeff_log, true)) {
            datafile << T << ' ' << rcoeff_lin << ' ' << rcoeff_log << std::endl;
        }
    }
    datafile.close();

    datafile.open("total.csv");
    for (size_t i = 0; i < nSteps; ++i) {
        double exponent = minexp + i * (maxexp - minexp) / nSteps;
        double T = std::pow(10, exponent);

        double rcoeff_lin, rcoeff_log;
        if (parser.getTotalRateCoefficient(T, 1, rcoeff_lin, false) &&
            parser.getTotalRateCoefficient(T, 1, rcoeff_log, true)) {
            datafile << T << ' ' << rcoeff_lin << ' ' << rcoeff_log << std::endl;
        }
    }
    datafile.close();
}

