#include <iostream>
#include <fstream>

#include "constants.h"
#include "util.h"

int main() {
    double vmean = Util::getMBAverage(ROOM_TEMPERATURE_EV, ELECTRON_MASS_EV);
    // Test convergence
    double expmin = 1.0, expmax = 7.0;
    size_t n_points = 20;
    std::ofstream fout("data.csv");
    for (size_t i = 0; i < n_points; ++i) {
        size_t samplesize =
            std::pow(10.0, expmin + (expmax - expmin) * (double)i/n_points);
        fout << samplesize << "," << Util::calculateMBRelativeSpeed(vmean, ROOM_TEMPERATURE_EV,
            ELECTRON_MASS_EV, samplesize) << std::endl;
    }
    fout.close();

    // Test correspondence of results with mean
    double speeds[] = { 0.0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0 };
    std::cout << "vmean: " << vmean << std::endl;

    for (double speed : speeds) {
        double v = speed * vmean;
        double vrel = Util::calculateMBRelativeSpeed(v, ROOM_TEMPERATURE_EV,
            ELECTRON_MASS_EV, 100000);
        std::cout << "v: " << v << ", vrel: " << vrel << std::endl;
    }
}

