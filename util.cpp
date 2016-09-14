#include "util.h"

void Util::printPoint(Point p) {
    std::cout << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
}

void Util::printVector(Vector v) {
    std::cout << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
}

double Util::getMBSpeed(Rng& rng, double T, double molarMass) {
    // Maxwell-Boltzmann distribution implemented by drawing the three velocity
    // components from normal distribution and taking the resulting magnitude
    const int N_DIM = 3;
    const double a = std::sqrt(GAS_CONSTANT * (T / molarMass));
    std::normal_distribution<double> d(0.0, a);

    double v = 0.0;
    for (int i = 0; i < N_DIM; i++) {
        double tmp = d(rng.engine());
        v += tmp*tmp;
    }
    return std::sqrt(v);
}

double Util::getMBAverage(double T, double molarMass) {
    return std::sqrt(8.0 * GAS_CONSTANT / M_PI * (T / molarMass));
}

