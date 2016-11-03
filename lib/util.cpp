#include "util.h"

void Util::printPoint(Point p) {
    std::cout << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
}

void Util::printVector(Vector v) {
    std::cout << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
}

double Util::getMBSpeed(Rng& rng, double T_eV, double mass_eV) {
    // Maxwell-Boltzmann distribution implemented by drawing the three velocity
    // components from normal distribution and taking the resulting magnitude
    const int N_DIM = 3;
    const double a = std::sqrt(T_eV / mass_eV) * SPEED_OF_LIGHT;
    std::normal_distribution<double> d(0.0, a);

    double v = 0.0;
    for (int i = 0; i < N_DIM; i++) {
        double tmp = d(rng);
        v += tmp*tmp;
    }
    return std::sqrt(v);
}

double Util::getMBAverage(double T_eV, double mass_eV) {
    return std::sqrt(8.0 / M_PI * (T_eV / mass_eV)) * SPEED_OF_LIGHT;
}

