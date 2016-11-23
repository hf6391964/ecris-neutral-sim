#include "integration/gsl_integration.h"
#include "util.h"

void Util::printPoint(const Point &p) {
    std::cout << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
}

void Util::printVector(const Vector &v) {
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

double Util::mbHelperFn(double speed, void *p) {
    mbparams *params = (mbparams *)p;

    return evaluateMBDistribution(params->T_eV, params->mass_eV, speed) *
        speed * (*(params->fn))(speed, params->fnParams);
}

double Util::evaluateMBDistribution(double T_eV, double mass_eV,
    double speed) {
    double mper2kT = mass_eV / (2.0 * T_eV * SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    double vsq = speed*speed;
    return 4.0 * std::sqrt(mper2kT*mper2kT*mper2kT / M_PI) * vsq *
        std::exp(-mper2kT * vsq);
}

double Util::calculateMBRateCoefficient(double T_eV, double mass_eV,
    double (*crossSectionFn)(double, void *), void *params) {
    double mbMean = getMBAverage(T_eV, mass_eV);
    gsl_integration_workspace *ws = gsl_integration_workspace_alloc(1000);
    gsl_function fn;
    mbparams par;
    par.T_eV = T_eV;
    par.mass_eV = mass_eV;
    par.fnParams = params;
    par.fn = crossSectionFn;
    fn.function = mbHelperFn;
    fn.params = &par;
    double result, error;
    gsl_integration_qag(&fn, 0, 8.0*mbMean, 1e-50, 1e-12, 1000, 6, ws, &result,
        &error);
    gsl_integration_workspace_free(ws);
    return result;
}

double Util::getMBAverage(double T_eV, double mass_eV) {
    return std::sqrt(8.0 / M_PI * (T_eV / mass_eV)) * SPEED_OF_LIGHT;
}

Vector Util::getMBVelocity(Rng& rng, double T_eV, double mass_eV) {
    const double a = std::sqrt(T_eV / mass_eV) * SPEED_OF_LIGHT;
    std::normal_distribution<double> d(0.0, a);
    return Vector(d(rng), d(rng), d(rng));
}

Direction Util::getIsotropicSphereDirection(Rng &rng) {
    double u = 2.0 * uni01(rng) - 1.0;
    double v = std::sqrt(1 - u*u);
    double theta = 2.0 * M_PI * uni01(rng);
    return Direction(v * std::cos(theta), v * std::sin(theta), u);
}

