#include "integration/gsl_integration.h"
#include "monte/gsl_monte.h"
#include "monte/gsl_monte_vegas.h"
#include "rng/gsl_rng.h"
#include "util.h"

void Util::printPoint(const Point &p) {
    std::cout << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
}

void Util::printVector(const Vector &v) {
    std::cout << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
}

int Util::fastFloor(double x) {
    int i = x;
    return i - (i > x);
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
    double (*crossSectionFn)(double, void *), void *params,
    gsl_integration_workspace *ws) {
    double mbMean = getMBAverage(T_eV, mass_eV);
    gsl_function fn;
    mbparams par;
    par.T_eV = T_eV;
    par.mass_eV = mass_eV;
    par.fnParams = params;
    par.fn = crossSectionFn;
    fn.function = mbHelperFn;
    fn.params = &par;
    double result, error;
    gsl_integration_qag(&fn, 0.0, 8.0*mbMean, 1e-12, 1e-6,
        RATE_COEFF_WORKSPACE_SIZE, 6, ws, &result, &error);
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

double Util::relativeSpeedHelper(double *v, size_t, void *params) {
    mbrelativeparams *par = (mbrelativeparams *)params;
    double mper2kT = par->mass_eV /
        (2.0 * par->T_eV * SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    double vmean2 = par->vmean * par->vmean;
    double vx2 = v[0]*v[0];
    double vy2 = v[1]*v[1];
    double vz2 = v[2]*v[2];
    double vzdiff2 = v[2] - par->particleSpeed;
    vzdiff2 *= vzdiff2;
    double vsq = (vx2 + vy2 + vz2) * vmean2;
    double vdiffsq = (vx2 + vy2 + vzdiff2) * vmean2;
    double vdiff = std::sqrt(vdiffsq);
    return std::pow(mper2kT / M_PI, 1.5) * vdiff *
        par->fn(vdiff, par->fnParams) * std::exp(-mper2kT * vsq);
}

double Util::calculateMBRelativeRateCoeff(double particleSpeed, double T_eV,
    double mass_eV, gsl_rng *rng, monte_state *ms,
    double (*fn)(double, void *), void *fnParams, size_t N_calls) {
    double vmean = getMBAverage(T_eV, mass_eV);
    mbrelativeparams par;
    par.T_eV = T_eV;
    par.mass_eV = mass_eV;
    par.particleSpeed = particleSpeed / vmean;
    par.vmean = vmean;
    par.fn = fn;
    par.fnParams = fnParams;

    gsl_monte_function mfn;
    mfn.f = relativeSpeedHelper;
    mfn.dim = 3;
    mfn.params = &par;

    double lim = 8.0;
    double xl[] = { 0.0, 0.0, -lim };
    double xu[] = { lim, lim, lim };
    double result, absErr;
    gsl_monte_vegas_init(ms);
    gsl_monte_vegas_integrate(&mfn, xl, xu, 3, N_calls, rng, ms, &result,
        &absErr);

    return 4.0*vmean*vmean*vmean*result;
}

double Util::calculateMBRelativeSpeed(double particleSpeed, double T_eV,
    double mass_eV, gsl_rng *rng, monte_state *ms, size_t N_calls) {
    return calculateMBRelativeRateCoeff(particleSpeed, T_eV, mass_eV, rng,
        ms, fOne, NULL, N_calls);
}

simthreadresources *Util::allocateThreadResources(uint_least32_t seed) {
    simthreadresources *thread_res = new simthreadresources;

    thread_res->rng = Rng(seed);
    thread_res->ms = gsl_monte_vegas_alloc(3);
    thread_res->ws = gsl_integration_workspace_alloc(RATE_COEFF_WORKSPACE_SIZE);
    thread_res->gslrng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(thread_res->gslrng, seed);

    return thread_res;
}

void Util::deallocateThreadResources(simthreadresources *thread_res) {
    if (thread_res == NULL) return;
    gsl_monte_vegas_free(thread_res->ms);
    gsl_integration_workspace_free(thread_res->ws);
    gsl_rng_free(thread_res->gslrng);
    delete thread_res;
}
