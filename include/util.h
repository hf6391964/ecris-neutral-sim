#pragma once

#include "cgal_and_typedefs.h"
#include "constants.h"

#include "rng/gsl_rng.h"
#include "integration/gsl_integration.h"
#include "monte/gsl_monte_vegas.h"

#define RATE_COEFF_WORKSPACE_SIZE 1000

struct mbparams {
    double mass_eV;
    double T_eV;
    void *fnParams;
    double (*fn)(double, void *);
};

struct mbrelativeparams {
    double mass_eV;
    double T_eV;
    double particleSpeed;
    double vmean;
    void *fnParams;
    double (*fn)(double, void *);
};

struct simthreadresources {
    Rng rng;
    monte_state *ms;
    gsl_integration_workspace *ws;
    gsl_rng *gslrng;

    simthreadresources(uint_least32_t seed) {
        rng = Rng(seed);
        ms = gsl_monte_vegas_alloc(3);
        ws = gsl_integration_workspace_alloc(RATE_COEFF_WORKSPACE_SIZE);
        gslrng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(gslrng, seed);
    }

    ~simthreadresources() {
        gsl_monte_vegas_free(ms);
        gsl_integration_workspace_free(ws);
        gsl_rng_free(gslrng);
    }
};

class Util {
    private:
        static double mbHelperFn(double speed, void *params);
        static double relativeSpeedHelper(double *v, size_t dim, void *params);

        static double fOne(double, void *) { return 1.0; }

    public:
        static void printPoint(const Point &p);
        static void printVector(const Vector &v);

        static double getMBSpeed(Rng &rng, double T_eV, double mass_eV);
        static double getMBAverage(double T, double mass_eV);
        static Vector getMBVelocity(Rng &rng, double T_eV, double mass_eV);
        static double evaluateMBDistribution(double T_eV, double mass_eV,
            double speed);
        static double calculateMBRateCoefficient(double T_eV, double mass_eV,
            double (*crossSectionFn)(double, void *), void *params,
            gsl_integration_workspace *ws);

        static Direction getIsotropicSphereDirection(Rng &rng);

        static int fastFloor(double x);

        static double calculateMBRelativeRateCoeff(double particleSpeed,
            double T_eV, double mass_eV, gsl_rng *rng, monte_state *ms,
            double (*fn)(double, void *), void *fnParams,
            size_t N_calls = 10000);

        static double calculateMBRelativeSpeed(double particleSpeed,
            double T_eV, double mass_eV, gsl_rng *rng, monte_state *ms,
            size_t N_calls = 10000);
};
