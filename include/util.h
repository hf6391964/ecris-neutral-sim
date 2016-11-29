#pragma once

#include "cgal_and_typedefs.h"
#include "constants.h"

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
};

class Util {
    private:
        static double mbHelperFn(double speed, void *params);
        static double relativeSpeedHelper(double *v, size_t dim, void *params);

    public:
        static void printPoint(const Point &p);
        static void printVector(const Vector &v);

        static double getMBSpeed(Rng &rng, double T_eV, double mass_eV);
        static double getMBAverage(double T, double mass_eV);
        static Vector getMBVelocity(Rng &rng, double T_eV, double mass_eV);
        static double evaluateMBDistribution(double T_eV, double mass_eV,
            double speed);
        static double calculateMBRateCoefficient(double T_eV, double mass_eV,
            double (*crossSectionFn)(double, void *), void *params);

        static Direction getIsotropicSphereDirection(Rng &rng);

        static int fastFloor(double x) {
            int i = x;
            return i - (i > x);
        }

        static double calculateMBRelativeSpeed(double particleSpeed,
            double T_eV, double mass_eV, size_t N_calls = 10000);
};

