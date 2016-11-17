#pragma once

#include "cgal_and_typedefs.h"
#include "constants.h"

class Util {
    public:
        static void printPoint(const Point &p);
        static void printVector(const Vector &v);

        static double getMBSpeed(Rng &rng, double T_eV, double mass_eV);
        static double getMBAverage(double T, double mass_eV);
        static Vector getMBVelocity(Rng &rng, double T_eV, double mass_eV);

        static Direction getIsotropicSphereDirection(Rng &rng);

        static int fastFloor(double x) {
            int i = x;
            return i - (i > x);
        }
};

