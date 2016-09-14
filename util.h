#ifndef UTIL_H
#define UTIL_H

#include "cgal_and_typedefs.h"
#include "constants.h"

class Util {
    public:
        static void printPoint(Point p);
        static void printVector(Vector v);

        static double getMBSpeed(Rng& rng, double T, double molarMass);
        static double getMBAverage(double T, double molarMass);
};

#endif
