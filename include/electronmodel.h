#ifndef ELECTRONMODEL_H
#define ELECTRONMODEL_H

#include <cmath>

#include "cgal_and_typedefs.h"

class ElectronModel {
    static constexpr double Ai[7] = { 0.341949, -0.102496, 57.3281, 20.3848,
        -1175.09, -4300.16, 23335.8 };

    private:
        Vector Bfield(const Vector vx, const double B0, const double r0,
            const double a[] = Ai) const;
};

#endif

