#include <iostream>

#include "cgal_and_typedefs.h"
#include "plasmadensities.h"

int main() {
    std::vector<double> ionDensities;
    PlasmaDensities densities("electron_density_hiisi.csv", 1.0, ionDensities);

    return 0;
}

