#include <iostream>

#include "cgal_and_typedefs.h"
#include "electronmodel.h"

int main() {
    const double z1 = -203.8e-3,
                 z2 = 196.2e-3,
                 r0 = 55e-3,
                 B0 = 1.188,
                 dt = 5e-12,
                 confinementTime = 1e-6,
                 gridSize = 0.005,
                 Becr = 0.643;
    const double Ai[] = { 0.413286, 0.954437, 52.1123, -74.7674, -1258.17, 2994.75, 22469.9 };
    const unsigned long N_PARTICLES = 100000;

    // B0, r0, dt, z1, z2
    ElectronModel model(B0, r0, dt, z1, z2, gridSize, confinementTime, Ai);

    Rng rng(2016);
    model.runMonoenergeticSimulation(N_PARTICLES, 10e3, Becr, rng);

    model.writeDensityToFile("electron_density_hiisi.csv");
    model.writeElectronEndpoints("z1_collision_points.csv",
        "z2_collision_points.csv", "cylinder_collision_points.csv");

    return 0;
}

