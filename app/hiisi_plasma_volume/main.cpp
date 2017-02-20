#include <iostream>

#include "gsl/gsl_odeiv2.h"
#include "cgal_and_typedefs.h"
#include "electronmodel.h"


Vector arr_to_Vector(const double x[]) {
    return Vector(x[0], x[1], x[2]);
}

struct Params {
    double B0;
    double r0;
    std::vector<double> a;
};

int func(double t, const double y[], double f[], void *par) {
    Vector x = arr_to_Vector(y);
    const Params *params = (const Params *)par;
    Vector B = ElectronModel::totalBfield(x, params->B0, params->r0, params->a);
    double Bnorm = std::sqrt(B.squared_length());
    double c = -1.0 / Bnorm;
    f[0] = c*B.x();
    f[1] = c*B.y();
    f[2] = c*B.z();
    return GSL_SUCCESS;
}

int main() {
    const double z1 = -203.8e-3,
                 z2 = 196.2e-3,
                 r0 = 55e-3,
                 B0 = 1.188,
                 gridSize = 0.001,
                 extractionRadius = 0.004,
                 stepLength = 0.5 * gridSize;
    const std::vector<double> Ai
        { 0.413286, 0.954437, 52.1123, -74.7674, -1258.17, 2994.75, 22469.9 };

    Params params {
        .B0 = B0,
        .r0 = r0,
        .a = Ai
    };
    const unsigned long N_PARTICLES = 1000000;

    Rng rng(20022017);

    Grid grid(Bbox(-r0, -r0, z1, r0, r0, z2), gridSize);
    const double cellVolume = gridSize*gridSize*gridSize;
    size_t arraySize = grid.arraySize();
    std::vector<bool> visited(arraySize, false);

    gsl_odeiv2_system sys = {func, NULL, 3, &params};

    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new (&sys,
        gsl_odeiv2_step_msadams, 1e-6, 1e-9, 0.0);

    for (unsigned long i = 0; i < N_PARTICLES; ++i) {
        // Uniformly sample a position in the extraction hole
        double sqrtr = extractionRadius * std::sqrt(uni01(rng)),
               theta = 2.0 * M_PI * uni01(rng);
        double position[3] = {
            sqrtr * std::cos(theta),
            sqrtr * std::sin(theta),
            z2
        };
        double t = 0.0;

        gsl_odeiv2_driver_reset(driver);

        do {
            // Trace the fieldline and mark the visited cells
            if (gsl_odeiv2_driver_apply(driver, &t, t + stepLength,
                position) != GSL_SUCCESS) {
                throw std::runtime_error("Integration failed");
            }

            size_t idx;
            if (!grid.arrayIndex(position[0], position[1], position[2], idx)) {
                break;
            }
            visited[idx] = true;
        } while (position[2] > z1 && position[2] <= z2 &&
            position[0]*position[0] + position[1]*position[1] < r0*r0);
    }

    gsl_odeiv2_driver_free(driver);

    unsigned long sum = 0;
    for (bool cell_visited : visited) {
        if (cell_visited) {
            sum += 1;
        }
    }

    double volume = (double)sum * cellVolume;
    double chamberVolume = (z2 - z1) * M_PI * r0*r0;
    std::cout << "Plasma volume = " << volume << " m^3 = " <<
        (volume * 1e9) << " mm^3\n";
    std::cout << (100.0 * volume / chamberVolume) << " % of chamber volume\n";

    std::ofstream f("plasma_volume.txt");
    f << volume << '\n';

    return EXIT_SUCCESS;
}

