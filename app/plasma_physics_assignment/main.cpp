#include <iostream>
#include <fstream>
#include <string>

#include "cgal_and_typedefs.h"
#include "electronmodel.h"

int main() {
    const double z1 = -160e-3,
                 z2 = 170e-3,
                 r0 = 39e-3,
                 B0 = 1.07,
                 dt = 1e-9;
    // B0, r0, dt, z1, z2
    ElectronModel model(B0, r0, dt, z1, z2);

    // Plot solenoid Bz(r=0, z)
    int nPoints = 100;
    std::ofstream fout("output/solenoid_B_onaxis.csv");
    for (int i = 0; i <= nPoints; i++) {
        double z = z1 + (z2 - z1) * (double)i/nPoints;
        Vector B = ElectronModel::solenoidBfield(0.0, 0.0, z);
        fout << z << CSV_SEP << B.z() << std::endl;
    }
    fout.close();

    // Plot hexapole Bxy(x, y)
    fout.open("output/hexapole_B.csv");
    nPoints = 15;
    for (int i = 0; i <= nPoints; i++) {
        double y = -r0 + 2.0*r0 * (double)i/nPoints;
        for (int j = 0; j <= nPoints; j++) {
            double x = -r0 + 2.0*r0 * (double)j/nPoints;
            Vector B = ElectronModel::hexapoleBfield(x, y, B0, r0);
            fout << x << CSV_SEP << y << CSV_SEP << B.x() << CSV_SEP <<
                B.y() << std::endl;
        }
    }
    fout.close();

    return 0;
}

