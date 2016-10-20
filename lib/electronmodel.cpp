#include "electronmodel.h"

Vector ElectronModel::solenoidBfield(const double x, const double y,
    const double z, const double a[]) {
    double r2 = x*x + y*y;
    // Horner scheme w.r.t. z
    double solenoid_Bz =
           a[0] + r2*(0.5*a[2] + r2*((3.0/8.0)*a[4] - r2*(5.0/16.0)*a[6])) +
        z*(a[1] + r2*(-1.5*a[3] + r2*(15.0/8.0)*a[5]) +
        z*(a[2] + r2*(-3.0*a[4] + r2*(45.0/8.0)*a[6]) +
        z*(a[3] - r2*5.0*a[5] +
        z*(a[4] - r2*(15.0/2.0)*a[6] +
        z*(a[5] + z*a[6])))));
    double solenoid_Br_divided_by_r =
           -0.5*a[1] + r2*((3.0/8.0) - r2*(5.0/16.0)*a[5]) +
        z*(-a[2] + r2*(1.5*a[4] - r2*(15.0/8.0)*a[6]) +
        z*(-1.5*a[3] + r2*(15.0/4.0)*a[5] +
        z*(-2.0*a[4] + r2*(15.0/2.0)*a[6] +
        z*(-2.5*a[5] + z*(-3.0*a[6])))));

    return Vector(
        solenoid_Br_divided_by_r * x,
        solenoid_Br_divided_by_r * y,
        solenoid_Bz
    );
}

Vector ElectronModel::hexapoleBfield(const double x, const double y,
    const double B0, const double r0) {
    double C = -2.0*B0 / (r0*r0);
    return Vector(C*x*y, C*(x*x - y*y) / 2.0, 0.0);
}

Vector ElectronModel::totalBfield(const Vector vx, const double B0,
    const double r0, const double a[]) {
    double x = vx.x();
    double y = vx.y();
    double z = vx.z();
    return solenoidBfield(x, y, z, a) + hexapoleBfield(x, y, B0, r0);
}

