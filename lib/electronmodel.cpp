#include "electronmodel.h"

ElectronModel::ElectronModel(const double B0, const double r0, const double dt,
    const double z1, const double z2, double a[]) : B0_(B0), r0_(r0), dt_(dt),
    z1_(z1), z2_(z2), a_(a) {
}

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

void ElectronModel::newParticle(const double energy, Rng& rng) {
    double speed = std::sqrt(2.0*energy / ELECTRON_MASS_EV) * SPEED_OF_LIGHT;

    double u = 2.0 * uni01(rng) - 1.0;
    double theta = 2.0*M_PI*uni01(rng);
    double sqrtu = std::sqrt(1.0 - u*u);

    // Pick a point uniformly on a cylinder
    double z = z1_ + (z2_ - z1_) * uni01(rng);
    double r = r0_ * sqrt(uni01(rng));
    theta = 2.0*M_PI * uni01(rng);
    particlePosition_ = Vector(r*std::cos(theta), r*std::sin(theta), z);

    Vector v = speed * Vector(
        sqrtu * std::cos(theta),
        sqrtu * std::sin(theta),
        u
    );
    // Make a half step backwards to get the "half timestep" velocity for
    // leapfrog to get a start
    particleVelocity_ = v - 0.5*dt_*(-ELECTRON_CHARGE_MASS_RATIO) *
        CGAL::cross_product(v, totalBfield(particlePosition_, B0_, r0_, a_));
}

bool ElectronModel::moveParticle() {
    Vector B = totalBfield(particlePosition_, B0_, r0_, a_);

    // Boris integrator with zero E field

    double Bsqnorm = B.squared_length();
    double Bnorm = std::sqrt(Bsqnorm);
    double f1 = std::tan(-ELECTRON_CHARGE_MASS_RATIO * dt_ * Bnorm / 2.0) /
        Bnorm;
    double f2 = 2.0 * f1 / (1.0 + f1*f1 * Bsqnorm);

    Vector v = particleVelocity_ + f1*CGAL::cross_product(particleVelocity_, B);
    particleVelocity_ = particleVelocity_ + f2*CGAL::cross_product(v, B);
    particlePosition_ = particlePosition_ + dt_ * particleVelocity_;

    return false;
}

