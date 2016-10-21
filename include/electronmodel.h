#ifndef ELECTRONMODEL_H
#define ELECTRONMODEL_H

#include <cmath>

#include "cgal_and_typedefs.h"
#include "constants.h"

class ElectronModel {
    public:
        ElectronModel(const double B0, const double r0, const double dt,
            const double z1, const double z2, double a[] = SOLENOID_FIELD_AI);

        void newParticle(const double energy, Rng& rng);
        Vector velocity() const;
        Vector position() const;
        Vector energy() const;
        bool moveParticle();

        static Vector totalBfield(const Vector vx, const double B0,
            const double r0, const double a[] = SOLENOID_FIELD_AI);
        static Vector solenoidBfield(const double x, const double y,
            const double z, const double a[]);
        static Vector hexapoleBfield(const double x, const double y,
            const double B0, const double r0);

    private:
        Vector particleVelocity_;
        Vector particlePosition_;
        double B0_, r0_, dt_, t_, z1_, z2_;
        double* a_;
};

#endif

