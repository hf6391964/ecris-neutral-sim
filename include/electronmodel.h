#ifndef ELECTRONMODEL_H
#define ELECTRONMODEL_H

#include <cmath>
#include <functional>

#include "cgal_and_typedefs.h"
#include "constants.h"

class ElectronModel {
    public:
        ElectronModel(const double B0, const double r0, const double dt,
            const double z1, const double z2,
            const double a[] = SOLENOID_FIELD_AI);

        void newParticle(const double energy,
            std::function<bool(Vector)> criterion, Rng& rng);

        Vector velocity() const {
            return particleVelocity_;
        }

        Vector position() const {
            return particlePosition_;
        }

        double energy() const {
            return 0.5 * ELECTRON_MASS_EV * velocity().squared_length() /
                SPEED_OF_LIGHT*SPEED_OF_LIGHT;
        }

        bool moveParticle();

        static Vector totalBfield(const Vector vx, const double B0,
            const double r0, const double a[] = SOLENOID_FIELD_AI);
        static Vector solenoidBfield(const double x, const double y,
            const double z, const double a[] = SOLENOID_FIELD_AI);
        static Vector hexapoleBfield(const double x, const double y,
            const double B0, const double r0);

    private:
        Vector particleVelocity_;
        Vector particlePosition_;
        double B0_, r0_, dt_, t_, z1_, z2_;
        const double* a_;
};

#endif

