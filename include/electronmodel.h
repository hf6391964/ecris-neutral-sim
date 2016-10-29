#ifndef ELECTRONMODEL_H
#define ELECTRONMODEL_H

#include <cmath>
#include <functional>
#include <numeric>

#include "cgal_and_typedefs.h"
#include "constants.h"
#include "grid.h"

class ElectronModel {
    public:
        ElectronModel(const double B0, const double r0, const double dt,
            const double z1, const double z2, const double gridSize,
            const double confinementTime,
            const double a[] = SOLENOID_FIELD_AI);

        ~ElectronModel();

        void resetCounters();

        void newParticle(const double energy,
            std::function<bool(Vector)> criterion, Rng& rng);

        Vector velocity() const {
            return particleVelocity_;
        }

        Vector position() const {
            return particlePosition_;
        }

        std::vector<Vector> z1CollisionPoints() const {
            return z1CollisionPoints_;
        }

        std::vector<Vector> z2CollisionPoints() const {
            return z2CollisionPoints_;
        }

        std::vector<Vector> cylinderCollisionPoints() const {
            return cylinderCollisionPoints_;
        }

        std::vector<std::tuple<double, double>> lostVelocities() const {
            return lostVelocities_;
        }

        std::vector<std::tuple<double, double>> nonLostVelocities() const {
            return nonLostVelocities_;
        }

        double energy() const {
            return 0.5 * ELECTRON_MASS_EV * velocity().squared_length() /
                (SPEED_OF_LIGHT*SPEED_OF_LIGHT);
        }

        double meanEnergy() const;

        double energyStdDev() const;

        bool moveParticle();

        void runMonoenergeticSimulation(const unsigned long nParticles,
            const double energy, const double BnormLimit, Rng& rng);

        static Vector totalBfield(const Vector vx, const double B0,
            const double r0, const double a[] = SOLENOID_FIELD_AI);
        static Vector solenoidBfield(const double x, const double y,
            const double z, const double a[] = SOLENOID_FIELD_AI);
        static Vector hexapoleBfield(const double x, const double y,
            const double B0, const double r0);

    private:
        Vector particleVelocity_;
        Vector particlePosition_;
        const double B0_, r0_, dt_, z1_, z2_, gridSize_, confinementTime_;
        const double* a_;
        double t_;
        Grid grid_;
        std::vector<Vector> z1CollisionPoints_;
        std::vector<Vector> z2CollisionPoints_;
        std::vector<Vector> cylinderCollisionPoints_;
        std::vector<double> finalEnergies_;
        std::vector<std::tuple<double, double>> lostVelocities_;
        std::vector<std::tuple<double, double>> nonLostVelocities_;
        unsigned long* particleCount_ = NULL;
};

#endif

