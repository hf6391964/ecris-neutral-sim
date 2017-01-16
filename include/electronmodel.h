#pragma once

#include <cmath>
#include <functional>
#include <numeric>
#include <string>
#include <fstream>

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

        template<typename Criterion>
        void newParticle(const double speed, Criterion criterion, Rng& rng) {
            double u = 2.0 * uni01(rng) - 1.0;
            double theta = 2.0*M_PI*uni01(rng);
            double sqrtu = std::sqrt(1.0 - u*u);

            Vector B;
            do {
                // Pick a point uniformly on a cylinder
                double z = z1_ + (z2_ - z1_) * uni01(rng);
                double r = r0_ * sqrt(uni01(rng));
                theta = 2.0*M_PI * uni01(rng);
                particlePosition_ = Vector(r*std::cos(theta), r*std::sin(theta), z);
                B = totalBfield(particlePosition_, B0_, r0_, a_);
            } while (!criterion(B));

            Vector v = speed * Vector(
                sqrtu * std::cos(theta),
                sqrtu * std::sin(theta),
                u
            );

            int nSteps = 10000;
            double dt2 = 0.5/nSteps * dt_;
            // Make a half step backwards to get the "half timestep" velocity for
            // leapfrog to get a start
            Vector pos = particlePosition_;
            for (int i = 0; i < nSteps; i++) {
                v = v - dt2*(-ELECTRON_CHARGE_MASS_RATIO) *
                    CGAL::cross_product(v, B);
                pos = pos - dt2*v;
            }

            particleVelocity_ = v;

            t_ = 0.0;
        };

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

        void runSimulation(const unsigned long nParticles,
            std::function<double(Rng&)> getSpeed, const double Becr, Rng& rng);

        void runMonoenergeticSimulation(const unsigned long nParticles,
            const double energy, const double Becr, Rng& rng);

        void runMaxwellSimulation(const unsigned long nParticles,
            const double T_eV, const double Becr, Rng& rng);

        void writeDensityToFile(const std::string &filename,
            bool normalize = true) const;

        static Vector totalBfield(const Vector vx, const double B0,
            const double r0, const double a[] = SOLENOID_FIELD_AI);
        static Vector solenoidBfield(const double x, const double y,
            const double z, const double a[] = SOLENOID_FIELD_AI);
        static Vector hexapoleBfield(const double x, const double y,
            const double B0, const double r0);

        void writeElectronEndpoints(const std::string &z1Filename,
            const std::string &z2Filename,
            const std::string &radialFilename) const;

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

