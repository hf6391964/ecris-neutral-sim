#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <limits>

#include "cgal_and_typedefs.h"
#include "grid.h"

class PlasmaDensities {
    std::vector<double *> ionDensities_;
    std::vector<double> ionTemperatures_;
    Grid grid_;
    unsigned int maxChargeState_ = 0;
    double ionMassEv_ = 0.0;

    public:
        // Electron weight is the statistical weight of an electron in the
        // electron density distribution i.e.
        // # of real particles / test particle.
        //
        // Ion densities are relative to real electron density,
        // ionRelativeDensities[i] corresponds to charge state i+1
        // Quasineutrality should be satisfied:
        // sum(ionRelativeDensities[i] * (i+1)) = 1
        //
        // electron temperature is eV
        //
        // ion temperatures (eV) at indices corresponding to charge state
        //
        // ion mass should be in eV
        PlasmaDensities(std::string electronDensityFilename,
            double electronWeight, std::vector<double> ionRelativeDensities,
            double electronTemperature,
            std::vector<double> ionTemperatures, double ionMassEv);

        ~PlasmaDensities();

        void setCoordinateTransformation(const Aff_transformation &tf);

        // Here ion densities are indexed such that i = q.
        // A special case is i = 0 which corresponds to electron density.
        double getIonDensityAt(const Point &p, unsigned int chargeState) const;
        double getElectronDensityAt(const Point &p) const;

        double getIonTemperature(unsigned int chargeState) const {
            return ionTemperatures_[chargeState];
        }
        double getElectronTemperature() const {
            return getIonTemperature(0);
        }

        double getIonMass() const {
            return ionMassEv_;
        }

        Vector getIsotropicIonVelocity(Rng &rng, unsigned int chargeState)
            const;
        Vector getIsotropicElectronVelocity(Rng &rng) const;
};

