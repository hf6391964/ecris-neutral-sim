#ifndef PLASMADENSITIES_H
#define PLASMADENSITIES_H

#include <fstream>
#include <string>
#include <vector>
#include <limits>

#include "cgal_and_typedefs.h"
#include "grid.h"

class PlasmaDensities {
    std::vector<double *> ionDensities_;
    Grid grid_;
    unsigned int maxChargeState_ = 0;

    public:
        // Electron weight is the statistical weight of an electron in the
        // electron density distribution i.e.
        // # of real particles / test particle.
        // Ion densities are relative to real electron density,
        // ionRelativeDensities[i] corresponds to charge state i+1
        // Quasineutrality should be satisfied:
        // sum(ionRelativeDensities[i] * (i+1)) = 1
        PlasmaDensities(std::string electronDensityFilename,
            double electronWeight, std::vector<double> ionRelativeDensities);

        ~PlasmaDensities();

        void setCoordinateTransformation(const Aff_transformation &tf);

        // Here ion densities are indexed such that i = q.
        // A special case is i = 0 which corresponds to electron density.
        double getIonDensityAt(const Point &p, unsigned int chargeState) const;
        double getElectronDensityAt(const Point &p) const;
};

#endif

