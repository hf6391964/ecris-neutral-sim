#ifndef PLASMADENSITIES_H
#define PLASMADENSITIES_H

#include <fstream>
#include <string>
#include <vector>
#include <limits>

#include "cgal_and_typedefs.h"
#include "grid.h"

class PlasmaDensities {
    std::vector<double *> ionDensities;
    unsigned int maxChargeState = 0;
    Grid grid_;

    public:
        PlasmaDensities(std::string electronDensityFilename,
            double electronWeight, std::vector<double> ionDensities);

        ~PlasmaDensities();

        void setCoordinateTransformation();
        void removeCoordinateTransformation();

        double getIonDensityAt(Point p, unsigned int chargeState) const;
        double getElectronDensityAt(Point p) const;
};

#endif

