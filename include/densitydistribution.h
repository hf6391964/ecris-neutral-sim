#pragma once

#include "cgal_and_typedefs.h"
#include "grid.h"

class DensityDistribution {
    double *densityVector_ = NULL;
    Grid grid_;

    public:
        DensityDistribution(const DensityDistribution &src, double weight = 1.0);
        DensityDistribution(std::string filename, double weight = 1.0);
        ~DensityDistribution();

        double getDensityAt(const Point &p) const;

        void setCoordinateTransformation(const Aff_transformation
            &transformation);
        void removeCoordinateTransformation();
};

