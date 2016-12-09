#pragma once

#include "cgal_and_typedefs.h"
#include "grid.h"

template<typename T>
class SpatialDistribution {
    protected:
        T *valueVector_ = NULL;
        Grid grid_;

    public:
        SpatialDistribution() {}
        SpatialDistribution(const SpatialDistribution &src,
            double weight = 1.0);
        ~SpatialDistribution();

        T getValueAt(const Point &p) const;

        void setCoordinateTransformation(const Aff_transformation
            &transformation);
        void removeCoordinateTransformation();
};

class DensityDistribution : public SpatialDistribution<double> {
    public:
        DensityDistribution(std::string filename, double weight = 1.0);
};

typedef SpatialDistribution<Vector> VelocityDistribution;

