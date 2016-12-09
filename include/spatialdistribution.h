#pragma once

#include "cgal_and_typedefs.h"
#include "grid.h"

template<typename T>
class SpatialDistribution {
    protected:
        T *valueVector_ = NULL;
        Grid grid_;

    private:
        virtual T getNull() const = 0;

    public:
        SpatialDistribution() {}
        SpatialDistribution(const Grid &grid);
        template<typename T1>
        SpatialDistribution(const SpatialDistribution<T1> &src,
            double weight = 1.0);
        ~SpatialDistribution();

        T getValueAt(const Point &p) const;

        void setCoordinateTransformation(const Aff_transformation
            &transformation);
        void removeCoordinateTransformation();
};

class DensityDistribution : public SpatialDistribution<double> {
    private:
        double getNull() const { return 0.0; }

    public:
        DensityDistribution(std::string filename, double weight = 1.0);
        template<typename T>
        DensityDistribution(const SpatialDistribution<T> &src,
            double weight = 1.0);
};

class VelocityDistribution : public SpatialDistribution<Vector> {
    private:
        Vector getNull() const { return Vector(0.0, 0.0, 0.0); }
};

