#pragma once

#include "cgal_and_typedefs.h"
#include "grid.h"

template<typename T>
class SpatialDistribution {
    protected:
        T *valueVector_ = NULL;
        Grid grid_;

    private:
        T null_value_;

    public:
        SpatialDistribution(T null_value) : null_value_(null_value) {}
        SpatialDistribution(const Grid &grid, T null_value);
        ~SpatialDistribution();

        void initializeTo(T value);
        void initializeToNull();

        T getValueAt(const Point &p) const;

        void setCoordinateTransformation(const Aff_transformation
            &transformation);
        void removeCoordinateTransformation();
};

class DensityDistribution : public SpatialDistribution<double> {
    private:
        double sumDensity_;
        std::shared_ptr<std::vector<double>> cumulativeDensity_;
        void calculateCumulativeDensity();
        std::shared_ptr<DensityDistribution> sourceDistribution_;
        double sourceDistributionWeight_;

    public:
        DensityDistribution(std::string filename, double weight = 1.0);
        DensityDistribution(const DensityDistribution &src,
            double weight = 1.0);
        DensityDistribution(std::shared_ptr<DensityDistribution> src,
            double weight = 1.0);

        double getValueAt(const Point &p) const;

        Point getRandomPosition(Rng &rng) const;
};

template class SpatialDistribution<double>;

