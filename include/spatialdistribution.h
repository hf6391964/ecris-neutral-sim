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
        ~SpatialDistribution();

        void initializeTo(T value);
        void initializeToNull();

        virtual T getValueAt(const Point &p) const;

        void setCoordinateTransformation(const Aff_transformation
            &transformation);
        void removeCoordinateTransformation();
};

class DensityDistribution : public SpatialDistribution<double> {
    private:
        double getNull() const { return 0.0; }

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

class VelocityDistribution : public SpatialDistribution<Vector> {
    private:
        Vector getNull() const { return Vector(0.0, 0.0, 0.0); }
};

template class SpatialDistribution<double>;
template class SpatialDistribution<Vector>;
