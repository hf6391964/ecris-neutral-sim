#pragma once

#include <tuple>
#include <fstream>

#include "cgal_and_typedefs.h"
#include "util.h"

class Grid {
    Bbox bbox_;
    unsigned int intervalsX_, intervalsY_, intervalsZ_;
    double gridSize_, gridSizeInverse_;
    double xmin_, ymin_, zmin_;
    Aff_transformation coordTransformation_;
    bool doTransform_ = false;

    public:
        Grid() {};

        Grid(Bbox bbox, double gridSize);

        Grid(std::ifstream &fin);

        double getXatIndex(size_t i) const {
            return bbox_.xmin() + i*gridSize_;
        }

        double getYatIndex(size_t i) const {
            return bbox_.ymin() + i*gridSize_;
        }

        double getZatIndex(size_t i) const {
            return bbox_.zmin() + i*gridSize_;
        }

        std::tuple<unsigned int, unsigned int, unsigned int> dimensions() const {
            return std::make_tuple(intervalsX_, intervalsY_, intervalsZ_);
        }

        size_t arraySize() const {
            return intervalsX_ * intervalsY_ * intervalsZ_;
        }

        bool arrayIndex(const Point &p, size_t& i) const;

        bool arrayIndex(const double &x, const double &y, const double &z,
            size_t& i) const {
            return arrayIndex(Point(x, y, z), i);
        }

        void writeDimensions(std::ostream& os) const;

        // Set the coordinate transformation which maps coordinates from the
        // target geometry frame to the frame where the electron densities are
        // calculated. Typical uses would be translations and some 90 degree
        // rotations.
        // The default is identity where no transformation is performed.
        void setCoordinateTransformation(const Aff_transformation
            &transformation);
        // Removes the current coordinate transformation.
        void removeCoordinateTransformation();
};

