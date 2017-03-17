#pragma once

#include <tuple>
#include <fstream>

#include "cgal_and_typedefs.h"
#include "util.h"

typedef std::tuple<unsigned int, unsigned int, unsigned int> Index3D;

class Grid {
    Bbox bbox_;
    unsigned int intervalsX_, intervalsY_, intervalsZ_;
    double gridSize_, gridSizeInverse_;
    double xmin_, ymin_, zmin_;
    Aff_transformation coordTransformation_;
    bool doTransform_ = false;

    public:
        Grid() : coordTransformation_(CGAL::IDENTITY) {};
        Grid(Bbox bbox, double gridSize);
        Grid(std::ifstream &fin);

        std::tuple<unsigned int, unsigned int, unsigned int> dimensions() const;
        double cellVolume() const;
        double cellSideLength() const;
        size_t arraySize() const;
        bool arrayIndex3D(const Point &p, Index3D &i) const;
        bool arrayIndex3D(const double &x, const double &y, const double &z,
            Index3D &i) const;
        bool arrayIndex(const Point &p, size_t& i) const;
        bool arrayIndex(const double &x, const double &y, const double &z,
            size_t& i) const;
        bool getCellMidpoint(size_t i, Point &p) const;
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

