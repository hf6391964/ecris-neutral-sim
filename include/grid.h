#ifndef GRID_H
#define GRID_H

#include <tuple>
#include <fstream>

#include "cgal_and_typedefs.h"

class Grid {
    Bbox bbox_;
    size_t intervalsX_, intervalsY_, intervalsZ_;
    double gridSize_;

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

        std::tuple<size_t, size_t, size_t> dimensions() const {
            return std::make_tuple(intervalsX_, intervalsY_, intervalsZ_);
        }

        size_t arraySize() const {
            return intervalsX_ * intervalsY_ * intervalsZ_;
        }

        bool arrayIndex(Point p, size_t& i) const {
            return arrayIndex(p.x(), p.y(), p.z(), i);
        }

        bool arrayIndex(double x, double y, double z, size_t& i) const;

        void writeDimensions(std::ostream& os) const;
};

#endif
