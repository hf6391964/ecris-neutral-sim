#ifndef GRID_H
#define GRID_H

#include <tuple>

#include "cgal_and_typedefs.h"

class Grid {
    Bbox bbox_;
    size_t intervalsX_, intervalsY_, intervalsZ_;
    double xlen_, ylen_, zlen_;

    public:
        Grid(Bbox bbox, double gridSize) {
            xlen_ = bbox.xmax() - bbox.xmin(),
            ylen_ = bbox.ymax() - bbox.ymin(),
            zlen_ = bbox.zmax() - bbox.zmin();

            intervalsX_ = std::ceil(xlen_ / gridSize);
            intervalsY_ = std::ceil(ylen_ / gridSize);
            intervalsZ_ = std::ceil(zlen_ / gridSize);

            bbox_ = bbox + Point(
                bbox_.xmin() + intervalsX_ * gridSize,
                bbox_.ymin() + intervalsY_ * gridSize,
                bbox_.zmin() + intervalsZ_ * gridSize
            ).bbox();
        }

        std::tuple<size_t, size_t, size_t> dimensions() {
            return std::make_tuple(intervalsX_, intervalsY_, intervalsZ_);
        }

        size_t arraySize() {
            return intervalsX_ * intervalsY_ * intervalsZ_;
        }

        size_t arrayIndex(double x, double y, double z) {
            size_t ix = std::floor((x - bbox_.xmin()) / xlen_);
            size_t iy = std::floor((y - bbox_.ymin()) / ylen_);
            size_t iz = std::floor((z - bbox_.zmin()) / zlen_);

            return ix + intervalsX_ * (iy + intervalsY_ * iz);
        }

        size_t arrayIndex(Point p) {
            size_t ix = std::floor((p.x() - bbox_.xmin()) / xlen_);
            size_t iy = std::floor((p.y() - bbox_.ymin()) / ylen_);
            size_t iz = std::floor((p.z() - bbox_.zmin()) / zlen_);

            return ix + intervalsX_ * (iy + intervalsY_ * iz);
        }
};

#endif
