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
                bbox.xmin() + intervalsX_ * gridSize,
                bbox.ymin() + intervalsY_ * gridSize,
                bbox.zmin() + intervalsZ_ * gridSize
            ).bbox();
        }

        std::tuple<size_t, size_t, size_t> dimensions() {
            return std::make_tuple(intervalsX_, intervalsY_, intervalsZ_);
        }

        size_t arraySize() {
            return intervalsX_ * intervalsY_ * intervalsZ_;
        }

        bool arrayIndex(double x, double y, double z, size_t& i) {
            size_t ix = std::floor((x - bbox_.xmin()) / xlen_);
            size_t iy = std::floor((y - bbox_.ymin()) / ylen_);
            size_t iz = std::floor((z - bbox_.zmin()) / zlen_);

            if (ix < intervalsX_ && iy < intervalsY_ && iz < intervalsZ_) {
                i = ix + intervalsX_ * (iy + intervalsY_ * iz);
                return true;
            }

            return false;
        }

        bool arrayIndex(Point p, size_t& i) {
            return arrayIndex(p.x(), p.y(), p.z(), i);
        }

        void writeDimensions(std::ostream& os) {
            os << intervalsX_ << "," << intervalsY_ << "," << intervalsZ_ <<
                std::endl;
            os << bbox_.xmin() << "," << bbox_.ymin() << "," << bbox_.zmin() <<
                std::endl;
            os << bbox_.xmax() << "," << bbox_.ymax() << "," << bbox_.zmax() <<
                std::endl;
        }
};

#endif
