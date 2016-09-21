#ifndef GRID_H
#define GRID_H

#include <tuple>

#include "cgal_and_typedefs.h"

class Grid {
    Bbox bbox_;
    size_t intervalsX_, intervalsY_, intervalsZ_;
    double gridSize_;

    public:
        Grid(Bbox bbox, double gridSize) : gridSize_(gridSize) {
            intervalsX_ = std::ceil((bbox.xmax() - bbox.xmin()) / gridSize);
            intervalsY_ = std::ceil((bbox.ymax() - bbox.ymin()) / gridSize);
            intervalsZ_ = std::ceil((bbox.zmax() - bbox.zmin()) / gridSize);

            bbox_ = bbox + Point(
                bbox.xmin() + intervalsX_ * gridSize,
                bbox.ymin() + intervalsY_ * gridSize,
                bbox.zmin() + intervalsZ_ * gridSize
            ).bbox();
        }

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

        bool arrayIndex(double x, double y, double z, size_t& i) const {
            size_t ix = std::floor((x - bbox_.xmin()) / gridSize_);
            size_t iy = std::floor((y - bbox_.ymin()) / gridSize_);
            size_t iz = std::floor((z - bbox_.zmin()) / gridSize_);

            if (ix < intervalsX_ && iy < intervalsY_ && iz < intervalsZ_) {
                i = ix + intervalsX_ * (iy + intervalsY_ * iz);
                return true;
            }

            return false;
        }

        bool arrayIndex(Point p, size_t& i) const {
            return arrayIndex(p.x(), p.y(), p.z(), i);
        }

        void writeDimensions(std::ostream& os) const {
            os << "# No. of X intervals" << CSV_SEP << " Y intervals" <<
                CSV_SEP << " Z intervals" << std::endl;
            os << intervalsX_ << CSV_SEP << intervalsY_ << CSV_SEP <<
                intervalsZ_ << std::endl;
            os << "# Xmin" << CSV_SEP << " Ymin" << CSV_SEP << " Zmin" <<
                std::endl;
            os << bbox_.xmin() << CSV_SEP << bbox_.ymin() << CSV_SEP <<
                bbox_.zmin() << std::endl;
            os << "# Xmax" << CSV_SEP << " Ymax" << CSV_SEP << " Zmax" << std::endl;
            os << bbox_.xmax() << CSV_SEP << bbox_.ymax() << CSV_SEP <<
                bbox_.zmax() << std::endl;
        }
};

#endif
