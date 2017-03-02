#include "grid.h"
#include "constants.h"

Grid::Grid(Bbox bbox, double gridSize) : Grid() {
    gridSize_ = gridSize;
    intervalsX_ = std::ceil((bbox.xmax() - bbox.xmin()) / gridSize);
    intervalsY_ = std::ceil((bbox.ymax() - bbox.ymin()) / gridSize);
    intervalsZ_ = std::ceil((bbox.zmax() - bbox.zmin()) / gridSize);

    bbox_ = bbox + Point(
        bbox.xmin() + intervalsX_ * gridSize,
        bbox.ymin() + intervalsY_ * gridSize,
        bbox.zmin() + intervalsZ_ * gridSize
    ).bbox();

    xmin_ = bbox.xmin();
    ymin_ = bbox.ymin();
    zmin_ = bbox.zmin();

    removeCoordinateTransformation();
    gridSizeInverse_ = 1.0 / gridSize_;
}

Grid::Grid(std::ifstream &fin) : Grid() {
    double xmin, xmax, ymin, ymax, zmin, zmax;

    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    fin >> intervalsX_;
    fin.get();
    fin >> intervalsY_;
    fin.get();
    fin >> intervalsZ_;
    fin.get();

    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    fin >> xmin;
    fin.get();
    fin >> ymin;
    fin.get();
    fin >> zmin;
    fin.get();

    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    fin >> xmax;
    fin.get();
    fin >> ymax;
    fin.get();
    fin >> zmax;
    fin.get();

    bbox_ = Bbox(xmin, ymin, zmin, xmax, ymax, zmax);
    xmin_ = bbox_.xmin();
    ymin_ = bbox_.ymin();
    zmin_ = bbox_.zmin();
    gridSize_ = (xmax - xmin) / intervalsX_;
    gridSizeInverse_ = 1.0 / gridSize_;
}

std::tuple<unsigned int, unsigned int, unsigned int> Grid::dimensions() const {
    return std::make_tuple(intervalsX_, intervalsY_, intervalsZ_);
}

double Grid::cellVolume() const {
    return gridSize_ * gridSize_ * gridSize_;
}

double Grid::cellSideLength() const {
    return gridSize_;
}

size_t Grid::arraySize() const {
    return intervalsX_ * intervalsY_ * intervalsZ_;
}

bool Grid::arrayIndex(const double &x, const double &y, const double &z,
    size_t& i) const {
    return arrayIndex(Point(x, y, z), i);
}

void Grid::writeDimensions(std::ostream& os) const {
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

void Grid::setCoordinateTransformation(const Aff_transformation
    &transformation) {
    coordTransformation_ = transformation;
    doTransform_ = true;
}

void Grid::removeCoordinateTransformation() {
    coordTransformation_ = Aff_transformation(CGAL::IDENTITY);
    doTransform_ = false;
}

bool Grid::arrayIndex(const Point &p, size_t& i) const {
    const Point &pTransformed = doTransform_ ?
        coordTransformation_.transform(p) : p;
    unsigned int ix = Util::fastFloor((pTransformed.x() - xmin_) * gridSizeInverse_);
    unsigned int iy = Util::fastFloor((pTransformed.y() - ymin_) * gridSizeInverse_);
    unsigned int iz = Util::fastFloor((pTransformed.z() - zmin_) * gridSizeInverse_);

    if (ix < intervalsX_ && iy < intervalsY_ && iz < intervalsZ_) {
        i = ix + intervalsX_ * (iy + intervalsY_ * iz);
        return true;
    }

    return false;
}

bool Grid::getCellMidpoint(size_t i, Point &p) const {
    if (i >= arraySize()) return false;

    unsigned int iz = i / (intervalsX_*intervalsY_);
    unsigned int iy = (i - iz*intervalsX_*intervalsY_) / intervalsX_;
    unsigned int ix = i - intervalsX_ * (intervalsY_*iz + iy);

    p = Point(
        bbox_.xmin() + (ix + 0.5) * gridSize_,
        bbox_.ymin() + (iy + 0.5) * gridSize_,
        bbox_.zmin() + (iz + 0.5) * gridSize_
    );

    return true;
}
