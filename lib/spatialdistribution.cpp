#include "spatialdistribution.h"

template<typename T>
SpatialDistribution<T>::SpatialDistribution(const Grid &grid) {
    size_t n = grid.arraySize();
    valueVector_ = new T[n];
    grid_ = grid;
}

template<typename T>
SpatialDistribution<T>::~SpatialDistribution() {
    if (valueVector_ != NULL) {
        delete[] valueVector_;
        valueVector_ = NULL;
    }
}

template<typename T>
void SpatialDistribution<T>::initializeTo(T value) {
    size_t n = grid_.arraySize();
    for (size_t i = 0; i < n; ++i) {
        valueVector_[i] = value;
    }
}

template<typename T>
void SpatialDistribution<T>::initializeToNull() {
    initializeTo(getNull());
}

template<typename T>
T SpatialDistribution<T>::getValueAt(const Point &p) const {
    size_t index;
    bool hasValue = grid_.arrayIndex(p, index);

    if (hasValue) {
        return valueVector_[index];
    }

    return getNull();
}

template<typename T>
void SpatialDistribution<T>::setCoordinateTransformation(
    const Aff_transformation &transformation) {
    grid_.setCoordinateTransformation(transformation);
}

template<typename T>
void SpatialDistribution<T>::removeCoordinateTransformation() {
    grid_.removeCoordinateTransformation();
}

DensityDistribution::DensityDistribution(const DensityDistribution &src,
    double weight) : SpatialDistribution(src.grid_) {
    size_t n = grid_.arraySize();
    for (size_t i = 0; i < n; ++i) {
        valueVector_[i] = src.valueVector_[i] * weight;
    }
    calculateCumulativeDensity();
}

DensityDistribution::DensityDistribution(std::shared_ptr<DensityDistribution> src,
    double weight)
    : SpatialDistribution(), sourceDistribution_(src),
      sourceDistributionWeight_(weight) {
    grid_ = src->grid_;
    calculateCumulativeDensity();
}

DensityDistribution::DensityDistribution(std::string filename,
    double weight) {
    std::ifstream fin(filename);
    // Read grid dimensions from given file
    grid_ = Grid(fin);
    // Read density data from that file
    size_t n = grid_.arraySize();
    valueVector_ = new double[n];
    for (size_t i = 0; i < n; ++i) {
        double nParticles;
        fin >> nParticles;
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        valueVector_[i] = nParticles * weight;
    }
    fin.close();
    calculateCumulativeDensity();
}

void DensityDistribution::calculateCumulativeDensity() {
    if (sourceDistribution_) {
        sumDensity_ = sourceDistribution_->sumDensity_;
        cumulativeDensity_ = sourceDistribution_->cumulativeDensity_;
    } else {
        sumDensity_ = 0.0;
        size_t n = grid_.arraySize();
        cumulativeDensity_ =
            std::shared_ptr<std::vector<double>>(new std::vector<double>(n));

        for (size_t i = 0; i < n; ++i) {
            sumDensity_ += valueVector_[i];
            (*cumulativeDensity_)[i] = sumDensity_;
        }
    }
}

Point DensityDistribution::getRandomPosition(Rng &rng) const {
    double x = uni01(rng) * sumDensity_;

    size_t i = std::lower_bound(cumulativeDensity_->begin(),
        cumulativeDensity_->end(), x) - cumulativeDensity_->begin();

    double l = grid_.cellSideLength();

    Point midp;
    if (!grid_.getCellMidpoint(i, midp)) {
        throw std::out_of_range("grid index out of range");
    }

    return Point(
        midp.x() + l * (uni01(rng) - 0.5),
        midp.y() + l * (uni01(rng) - 0.5),
        midp.z() + l * (uni01(rng) - 0.5)
    );
}

double DensityDistribution::getValueAt(const Point &p) const {
    if (sourceDistribution_) {
        return sourceDistributionWeight_ * sourceDistribution_->getValueAt(p);
    }

    return SpatialDistribution::getValueAt(p);
}

