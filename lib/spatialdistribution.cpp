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
}

