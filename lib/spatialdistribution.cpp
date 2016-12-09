#include "spatialdistribution.h"

template<typename T>
SpatialDistribution<T>::SpatialDistribution(const Grid &grid) {
    size_t n = grid_.arraySize();
    valueVector_ = new T[n];
    grid_ = grid;
    for (size_t i = 0; i < n; ++i) {
        valueVector_[i] = getNull();
    }
}

template<typename T>
SpatialDistribution<T>::SpatialDistribution(const SpatialDistribution<T> &src,
    double weight) {
    grid_ = src.grid_;
    size_t n = grid_.arraySize();
    valueVector_ = new T[n];
    for (size_t i = 0; i < n; ++i) {
        valueVector_[i] = src.valueVector_[i] * weight;
    }
}

template<typename T>
SpatialDistribution<T>::~SpatialDistribution() {
    if (valueVector_ != NULL) {
        delete[] valueVector_;
    }
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

DensityDistribution::DensityDistribution(std::string filename,
    double weight) {
    std::ifstream fin(filename);
    // Read grid dimensions from given file
    grid_ = Grid(fin);

    // Read electron data from that file
    size_t n = grid_.arraySize();
    valueVector_ = new double[n];
    // Populate electron desities first.
    for (size_t i = 0; i < n; ++i) {
        double nParticles;
        fin >> nParticles;
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        valueVector_[i] = nParticles * weight;
    }
    fin.close();
}

