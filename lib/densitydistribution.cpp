#include "densitydistribution.h"

DensityDistribution::DensityDistribution(const DensityDistribution &src,
    double weight) {
    grid_ = src.grid_;
    size_t n = grid_.arraySize();
    densityVector_ = new double[n];
    for (size_t i = 0; i < n; ++i) {
        densityVector_[i] = src.densityVector_[i] * weight;
    }
}

DensityDistribution::DensityDistribution(std::string filename,
    double weight) {
    std::ifstream fin(filename);
    // Read grid dimensions from given file
    grid_ = Grid(fin);

    // Read electron data from that file
    size_t n = grid_.arraySize();
    densityVector_ = new double[n];
    // Populate electron desities first.
    for (size_t i = 0; i < n; ++i) {
        unsigned long nParticles;
        fin >> nParticles;
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        densityVector_[i] = nParticles * weight;
    }
    fin.close();
}

DensityDistribution::~DensityDistribution() {
    if (densityVector_ != NULL) {
        delete[] densityVector_;
    }
}

double DensityDistribution::getDensityAt(const Point &p) const {
    size_t index;
    bool hasValue = grid_.arrayIndex(p, index);

    if (hasValue) {
        return densityVector_[index];
    }

    return 0.0;
}

void DensityDistribution::setCoordinateTransformation(const Aff_transformation
    &transformation) {
    grid_.setCoordinateTransformation(transformation);
}

void DensityDistribution::removeCoordinateTransformation() {
    grid_.removeCoordinateTransformation();
}

