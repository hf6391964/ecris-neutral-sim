#include "plasmadensities.h"

PlasmaDensities::PlasmaDensities(std::string electronDensityFilename,
    double electronWeight, std::vector<double> ionRelativeDensities) {
    std::ifstream fin(electronDensityFilename);
    // Read grid dimensions from given file
    grid_ = Grid(fin);

    // Read electron data from that file
    size_t n = grid_.arraySize();
    ionDensities_.reserve(ionRelativeDensities.size());
    ionDensities_.push_back(new double[n]);
    // Populate electron desities first.
    for (size_t i = 0; i < n; i++) {
        unsigned long nElectrons;
        fin >> nElectrons;
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        ionDensities_[0][i] = nElectrons * electronWeight;
    }
    fin.close();

    // Populate the ion densities
    for (double relDensity : ionRelativeDensities) {
        double *densityArr = new double[n];
        for (size_t i = 0; i < n; i++) {
            densityArr[i] = ionDensities_[0][i] * relDensity;
        }
        ionDensities_.push_back(densityArr);
    }
}

PlasmaDensities::~PlasmaDensities() {
    for (double *arr : ionDensities_) {
        if (arr != NULL) {
            delete[] arr;
        }
    }
}
