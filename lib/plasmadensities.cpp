#include "plasmadensities.h"

PlasmaDensities::PlasmaDensities(std::string electronDensityFilename,
    double electronWeight, std::vector<double> ionWeights) {
    std::ifstream fin(electronDensityFilename);

    grid_ = Grid(fin);

    fin.close();
}

PlasmaDensities::~PlasmaDensities() {
}
