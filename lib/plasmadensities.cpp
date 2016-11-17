#include "plasmadensities.h"

PlasmaDensities::PlasmaDensities(std::string electronDensityFilename,
    double electronWeight, std::vector<double> ionRelativeDensities,
    double electronTemperature, std::vector<double> ionTemperatures,
    double ionMassEv) : ionMassEv_(ionMassEv) {
    std::ifstream fin(electronDensityFilename);
    // Read grid dimensions from given file
    grid_ = Grid(fin);

    // Read electron data from that file
    size_t n = grid_.arraySize();
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

    ionTemperatures_.push_back(electronTemperature);
    ionTemperatures_.insert(ionTemperatures_.end(), ionTemperatures.begin(),
        ionTemperatures.end());

    maxChargeState_ = ionDensities_.size() - 1;
}

PlasmaDensities::~PlasmaDensities() {
    for (double *arr : ionDensities_) {
        if (arr != NULL) {
            delete[] arr;
        }
    }
}

void PlasmaDensities::setCoordinateTransformation(
    const Aff_transformation &tf) {
    grid_.setCoordinateTransformation(tf);
}

double PlasmaDensities::getIonDensityAt(const Point &p,
    unsigned int chargeState) const {
    size_t index;
    bool hasValue = grid_.arrayIndex(p, index);

    if (hasValue && chargeState <= maxChargeState_) {
        return ionDensities_[chargeState][index];
    }

    return 0.0;
}

double PlasmaDensities::getElectronDensityAt(const Point &p) const {
    return getIonDensityAt(p, 0);
}

Vector PlasmaDensities::getIsotropicIonVelocity(Rng &rng,
    unsigned int chargeState) const {
    return Util::getMBVelocity(rng, getIonTemperature(chargeState), ionMassEv_);
}

Vector PlasmaDensities::getIsotropicElectronVelocity(Rng &rng) const {
    return Util::getMBVelocity(rng, getElectronTemperature(),
        ELECTRON_MASS_EV);
}

