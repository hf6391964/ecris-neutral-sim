#pragma once

#include <vector>
#include <string>
#include <map>

struct FlychkDataset {
    std::vector<double> temperatures;
    std::vector<std::vector<double>> rateCoefficients;
    int maxChargeState;
};

// Parse radiative and dielectronic recombination rate coefficients from the
// supplied FLYCHK data file.
//
// Compute total recombination rate coefficients by logarithmically
// interpolating and summing these.

typedef std::map<std::string, FlychkDataset> DatasetMap;

class FlychkParser {
    std::vector<std::string> datasetHeaders_;
    DatasetMap datasets_;

    bool parseFlychkData(std::string filename, std::string dataset);

    public:
        FlychkParser(const std::string &filename,
            std::vector<std::string> datasetHeaders = { "rr", "dr" });

        bool getRateCoefficient(double T, int chargeState,
            std::string header, double &result, bool logInterp = true) const;
        bool getTotalRateCoefficient(double T, int chargeState,
            double &result, bool logInterp = true) const;
};

