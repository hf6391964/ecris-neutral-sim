#include <fstream>
#include <algorithm>
#include <sstream>
#include <cstddef>
#include <cmath>

#include "flychkparser.h"

FlychkParser::FlychkParser(const std::string &filename,
    std::vector<std::string> datasetHeaders) : datasetHeaders_(datasetHeaders) {
    // Parse rr and dr datasets
    for (std::string header : datasetHeaders_) {
        parseFlychkData(filename, header);
    }
}

bool FlychkParser::parseFlychkData(std::string filename,
    std::string dataset) {
    std::ifstream datafile(filename);
    if (!datafile.is_open()) return false;

    // Seek the desired header in the file
    std::string header = dataset + " temp";
    std::string line;
    bool blockFound = false;
    while (!blockFound && std::getline(datafile, line)) {
        if (line.find(header) != std::string::npos) {
            blockFound = true;
        }
    }

    if (!blockFound) return false;

    std::stringstream linestream(line);
    linestream.ignore(header.length());

    // Count the number of fields
    int nFields = 1;
    while (linestream) {
        char c = linestream.get();
        if (c != ' ' && c != EOF) {
            nFields += 1;
            linestream.ignore(std::numeric_limits<int>::max(), ' ');
        }
    }

    // Parse lines as long as they have the desired number of fields
    int nParsedLines = 0;
    std::vector<std::vector<double>> rows;
    std::vector<double> temperatures;
    while (std::getline(datafile, line)) {
        std::vector<double> fields;
        fields.reserve(nFields - 1);
        std::stringstream linestream1(line);
        int nParsedFields = 0;
        double temperature = 0;
        while (linestream1) {
            float f;
            if (linestream1 >> f) {
                nParsedFields += 1;
                if (nParsedFields == 1) {
                    temperature = f;
                } else {
                    // The data are in cm3/s, normalize to m3/s
                    fields.push_back(f * 1e-6);
                }
            } else {
                break;
            }
        }

        if (nParsedFields == nFields) {
            nParsedLines += 1;
            rows.push_back(fields);
            temperatures.push_back(temperature);
        } else {
            break;
        }
    }

    FlychkDataset dset;
    dset.rateCoefficients = rows;
    dset.temperatures = temperatures;
    dset.maxChargeState = nFields;

    datasets_[dataset] = dset;

    datafile.close();

    return true;
}

bool FlychkParser::getRateCoefficient(double T, int chargeState,
    std::string header, double &result, bool logInterp) const {
    DatasetMap::const_iterator it = datasets_.find(header);
    if (it == datasets_.end()) {
        return false;
    }

    FlychkDataset dset = (*it).second;

    if (chargeState < 1 || chargeState > dset.maxChargeState) {
        return false;
    }

    std::vector<double>::const_iterator temp_it =
        std::upper_bound(dset.temperatures.begin(), dset.temperatures.end(), T);

    if (temp_it == dset.temperatures.begin() ||
        temp_it == dset.temperatures.end()) {
        return false;
    }

    ptrdiff_t idx = temp_it - dset.temperatures.begin();

    double tempLow = dset.temperatures[idx - 1],
        tempHigh = dset.temperatures[idx],
        rcoeffLow = dset.rateCoefficients[idx - 1][chargeState - 1],
        rcoeffHigh = dset.rateCoefficients[idx][chargeState - 1];

    if (logInterp) {
        tempLow = std::log(tempLow);
        tempHigh = std::log(tempHigh);
        T = std::log(T);
    }

    double k = (rcoeffHigh - rcoeffLow) / (tempHigh - tempLow);
    result = (T - tempLow) * k + rcoeffLow;

    return true;
}

bool FlychkParser::getTotalRateCoefficient(double T, int chargeState,
    double &result, bool logInterp) const {
    double res = 0.0;
    for (auto header : datasetHeaders_) {
        double tmp;
        if (getRateCoefficient(T, chargeState, header, tmp, logInterp)) {
            res += tmp;
        } else {
            return false;
        }
    }

    result = res;
    return true;
}

