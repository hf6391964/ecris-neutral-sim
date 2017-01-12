#include <fstream>

#include "flychkparser.h"

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

    int nFields = 1;
    while (linestream) {
        char c = linestream.get();
        if (c != ' ' && c != EOF) {
            nFields += 1;
            linestream.ignore(std::numeric_limits<int>::max(), ' ');
        }
    }

    int nParsedLines = 0;
    while (std::getline(datafile, line)) {
        linestream = std::stringstream(line);
        int nParsedFields = 0;
        while (linestream) {
            float f;
            if (linestream >> f) {
                nParsedFields += 1;
            } else {
                break;
            }
        }

        if (nParsedFields == nFields) {
            nParsedLines += 1;
        } else {
            break;
        }
    }

    datafile.close();

    return true;
}

