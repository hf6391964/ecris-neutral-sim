#pragma once

#include <string>

#include "neutralizationchannel.h"
#include "surfacecollection.h"

class WallNeutralization : NeutralizationChannel {
    std::vector<IntersectionPoint> wallPoints_, end1Points_, end2Points_;
    std::vector<std::vector<IntersectionPoint> *> pointSets_;
    double confinementTime_;

    public:
        WallNeutralization(double confinementTime,
            const std::string &wallFilename,
            const std::string &end1Filename,
            const std::string &end2Filename,
            const SurfaceCollection &surfaces);
};

