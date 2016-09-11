#ifndef SIMULATIONMODEL_H
#define SIMULATIONMODEL_H

#include "surface.h"
#include "particle.h"

class SimulationModel {
    std::vector<Surface> surfaces_;
    long nParticles_;
    double timestep_;
};

#endif
