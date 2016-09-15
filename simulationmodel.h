#ifndef SIMULATIONMODEL_H
#define SIMULATIONMODEL_H

#include <thread>

#include "surface.h"
#include "particle.h"
#include "constants.h"
#include "util.h"

class SimulationModel {
    std::vector<Surface*> surfaces_;
    long nParticles_;
    double timestep_;
    Bbox bbox_;

    void simulationThread(unsigned long nParticles, unsigned long maxSteps,
        double dt, Bbox bbox);

    public:
        void addSurface(Surface* surface) {
            surfaces_.push_back(surface);
            bbox_ += surface->bbox();
        }

        void runSimulation(unsigned long nParticles, double gridSize,
            int nThreads = -1);
};

#endif
