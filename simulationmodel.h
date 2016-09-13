#ifndef SIMULATIONMODEL_H
#define SIMULATIONMODEL_H

#include <thread>

#include "surface.h"
#include "particle.h"
#include "constants.h"

class SimulationModel {
    std::vector<Surface*> surfaces_;
    long nParticles_;
    double timestep_;
    Bbox bbox_;

    public:
        void addSurface(Surface* surface) {
            surfaces_.push_back(surface);
            bbox_ += surface->bbox();
        }

        void runSimulation(long nParticles, double timestep, double gridSize,
            int nThreads = -1);

        static double getMBSpeed(Rng& rng, double T, double molarMass);
        static double getMBAverage(double T, double molarMass);
};

#endif
