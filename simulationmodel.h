#ifndef SIMULATIONMODEL_H
#define SIMULATIONMODEL_H

#include <thread>
#include <random>
#include <iomanip>
#include <fstream>

#include "surface.h"
#include "particle.h"
#include "constants.h"
#include "grid.h"
#include "util.h"

class SimulationModel {
    std::vector<Surface*> surfaces_;
    long nParticles_;
    Bbox bbox_;

    void simulationThread(unsigned long nParticles, unsigned long maxSteps,
        double dt, Grid grid, uint_least32_t seed, Vector* velocity,
        unsigned long* count, bool stationary) const;

    void writeResults(std::string prefix, Vector* velocity,
        unsigned long* count, Grid& grid, double t, std::string suffix) const;

    public:
        void addSurface(Surface* surface) {
            surfaces_.push_back(surface);
            bbox_ += surface->bbox();
        }

        void runSimulation(unsigned long nParticles, std::string prefix,
            bool stationary = true,
            double gridSize = 0.1, double maxTime = 0.1,
            double timestepFactor = 2.0, int nThreads = -1);
};

#endif
