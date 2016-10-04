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
        unsigned long* count) const;

    void writeResults(std::string prefix, Vector* velocity,
        unsigned long* count, Grid& grid) const;

    public:
        void addSurface(Surface* surface) {
            surfaces_.push_back(surface);
            bbox_ += surface->bbox();
        }

        void runSimulation(unsigned long nParticles, double gridSize,
            std::string prefix, int nThreads = -1);
};

#endif
