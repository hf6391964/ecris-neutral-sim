#pragma once

#include <thread>
#include <random>
#include <iomanip>
#include <fstream>

#include "neutralsource.h"
#include "surface.h"
#include "particle.h"
#include "constants.h"
#include "argon_data.h"
#include "grid.h"
#include "util.h"
#include "collisiongenerator.h"

class SimulationModel {
    std::vector<NeutralSource*> sources_;
    std::vector<Surface*> surfaces_;
    long nParticles_;
    Bbox bbox_;

    void simulationThread(CollisionGenerator *collisionGenerator,
        unsigned long nParticles, unsigned long maxSteps,
        double dt, const Grid &grid, uint_least32_t seed, Vector* velocity,
        unsigned long* count, bool stationary) const;

    void writeResults(std::string prefix, Vector* velocity,
        unsigned long* count, Grid& grid, double t, std::string suffix) const;

    public:
        Grid getGrid(double gridSize) const {
            return Grid(bbox_, gridSize);
        }

        void addSource(NeutralSource* source) {
            sources_.push_back(source);
        }

        void addSurface(Surface* surface) {
            surfaces_.push_back(surface);
            bbox_ += surface->bbox();
        }

        void runSimulation(CollisionGenerator &collisionGenerator,
            unsigned long nParticles, std::string prefix,
            bool stationary = true,
            double gridSize = 0.1, double maxTime = 0.1,
            double timestepFactor = 2.0, int nThreads = -1);
};

