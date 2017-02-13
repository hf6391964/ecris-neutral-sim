#pragma once

#include <thread>
#include <random>
#include <iomanip>
#include <fstream>
#include <mutex>

#include "neutralsource.h"
#include "surfacecollection.h"
#include "particle.h"
#include "constants.h"
#include "element_data.h"
#include "grid.h"
#include "util.h"
#include "collisiongenerator.h"
#include "neutralizationgenerator.h"

class SimulationModel {
    std::vector<NeutralSource*> sources_;
    SurfaceCollection &surfaces_;
    long nParticles_;
    Bbox bbox_;

    void simulationThread(
        CollisionGenerator *collisionGenerator,
        NeutralizationGenerator *neutralizationGenerator,
        std::mutex &writeMutex,
        unsigned long nParticles, double samplingInterval, long nTimeSamples,
        const Grid &grid, uint_least32_t seed, Vector* velocity,
        unsigned long* count, bool stationary, double cutoffTime) const;

    void writeResults(std::string prefix, Vector* velocity,
        unsigned long* count, Grid& grid, double t, std::string suffix) const;

    public:
        SimulationModel(SurfaceCollection &surfaces);

        Grid getGrid(double gridSize) const;
        void addSource(NeutralSource* source);

        void runSimulation(
            CollisionGenerator &collisionGenerator,
            NeutralizationGenerator &neutralizationGenerator,
            unsigned long nParticles, std::string prefix,
            unsigned int nTimeSamples = 0,
            double gridSize = 0.1,
            double samplingInterval = 0.0005,
            double cutoffTime = 1.0,
            int nThreads = -1);
};
