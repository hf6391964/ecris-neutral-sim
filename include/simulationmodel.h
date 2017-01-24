#pragma once

#include <thread>
#include <random>
#include <iomanip>
#include <fstream>

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
        unsigned long nParticles, unsigned long maxSteps,
        double dt, const Grid &grid, uint_least32_t seed, Vector* velocity,
        unsigned long* count, bool stationary) const;

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
            bool stationary = true,
            double gridSize = 0.1, double maxTime = 0.1,
            double timestepFactor = 2.0, int nThreads = -1);
};
