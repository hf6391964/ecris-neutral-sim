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

struct Sample {
    Vector velocity = Vector(0.0, 0.0, 0.0);
    unsigned long count = 0;
};

typedef std::vector<Sample> SampleFrame;

class SimulationModel {
    std::vector<std::unique_ptr<NeutralSource>> sources_;
    SurfaceCollection &surfaces_;
    long nParticles_;
    Bbox bbox_;

    void simulationThread(
        const NeutralSource *pSource,
        const CollisionGenerator &collisionGenerator,
        const NeutralizationGenerator &neutralizationGenerator,
        std::mutex &writeMutex,
        unsigned long nParticles, double samplingInterval, long nTimeSamples,
        Grid grid, uint_least32_t seed, std::vector<SampleFrame> &frames,
        bool stationary, double cutoffTime) const;

    void writeResults(std::string prefix,
        const SampleFrame &frame,
        Grid grid, double t, std::string suffix) const;

    public:
        SimulationModel(SurfaceCollection &surfaces);

        Grid getGrid(double gridSize) const;
        void addSource(std::unique_ptr<NeutralSource> source);

        void runSimulation(
            const CollisionGenerator &collisionGenerator,
            const NeutralizationGenerator &neutralizationGenerator,
            unsigned long nParticles, std::string prefix,
            unsigned int nTimeSamples = 0,
            double gridSize = 0.1,
            double samplingInterval = 0.0005,
            double cutoffTime = 1.0,
            int nThreads = -1);
};
