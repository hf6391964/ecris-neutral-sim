#include "simulationmodel.h"

void SimulationModel::runSimulation(long nParticles, double timestep,
    int nThreads) {
    if (nThreads <= 0) {
        nThreads = std::thread::hardware_concurrency();
        if (nThreads == 0) {
            nThreads = 1;
        }
        nThreads = 2*nThreads + 1;
    }

    long particlesInChunk = nParticles / nThreads;

    for (Surface* pSurf : surfaces_) {
        if (!pSurf->isEmissive()) continue;

        Rng rng;

        for (long i = 0; i < nParticles; i++) {
            Particle particle;
        }
    }
}

