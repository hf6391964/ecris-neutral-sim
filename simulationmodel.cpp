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

            Point p;
            face_descriptor fd;
            std::tie(p, fd) = pSurf->getRandomPoint(rng);
            Direction d = pSurf->generateCosineLawDirection(fd, rng);
            double v = getMBSpeed(rng, particle.getTemperature(),
                particle.getMolarMass());
            particle.setVelocity(v, d);
        }
    }
}

double SimulationModel::getMBSpeed(Rng& rng, double T, double molarMass) {
    const int N_DIM = 3;
    const double a = std::sqrt(BOLTZMANN_CONSTANT * AVOGADRO_CONSTANT /
        (2.0 * 1e-3) * (T / molarMass));
    double v = 0.0;
    std::normal_distribution<double> d(0.0, a);
    for (int i = 0; i < N_DIM; i++) {
        double tmp = d(rng.engine());
        v += tmp*tmp;
    }
    return std::sqrt(v);
}
