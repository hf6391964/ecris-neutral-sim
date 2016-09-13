#include "simulationmodel.h"

void SimulationModel::runSimulation(long nParticles, double timestep,
    double gridSize, int nThreads) {
    if (nThreads <= 0) {
        nThreads = std::thread::hardware_concurrency();
        if (nThreads == 0) {
            nThreads = 1;
        }
        nThreads = 2*nThreads + 1;
    }

    long particlesInChunk = nParticles / nThreads;

    double xlen = bbox_.xmax() - bbox_.xmin(),
           ylen = bbox_.ymax() - bbox_.ymin(),
           zlen = bbox_.zmax() - bbox_.zmin();

    Bbox simBbox = bbox_ + Point(
        bbox_.xmin() + std::ceil(xlen / gridSize) * gridSize,
        bbox_.ymin() + std::ceil(ylen / gridSize) * gridSize,
        bbox_.zmin() + std::ceil(zlen / gridSize) * gridSize
    ).bbox();

    double maxTime = 50.0;
    double timeStepFactor = 2.0;
    // TODO figure out typical particle parameters (e.g. average of present
    // particle masses and surface temperatures)
    double averageSpeed = getMBAverage(273.0, 44.0);
    double dt = timeStepFactor * gridSize / averageSpeed;
    unsigned long maxSteps = maxTime / dt;

    for (Surface* pSurf : surfaces_) {
        if (!pSurf->isEmissive()) continue;

        Rng rng;

        // TODO run simulation separately for each particle species in the
        // spectrum

        for (long i = 0; i < nParticles; i++) {
            Particle particle;
            Point p;
            face_descriptor fd;
            std::tie(p, fd) = pSurf->getRandomPoint(rng);
            Direction d = pSurf->generateCosineLawDirection(fd, rng);
            double v = getMBSpeed(rng, particle.getTemperature(),
                particle.getMolarMass());
            particle.setVelocity(v, d);

            IntersectionPoint ip;
            double distance = 0.0;
            for (Surface* pTargetSurf : surfaces_) {
                IntersectionPoint curIp;
                double curDistance;
                if (pTargetSurf->computeIntersection(particle.getRay(), curIp,
                    curDistance)) {
                    if (ip.pSurface == NULL || curDistance < distance) {
                        distance = curDistance;
                        ip.pSurface = pTargetSurf;
                        ip.point = curIp.point;
                        ip.faceId = curIp.faceId;
                    }
                }
            }

            unsigned long step = 0;
            while (step < maxSteps) {
            }
        }
    }
}

double SimulationModel::getMBSpeed(Rng& rng, double T, double molarMass) {
    // Maxwell-Boltzmann distribution implemented by drawing the three velocity
    // components from normal distribution and taking the resulting magnitude
    const int N_DIM = 3;
    const double a = std::sqrt(GAS_CONSTANT * (T / molarMass));
    std::normal_distribution<double> d(0.0, a);

    double v = 0.0;
    for (int i = 0; i < N_DIM; i++) {
        double tmp = d(rng.engine());
        v += tmp*tmp;
    }
    return std::sqrt(v);
}

double SimulationModel::getMBAverage(double T, double molarMass) {
    return std::sqrt(8.0 * GAS_CONSTANT / M_PI * (T / molarMass));
}
