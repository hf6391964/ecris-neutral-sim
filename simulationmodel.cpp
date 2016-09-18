#include "simulationmodel.h"

void SimulationModel::runSimulation(unsigned long nParticles, double gridSize,
    int nThreads) {
    if (nThreads <= 0) {
        nThreads = std::thread::hardware_concurrency();
        if (nThreads == 0) {
            nThreads = 1;
        }
        nThreads = 2*nThreads + 1;
    }

    long particlesInChunk = nParticles / nThreads;

    double maxTime = 50.0;
    double timeStepFactor = 2.0;
    // TODO figure out typical particle parameters (e.g. average of present
    // particle masses and surface temperatures)
    double averageSpeed = Util::getMBAverage(273.0, 1.0);
    double dt = timeStepFactor * gridSize / averageSpeed;
    unsigned long maxSteps = maxTime / dt;

    Grid grid(bbox_, gridSize);

    size_t arraySize = (nThreads + 1) * grid.arraySize();
    Vector* velocity = new Vector[arraySize];
    unsigned long* count = new unsigned long[arraySize];

    std::cout << "Timestep: " << dt << std::endl;
    std::cout << "Running with " << nThreads << " thread(s)..." << std::endl;

    std::seed_seq sseq = {2016, 9, 19};
    std::vector<std::thread> threads;
    std::vector<uint_least32_t> seeds(nThreads);
    sseq.generate(seeds.begin(), seeds.end());
    for (int i = 0; i < nThreads; i++) {
        threads.push_back(std::thread(&SimulationModel::simulationThread,
            this, particlesInChunk, maxSteps, dt, grid, seeds[i]));
    }

    for (auto it = threads.begin(); it != threads.end(); ++it) {
        it->join();
    }

    delete[] velocity;
    delete[] count;
}

void SimulationModel::simulationThread(unsigned long nParticles,
    unsigned long maxSteps, double dt, Grid grid, uint_least32_t seed) {
    for (Surface* pSurf : surfaces_) {
        if (!pSurf->isEmissive()) continue;

        Rng rng(seed);

        // TODO run simulation separately for each particle species in the
        // spectrum

        for (unsigned long i = 0; i < nParticles; i++) {
            Particle particle;
            Point p;
            face_descriptor fd;

            std::tie(p, fd) = pSurf->getRandomPoint(rng);
            Direction d = pSurf->generateCosineLawDirection(fd, rng);
            double v = Util::getMBSpeed(rng, particle.getTemperature(),
                particle.getMolarMass());
            particle.setPosition(p);
            particle.setVelocity(v, d);

            if (!particle.findNextIntersection(surfaces_.begin(),
                surfaces_.end())) {
                std::cout << "Something wrong with intersections" << std::endl;
                continue;
            }

            for (unsigned long step = 0; step < maxSteps; step++) {
                double timeRemainder = dt;

                bool moving = true;
                while (moving) {
                    if (!particle.hasNextIntersection() ||
                        particle.getState() == Particle::Pumped) {
                        moving = false;
                    } else {
                        double isectDistance = particle.distanceToIntersection();
                        double speed = particle.getSpeed();
                        if (isectDistance <= speed*timeRemainder) {
                            timeRemainder -= isectDistance / particle.getSpeed();
                            particle.goToIntersection(rng);
                            particle.findNextIntersection(surfaces_.begin(),
                                surfaces_.end());
                        } else {
                            moving = false;
                            particle.goForward(timeRemainder);
                        }
                    }
                }

                if (!particle.hasNextIntersection() ||
                    particle.getState() == Particle::Pumped) {
                    break;
                }
            }
        }
    }
}
