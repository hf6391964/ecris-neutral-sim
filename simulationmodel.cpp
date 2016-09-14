#include "simulationmodel.h"

void SimulationModel::runSimulation(long nParticles, double gridSize,
    int nThreads) {
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
    double averageSpeed = Util::getMBAverage(273.0, 44.0);
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

                if (particle.getState() == Particle::Pumped) {
                    break;
                }
            }
        }
    }
}

