#include "simulationmodel.h"

void SimulationModel::runSimulation(CollisionGenerator &collisionGenerator,
    unsigned long nParticles,
    std::string prefix, bool stationary, double gridSize, double maxTime,
    double timestepFactor, int nThreads) {
    if (nThreads <= 0) {
        nThreads = std::thread::hardware_concurrency();
        if (nThreads == 0) {
            nThreads = 1;
        }
        nThreads = 2*nThreads + 1;
    }

    long particlesInChunk = nParticles / nThreads;

    // TODO figure out typical particle parameters (e.g. average of present
    // particle masses and surface temperatures)
    //
    // There will be two separate choices of timestep length:
    // the density sampling timestep (long, to reduce correlation between
    // sampled frames) and the particle simulation timestep (determined in the
    // particle loop based on reaction rate majorant estimate)
    double averageSpeed = Util::getMBAverage(ROOM_TEMPERATURE_EV,
        ARGON_DATA.mass * ATOMIC_MASS_TO_EV);
    double dt = timestepFactor * gridSize / averageSpeed;
    unsigned long nSteps = maxTime / dt;
    std::cout << "max steps: " << nSteps << std::endl;

    Grid grid(bbox_, gridSize);

    size_t gSize = grid.arraySize();
    size_t arraySize = gSize;
    if (!stationary) {
        arraySize *= nSteps;
    }

    std::cout << "Timestep: " << dt << std::endl;
    std::cout << "Running with " << nThreads << " thread(s)..." << std::endl;

    std::seed_seq sseq = {2016, 9, 19};
    std::vector<std::thread> threads;
    std::vector<uint_least32_t> seeds(nThreads);
    // TODO extend and use SpatialDistribution here
    std::vector<Vector*> velocityPointers;
    std::vector<unsigned long*> countPointers;
    sseq.generate(seeds.begin(), seeds.end());
    for (int i = 0; i < nThreads; ++i) {
        Vector* velocity = new Vector[arraySize];
        unsigned long* count = new unsigned long[arraySize];
        for (size_t k = 0; k < arraySize; ++k) {
            velocity[k] = Vector(0.0, 0.0, 0.0);
            count[k] = 0;
        }
        velocityPointers.push_back(velocity);
        countPointers.push_back(count);
        threads.push_back(std::thread(&SimulationModel::simulationThread,
            this, &collisionGenerator,
            particlesInChunk, nSteps, dt, grid, seeds[i],
            velocity, count, stationary));
    }

    for (auto it = threads.begin(); it != threads.end(); ++it) {
        it->join();
    }

    Vector* velocity = new Vector[arraySize];
    unsigned long* count = new unsigned long[arraySize];
    for (size_t j = 0; j < arraySize; ++j) {
        Vector v(0.0, 0.0, 0.0);
        unsigned long c = 0;
        for (int n = 0; n < nThreads; ++n) {
            v = v + velocityPointers[n][j];
            c += countPointers[n][j];
        }
        if (c > 0) {
            velocity[j] = v / c;
        }
        count[j] = c;
    }

    std::ofstream dim(prefix + "_dimensions.csv");
    grid.writeDimensions(dim);
    dim.close();

    if (stationary) {
        writeResults(prefix, velocity, count, grid, 0.0, "stationary");
    } else {
        int nDigits = std::log10(nSteps) + 1;
        for (unsigned long step = 0; step < nSteps; step++) {
            std::stringstream suffix;
            suffix << std::setw(nDigits) << std::setfill('0') << step;
            writeResults(prefix, velocity + step*gSize, count + step*gSize,
                grid, step * dt, suffix.str());
        }
    }

    for (auto p : velocityPointers) {
        delete[] p;
    }
    for (auto p : countPointers) {
        delete[] p;
    }
    delete[] velocity;
    delete[] count;
}

void SimulationModel::simulationThread(CollisionGenerator *collisionGenerator,
    unsigned long nParticles,
    unsigned long maxSteps, double dt, const Grid &grid, uint_least32_t seed,
    Vector* velocity, unsigned long* count, bool stationary) const {
    simthreadresources *thread_res = Util::allocateThreadResources(seed);

    size_t gridSize = grid.arraySize();

    for (NeutralSource* pSource : sources_) {
        // TODO run simulation separately for each particle species in the
        // spectrum

        for (unsigned long i = 0; i < nParticles; ++i) {
            Particle particle;

            pSource->generateParticle(particle, thread_res->rng);

            if (!particle.findNextIntersection(surfaces_.begin(),
                surfaces_.end())) {
                std::cout << "Something wrong with face normals / intersections"
                    << std::endl;
                continue;
            }

            bool destroyedInReaction = false;

            // Do this to spread the gas pulse evenly across the whole timestep
            // in time dependent simulations
            double timeRemainder = uni01(thread_res->rng)*dt;
            for (unsigned long step = 0; step < maxSteps; ++step) {
                bool moving = true;

                while (moving) {
                    if (!particle.hasNextIntersection() ||
                        particle.getState() == Particle::Pumped) {
                        moving = false;
                        destroyedInReaction = true;
                    } else {
                        double isectDistance = particle.distanceToIntersection();
                        double speed = particle.getSpeed();
                        double meanFreeTime =
                            collisionGenerator->getMeanFreeTime(speed);

                        double timestep = std::min(timeRemainder, meanFreeTime);

                        if (isectDistance <= speed*timestep) {
                            // TODO this might need some rethinking. If time to
                            // wall intersection is less than mean free time,
                            // currently no collision reactions will be sampled
                            // but the particle will be directly brought to the
                            // intersection.
                            timeRemainder -= isectDistance / speed;
                            particle.goToIntersection(thread_res->rng);
                            particle.findNextIntersection(surfaces_.begin(),
                                surfaces_.end());
                        } else {
                            particle.goForward(timestep);
                            timeRemainder -= timestep;
                            CollisionReaction *reaction =
                                collisionGenerator->sampleCollision(
                                    thread_res->rng, particle.getPosition(),
                                    speed, timestep);
                            if (reaction != NULL) {
                                std::vector<Particle> products =
                                    reaction->computeReactionProducts(thread_res->rng,
                                    particle.getPosition(), particle);

                                // TODO multiple reaction products?
                                if (products.size() >= 1) {
                                    particle = products[0];
                                    particle.findNextIntersection(surfaces_.begin(),
                                        surfaces_.end());
                                } else {
                                    // There are no reaction products, we are
                                    // done with this particle.
                                    destroyedInReaction = true;
                                }
                            }
                        }

                        if (timeRemainder <= 0.0) {
                            moving = false;
                        }
                    }
                }

                if (destroyedInReaction) {
                    break;
                }

                size_t arrIndex;
                if (grid.arrayIndex(particle.getPosition(), arrIndex)) {
                    if (!stationary) {
                        arrIndex += step*gridSize;
                    }
                    count[arrIndex] += 1;
                    velocity[arrIndex] = velocity[arrIndex] +
                        particle.getVelocity();
                }

                if (!particle.hasNextIntersection() ||
                    particle.getState() == Particle::Pumped) {
                    break;
                }

                timeRemainder = dt;
            }
        }
    }

    Util::deallocateThreadResources(thread_res);
}

void SimulationModel::writeResults(std::string prefix, Vector* velocity,
    unsigned long* count, Grid& grid, double t, std::string suffix) const {
    size_t nx, ny, nz;
    std::tie(nx, ny, nz) = grid.dimensions();
    std::stringstream filename;
    filename << prefix << '_' << suffix << ".csv";
    std::ofstream f(filename.str());
    f << "# vx" << CSV_SEP << "vy" << CSV_SEP << "vz" << CSV_SEP <<
        "count" << std::endl;
    f << "t = " << t << std::endl;

    size_t nxyz = nx*ny*nz;
    for (size_t j = 0; j < nxyz; j++) {
        Vector v = velocity[j];
        f << v.x() << CSV_SEP << v.y() << CSV_SEP << v.z() << CSV_SEP <<
            count[j] << std::endl;
    }

    f.close();
}

