#include "simulationmodel.h"
#include "logger.h"

SimulationModel::SimulationModel(SurfaceCollection &surfaces)
    : surfaces_(surfaces) {}

Grid SimulationModel::getGrid(double gridSize) const {
    return Grid(surfaces_.bbox(), gridSize);
}

void SimulationModel::addSource(NeutralSource* source) {
    sources_.push_back(source);
}

void SimulationModel::runSimulation(
    const CollisionGenerator &collisionGenerator,
    const NeutralizationGenerator &neutralizationGenerator,
    unsigned long nParticles,
    std::string prefix, unsigned int nTimeSamples, double gridSize,
    double samplingInterval, double cutoffTime, int nThreads) {
    if (logger.isLogging()) {
        nThreads = 1;
    } else if (nThreads <= 0) {
        nThreads = std::thread::hardware_concurrency() + 1;
    }

    bool stationary = nTimeSamples == 0;

    long particlesInChunk = nParticles / nThreads;

    Grid grid(surfaces_.bbox(), gridSize);

    size_t gSize = grid.arraySize();
    size_t arraySize = gSize;
    std::cout << "grid has " << gSize << " cells\n";
    std::cout << "Need " << ((sizeof(Vector) + sizeof(unsigned long)) * arraySize) << " bytes\n";

    std::cout << "Sampling interval: " << samplingInterval << std::endl;

    std::mutex writeMutex;
    std::seed_seq sseq = {2016, 9, 19};
    std::vector<std::thread> threads;
    std::vector<uint_least32_t> seeds(nThreads);
    sseq.generate(seeds.begin(), seeds.end());

    std::vector<SampleFrame> frames(stationary ? 1 : nTimeSamples,
        SampleFrame(arraySize, Sample()));

    std::cout << "Running with " << nThreads << " thread(s)..." << std::endl;
    // TODO extend and use SpatialDistribution here
    for (int i = 0; i < nThreads; ++i) {
        std::cout << "Starting thread " << i << std::endl;
        threads.push_back(std::thread(&SimulationModel::simulationThread,
            this, std::ref(collisionGenerator),
            std::ref(neutralizationGenerator),
            std::ref(writeMutex),
            particlesInChunk, samplingInterval, nTimeSamples, grid, seeds[i],
            std::ref(frames), stationary, cutoffTime));
    }

    for (auto &thread : threads) {
        thread.join();
    }

    for (auto &frame : frames) {
        for (auto &sample : frame) {
            if (sample.count > 0) {
                sample.velocity = sample.velocity / sample.count;
            } else {
                sample.velocity = Vector(0.0, 0.0, 0.0);
            }
        }
    }

    std::ofstream dim(prefix + "_dimensions.csv");
    grid.writeDimensions(dim);
    dim.close();

    if (stationary) {
        writeResults(prefix, frames[0], grid, 0.0, "stationary");
    } else {
        int nDigits = std::log10(nTimeSamples) + 1;
        for (unsigned long step = 0; step < nTimeSamples; ++step) {
            std::stringstream suffix;
            suffix << std::setw(nDigits) << std::setfill('0') << step;
            writeResults(prefix, frames[step],
                grid, step * samplingInterval, suffix.str());
        }
    }
}

void SimulationModel::simulationThread(
    const CollisionGenerator &collisionGenerator,
    const NeutralizationGenerator &neutralizationGenerator,
    std::mutex &writeMutex,
    unsigned long nParticles, double samplingInterval,
    long nTimeSamples, Grid grid, uint_least32_t seed,
    std::vector<SampleFrame> &frames,
    bool stationary, double cutoffTime) const {

    simthreadresources *thread_res = Util::allocateThreadResources(seed);

    const size_t gridSize = grid.arraySize();

    for (NeutralSource* pSource : sources_) {
        // TODO run simulation separately for each particle species in the
        // spectrum

        for (unsigned long i = 0; i < nParticles; ++i) {
            // Spread the gas pulse evenly across the whole sample interval
            // in time dependent simulations
            double initialTime = -uni01(thread_res->rng) * samplingInterval;
            Particle particle =
                pSource->generateParticle(thread_res->rng, initialTime);
            long nextSampleIndex = 0;

            logger << "\nParticle number " << i <<
                " generated with initial time = " << initialTime <<
                ", mass = " << particle.getMass_eV() << " eV\n";

            if (!particle.findNextIntersection(surfaces_, true)) {
                logger << "ERR: generated particle doesn't have intersections\n";
                continue;
            }

            IntersectionPoint ip = particle.getNextIntersection();
            if (ip.pSurface != NULL) {
                logger << "Next intersected surface: " <<
                    ip.pSurface->getLabel() << ", point is (" <<
                    ip.point.x() << ", " << ip.point.y() <<
                    ", " << ip.point.z() << ")\n";
            }

            while (particle.getTime() < cutoffTime &&
                   (stationary || nextSampleIndex < nTimeSamples)) {
                double isectDistance = particle.distanceToIntersection();
                double speed = particle.getSpeed();
                double timeToIntersection = isectDistance / speed;
                // This mean free time is determined by the reaction with
                // the greatest rate, so in principle it should suffice
                // for accurate enough sampling of collision reactions
                double meanFreeTime =
                    collisionGenerator.getMeanFreeTime(speed);

                double timeToSampleBoundary = nextSampleIndex * samplingInterval -
                    particle.getTime();

                // If temporally close to the sample boundary, ensure that we step
                // right there
                double timestep = std::min(meanFreeTime, timeToSampleBoundary +
                    0.001 * samplingInterval);

                logger << "Moving forward by step = " << timestep <<
                    ", isectDistance = " << isectDistance <<
                    ", speed = " << speed <<
                    ", mft = " << meanFreeTime <<
                    ", particle time = " << particle.getTime() << '\n';

                if (timeToIntersection < timestep) {
                    // Approximation: if time to intersection is smaller than
                    // the mean free time, just go directly to the intersection
                    particle.goToIntersection(thread_res->rng);
                    particle.findNextIntersection(surfaces_, true);

                    IntersectionPoint ip = particle.getNextIntersection();
                    if (ip.pSurface != NULL) {
                        logger << "Next intersected surface: " <<
                            ip.pSurface->getLabel() << ", point is (" <<
                            ip.point.x() << ", " << ip.point.y() <<
                            ", " << ip.point.z() << ")\n";
                    }
                } else {
                    particle.goForward(timestep);

                    Point p = particle.getPosition();
                    logger << "Stepped forward by " << timestep <<
                        ", position is (" << p.x() << ", " << p.y() <<
                        ", " << p.z() << ")\n";

                    CollisionReaction *reaction =
                        collisionGenerator.sampleCollision(thread_res->rng,
                            p, speed, timestep);

                    if (reaction != nullptr) {
                        CollisionProducts products =
                            reaction->computeReactionProducts(thread_res->rng,
                            p, particle);

                        std::vector<Particle> neutralProducts =
                            products.first;
                        unsigned int ionProducts = products.second;

                        logger << "sampled a collision reaction, label = " <<
                            reaction->getLabel() << ", " <<
                            neutralProducts.size() <<
                            " neutral products and " <<
                            ionProducts << " ionized products\n";

                        if (neutralProducts.size() > 0) {
                            particle = neutralProducts[0];
                        } else if (ionProducts == 1) {
                            // One ionized product produced, sample a neutralization reaction
                            particle =
                                neutralizationGenerator.sampleNeutralizationReaction(thread_res->rng, particle);
                        } else {
                            // There are no reaction products, we are
                            // done with this particle.
                            particle.setState(Particle::DestroyedInReaction);
                            logger << "No reaction products, particle marked as destroyed\n";
                        }

                        if (particle.getState() != Particle::DestroyedInReaction) {
                            // The particle might have transitioned in time, so
                            // adjust nextSampleIndex so that the next sampling
                            // will be aligned again to sample boundaries
                            // (otherwise the particle could get sampled immediately)
                            long sampleIndex = particle.getTime() / samplingInterval;
                            nextSampleIndex = sampleIndex + 1;

                            particle.findNextIntersection(surfaces_, false);
                            IntersectionPoint ip = particle.getNextIntersection();
                            if (ip.pSurface != NULL) {
                                logger << "Next intersected surface: " <<
                                    ip.pSurface->getLabel() << ", point is (" <<
                                    ip.point.x() << ", " << ip.point.y() <<
                                    ", " << ip.point.z() << ")\n";
                            }
                        }
                    }
                }

                Particle::State particleState = particle.getState();
                if (!particle.hasNextIntersection() ||
                    particleState != Particle::Free) {
                    std::string logStr;
                    if (particleState == Particle::Pumped) {
                        logStr = "Particle pumped out\n";
                    } else if (particleState == Particle::DestroyedInReaction) {
                        logStr = "Particle was destroyed in a reaction\n";
                    } else {
                        logStr = "Particle doesn't have intersections\n";
                    }
                    logger << logStr;
                    break;
                }

                // Check if we should sample
                long sampleIndex = Util::fastFloor(particle.getTime() / samplingInterval);
                if (sampleIndex >= nextSampleIndex) {
                    nextSampleIndex = sampleIndex + 1;
                    logger << "Sampled particle at t = " << particle.getTime() << '\n';
                    size_t arrIndex;
                    if (grid.arrayIndex(particle.getPosition(), arrIndex)) {
                        size_t frameIndex = stationary ? 0 : sampleIndex;
                        std::lock_guard<std::mutex> writeGuard(writeMutex);
                        Sample &sample = frames[frameIndex][arrIndex];
                        sample.count += 1;
                        sample.velocity = sample.velocity + particle.getVelocity();
                    }
                }
            }
        }
    }

    Util::deallocateThreadResources(thread_res);
}

void SimulationModel::writeResults(std::string prefix,
    const SampleFrame &frame,
    Grid grid, double t, std::string suffix) const {
    size_t nx, ny, nz;
    std::tie(nx, ny, nz) = grid.dimensions();
    std::stringstream filename;
    filename << prefix << '_' << suffix << ".csv";
    std::ofstream f(filename.str());
    f << "# vx" << CSV_SEP << "vy" << CSV_SEP << "vz" << CSV_SEP <<
        "count" << std::endl;
    f << "t = " << t << std::endl;

    size_t nxyz = nx*ny*nz;
    for (size_t j = 0; j < nxyz; ++j) {
        Vector v = frame[j].velocity;
        f << v.x() << CSV_SEP << v.y() << CSV_SEP << v.z() << CSV_SEP <<
            frame[j].count << '\n';
    }

    f.close();
}

