#include "simulationmodel.h"

void SimulationModel::runSimulation(unsigned long nParticles,
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
    double averageSpeed = Util::getMBAverage(273.0, 1.0);
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
    std::vector<Vector*> velocityPointers;
    std::vector<unsigned long*> countPointers;
    sseq.generate(seeds.begin(), seeds.end());
    for (int i = 0; i < nThreads; i++) {
        Vector* velocity = new Vector[arraySize];
        unsigned long* count = new unsigned long[arraySize];
        for (size_t k = 0; k < arraySize; k++) {
            velocity[k] = Vector(0.0, 0.0, 0.0);
            count[k] = 0;
        }
        velocityPointers.push_back(velocity);
        countPointers.push_back(count);
        threads.push_back(std::thread(&SimulationModel::simulationThread,
            this, particlesInChunk, nSteps, dt, grid, seeds[i],
            velocity, count, stationary));
    }

    for (auto it = threads.begin(); it != threads.end(); ++it) {
        it->join();
    }

    Vector* velocity = new Vector[arraySize];
    unsigned long* count = new unsigned long[arraySize];
    for (size_t j = 0; j < arraySize; j++) {
        Vector v(0.0, 0.0, 0.0);
        unsigned long c = 0;
        for (int n = 0; n < nThreads; n++) {
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

void SimulationModel::simulationThread(unsigned long nParticles,
    unsigned long maxSteps, double dt, Grid grid, uint_least32_t seed,
    Vector* velocity, unsigned long* count, bool stationary) const {
    size_t gridSize = grid.arraySize();

    for (NeutralSource* pSource : sources_) {
        Rng rng(seed);

        // TODO run simulation separately for each particle species in the
        // spectrum

        for (unsigned long i = 0; i < nParticles; i++) {
            Particle particle;

            pSource->generateParticle(particle, rng);

            if (!particle.findNextIntersection(surfaces_.begin(),
                surfaces_.end())) {
                /* std::cout << "Something wrong with intersections" << std::endl; */
                continue;
            }

            // Do this to spread the gas pulse evenly across the whole timestep
            // in time dependent simulations
            double timeRemainder = uni01(rng)*dt;
            for (unsigned long step = 0; step < maxSteps; step++) {

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
