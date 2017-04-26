#include "simpleplasmamodel.h"
#include "simulationmodel.h"
#include "surface.h"
#include "surfaceemission.h"
#include "element_data.h"
#include "logger.h"

const double CONTINUOUS_FEED_RATE = 1e-6 / 3600.0 * NTP_VOLUME_TO_PARTICLES;

void doRun(size_t N_PARTICLES, double PULSE_LENGTH = 1.0e-3,
    double PARTICLES_PER_SECOND = CONTINUOUS_FEED_RATE,
    double SAMPLING_INTERVAL = 0.005, size_t N_TIME_SAMPLES = 0,
    uint_least32_t seed = 13122016) {

    // Set this to true to get detailed data about particle trajectories.
    // Useful for debugging with small particle numbers (< 100).
    const bool PARTICLE_LOOP_LOGGING = false;

    const double ELECTRON_DENSITY = 1e17;  // 1 / m^3
    const double ELECTRON_TEMPERATURE = 10e3;  // eV
    const Element ELEMENT = ARGON;
    const double ION_TEMPERATURE = 1.0;  // eV
    std::vector<double> ION_RELATIVE_DENSITIES {
        0.0729343, 0.04383817, 0.03905347, 0.03666112, 0.0586578, 0.05133851,
        0.06074164, 0.13743072, 0.1514723, 0.13927348, 0.10523982, 0.06349422,
        0.02817102, 0.00892281, 0.00243114, 0.00033945
    };
    const double GRID_SIZE = 0.001;  // metre
    const double CONFINEMENT_TIME = 1.9e-3;  // seconds
    const double EXTRACTION_EFFICIENCY = 0.1;  // Extraction efficiency factor, i.e. the pumping factor of the extraction surface
    const double CUTOFF_TIME = 5.0;  // A time at which the tracking of an individual particle is terminated

    // Create a vector of ion temperatures, the same value for every population
    // by default
    std::vector<double> ION_TEMPERATURES(ION_RELATIVE_DENSITIES.size(),
        ION_TEMPERATURE);
    // Create the plasma model, which is responsible for populating the
    // collision reaction data with the right parameters
    SimplePlasmaModel plasmamodel(
        "electron_data/electron_density_hiisi.csv",
        ELECTRON_DENSITY, ION_RELATIVE_DENSITIES, ELECTRON_TEMPERATURE,
        ION_TEMPERATURES, ELEMENT);

    // Surfaces are shared between various objects, so they are allocated by
    // make_shared
    SurfacePtr chamber = std::make_shared<Surface>(
        "real-model/chamber.stl",  // path to the STL
        0.0,  // pumping probability
        ROOM_TEMPERATURE_EV,  // surface temperature
        "chamber wall",  // surface label
        false,  // flip normals?
        0.001  // scaling factor. use 0.001 for STLs that are in millimeters
    );
    SurfacePtr surrounding_cylinder = std::make_shared<Surface>(
        "real-model/sheath.stl",
        1.0, ROOM_TEMPERATURE_EV, "surrounding cylinder", true, 0.001);
    SurfacePtr injection_surface = std::make_shared<Surface>(
        "real-model/injection_surface.stl",
        1.0, ROOM_TEMPERATURE_EV, "injection surface", true, 0.001);
    SurfacePtr extraction_surface = std::make_shared<Surface>(
        "real-model/extraction_surface.stl",
        EXTRACTION_EFFICIENCY, ROOM_TEMPERATURE_EV, "extraction surface", false,
        0.001);
    // The whole simulation geometry is defined as a surface collection. Add
    // all the loaded surfaces to it.
    SurfaceCollection surfaces;
    surfaces.addSurface(chamber);
    surfaces.addSurface(surrounding_cylinder);
    surfaces.addSurface(injection_surface);
    surfaces.addSurface(extraction_surface);
    // Create the simulation model object which is responsible for the control
    // flow.
    SimulationModel simModel(surfaces);
    // Add a particle source. Multiple sources can be defined, resulting in
    // different output files tagged with the label to be generated.
    simModel.addSource(std::make_unique<SurfaceEmission>(
        injection_surface,  // the surface object
        PARTICLES_PER_SECOND,  // particle emission rate
        ELEMENT,  // the element of the particles
        "injection"  // the label of the neutral source
    ));

    // Initialize a collision generator and populate the collision reactions
    CollisionGenerator generator(simModel.getGrid(GRID_SIZE));
    plasmamodel.populateCollisionReactions(generator,
        14032017,  // seed for Monte Carlo integration of the rate coefficients
        0.01  // Resolution of the rate coefficient LUT. Fraction of the maximal ion population temperature
    );
    // Initialize a neutralization generator and populate the neutralization
    // reactions
    NeutralizationGenerator ngenerator;
    plasmamodel.populateNeutralizationReactions(ngenerator,
        CONFINEMENT_TIME,
        "electron_data/cylinder_collision_points.csv",
        "electron_data/z1_collision_points.csv",
        "electron_data/z2_collision_points.csv",
        "rt.018.dat", surfaces, ELECTRON_DENSITY);

    std::string name = "test_" + std::to_string(N_PARTICLES);
    logger.setLogging(PARTICLE_LOOP_LOGGING);
    clock_t start_clock = clock();
    time_t start_time = time(NULL);
    simModel.runSimulation(generator, ngenerator,
        N_PARTICLES,  // number of test particles
        name,  // the file name prefix
        N_TIME_SAMPLES,  // number of sampled frames in the case of time-dependent simulation. 0 for stationary
        GRID_SIZE,  // grid size (metres)
        SAMPLING_INTERVAL,  // sampling interval (seconds)
        PULSE_LENGTH,  // the length of the injected gas pulse
        CUTOFF_TIME,
        -1,  // number of threads. -1 for automatic: 1 + 2 * hardware concurrency
        seed  // random seed
    );
    clock_t end_clock = clock();
    time_t end_time = time(NULL);

    std::cout << "System time: " << ((end_clock - start_clock) / CLOCKS_PER_SEC) << " sec\n";
    std::cout << "Wall clock time: " << (end_time - start_time) << " sec\n";

    logger.setLogging(true);
    surfaces.writeStatistics(logger);
    generator.writeStatistics(logger);
    ngenerator.writeStatistics(logger);
}

void run_stationary(size_t N_PARTICLES = 10000000,
    uint_least32_t seed = 1204146) {
    doRun(N_PARTICLES, 1.0e-3, CONTINUOUS_FEED_RATE, 0.0005, 0, seed);
}

void run_convergence_test() {
    std::vector<size_t> runs = {
        100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000
    };
    std::seed_seq sseq { 2017 };
    std::vector<uint_least32_t> seeds(runs.size());
    sseq.generate(seeds.begin(), seeds.end());

    for (size_t i = 0; i < runs.size(); ++i) {
        run_stationary(runs[i], seeds[i]);
    }
}

void run_time_dependent() {
    doRun(500000000, 1.0e-3, 3e12 / 1.0e-3, 0.0005, 20);
}

int main() {
    run_stationary(100000);
    //run_time_dependent();

    return EXIT_SUCCESS;
}

