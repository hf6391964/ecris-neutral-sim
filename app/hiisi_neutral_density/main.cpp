#include "simpleplasmamodel.h"
#include "simulationmodel.h"
#include "surface.h"
#include "surfaceemission.h"
#include "element_data.h"
#include "logger.h"

void doRun(size_t N_PARTICLES) {
    const bool PARTICLE_LOOP_LOGGING = false;
    const double ELECTRON_DENSITY = 1e17;
    const double ELECTRON_TEMPERATURE = 10e3;
    const Element ELEMENT = ARGON;
    const double ION_TEMPERATURE = 4.0;
    std::vector<double> ION_RELATIVE_DENSITIES {
        0.0729343, 0.04383817, 0.03905347, 0.03666112, 0.0586578, 0.05133851,
        0.06074164, 0.13743072, 0.1514723, 0.13927348, 0.10523982, 0.06349422,
        0.02817102, 0.00892281, 0.00243114, 0.00033945
    };
    const double GRID_SIZE = 0.002;
    const double SAMPLING_INTERVAL = 0.0005;
    const size_t N_TIME_SAMPLES = 0;  // 0 for stationary
    const double ESCAPE_TIME = 1e-5;

    std::vector<double> ION_TEMPERATURES(ION_RELATIVE_DENSITIES.size(),
        ION_TEMPERATURE);
    SimplePlasmaModel plasmamodel(
        "electron_density_hiisi.csv",
        ELECTRON_DENSITY, ION_RELATIVE_DENSITIES, ELECTRON_TEMPERATURE,
        ION_TEMPERATURES, ELEMENT);

    simthreadresources *thread_res = Util::allocateThreadResources(13122016);

    /*Surface radial_wall("model/radial_wall.stl",
        0.0, ROOM_TEMPERATURE_EV, "radial_wall", false);
    Surface end1("model/end1.stl",
        1.0, ROOM_TEMPERATURE_EV, "end1", true);
    Surface end2("model/end2.stl",
        1.0, ROOM_TEMPERATURE_EV, "end2", false);
    SurfaceCollection surfaces;
    surfaces.addSurface(&radial_wall);
    surfaces.addSurface(&end1);
    surfaces.addSurface(&end2);
    SurfaceEmission gasFeed(&radial_wall, 1e-12, ELEMENT);*/

    Surface chamber("real-model/stl-external/plasmannakemat2-trimmed.stl",
        0.0, ROOM_TEMPERATURE_EV, "chamber wall", false, 0.001);
    Surface surrounding_cylinder("real-model/surrounding_cylinder.stl",
        1.0, ROOM_TEMPERATURE_EV, "surrounding cylinder", true, 0.001);
    Surface injection_surface("real-model/injection_surface.stl",
        1.0, ROOM_TEMPERATURE_EV, "injection surface", false, 0.001);
    SurfaceCollection surfaces;
    surfaces.addSurface(&chamber);
    surfaces.addSurface(&surrounding_cylinder);
    surfaces.addSurface(&injection_surface);
    SurfaceEmission gasFeed(&injection_surface, 1e-12, ELEMENT);

    SimulationModel simModel(surfaces);
    simModel.addSource(&gasFeed);

    CollisionGenerator generator(simModel.getGrid(GRID_SIZE));
    NeutralizationGenerator ngenerator;
    plasmamodel.populateNeutralizationReactions(ngenerator, ESCAPE_TIME,
        "cylinder_collision_points.csv",
        "z1_collision_points.csv", "z2_collision_points.csv",
        "rt.018.dat", surfaces, ELECTRON_DENSITY);
    plasmamodel.populateCollisionReactions(generator, thread_res, 2.0);

    std::string name = "test" + std::to_string(N_PARTICLES);
    logger.setLogging(PARTICLE_LOOP_LOGGING);
    clock_t start_clock = clock();
    simModel.runSimulation(generator, ngenerator,
        N_PARTICLES, name, N_TIME_SAMPLES, GRID_SIZE, SAMPLING_INTERVAL);
    clock_t end_clock = clock();
    Util::deallocateThreadResources(thread_res);

    std::cout << "Time: " << ((end_clock - start_clock) / CLOCKS_PER_SEC) << " sec\n";

    logger.setLogging(true);
    surfaces.writeStatistics(logger);
    generator.writeStatistics(logger);
    ngenerator.writeStatistics(logger);
}

int main() {
    /*std::vector<size_t> runs = {
        1000000, 2000000, 5000000, 10000000, 20000000, 50000000,
        100000000, 200000000
    };

    for (size_t size : runs) {
        std::cout << "Running with " << size << " particles" << std::endl;
        doRun(size);
    }*/

    doRun(10000);

    return 0;
}

