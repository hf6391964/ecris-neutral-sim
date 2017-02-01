#include "simpleplasmamodel.h"
#include "simulationmodel.h"
#include "surface.h"
#include "surfaceemission.h"
#include "element_data.h"
#include "logger.h"

int main() {
    logger.setLogging(true);

    const double ELECTRON_DENSITY = 1e17;
    const double ELECTRON_TEMPERATURE = 10e3;
    const Element ELEMENT = ARGON;
    const double ION_TEMPERATURE = 4.0;
    std::vector<double> ION_RELATIVE_DENSITIES = {
        1.0/4096.0, 1.0/1024.0, 1.0/256.0, 1.0/64.0, 1.0/16.0, 1.0/8.0, 1.0/2.0,
        1.0/8.0, 1.0/16.0, 1.0/64.0, 1.0/256.0, 1.0/1024.0, 1.0/4096.0
    };
    const double GRID_SIZE = 0.005;
    const size_t N_PARTICLES = 20;//100000000;

    std::vector<double> ION_TEMPERATURES(ION_RELATIVE_DENSITIES.size(),
        ION_TEMPERATURE);
    SimplePlasmaModel plasmamodel(
        "electron_density_hiisi.csv",
        ELECTRON_DENSITY, ION_RELATIVE_DENSITIES, ELECTRON_TEMPERATURE,
        ION_TEMPERATURES, ELEMENT);

    simthreadresources *thread_res = Util::allocateThreadResources(13122016);

    Surface radial_wall("model/radial_wall.stl",
        0.0, ROOM_TEMPERATURE_EV, "radial_wall", false, 0.1);
    Surface end1("model/end1.stl",
        1.0, ROOM_TEMPERATURE_EV, "end1", true, 0.1);
    Surface end2("model/end2.stl",
        1.0, ROOM_TEMPERATURE_EV, "end2", false, 0.1);
    SurfaceCollection surfaces;
    surfaces.addSurface(&radial_wall);
    surfaces.addSurface(&end1);
    surfaces.addSurface(&end2);
    SurfaceEmission gasFeed(&radial_wall, 1e-12, ELEMENT);

    SimulationModel simModel(surfaces);
    simModel.addSource(&gasFeed);

    CollisionGenerator generator(simModel.getGrid(GRID_SIZE));
    NeutralizationGenerator ngenerator;
    plasmamodel.populateNeutralizationReactions(ngenerator, 1e-6,
        "cylinder_collision_points.csv",
        "z1_collision_points.csv", "z2_collision_points.csv",
        "rt.018.dat", surfaces, ELECTRON_DENSITY);
    plasmamodel.populateCollisionReactions(generator, thread_res, 0.1);

    simModel.runSimulation(generator, ngenerator,
        N_PARTICLES, "test", true, GRID_SIZE, 0.1, 2.0, 4);

    Util::deallocateThreadResources(thread_res);

    return 0;
}

