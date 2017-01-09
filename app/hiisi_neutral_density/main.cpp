#include "simpleplasmamodel.h"
#include "simulationmodel.h"
#include "surface.h"
#include "surfaceemission.h"
#include "argon_data.h"

int main() {
    const double ELECTRON_DENSITY = 1e17;
    const double ELECTRON_TEMPERATURE = 10e3;
    const ElementData ELEM_DATA = ARGON_DATA;
    const double ION_TEMPERATURE = 4.0;
    std::vector<double> ION_RELATIVE_DENSITIES = {
        1.0/4096.0, 1.0/1024.0, 1.0/256.0, 1.0/64.0, 1.0/16.0, 1.0/8.0, 1.0/2.0, 1.0/8.0, 1.0/16.0, 1.0/64.0, 1.0/256.0, 1.0/1024.0, 1.0/4096.0
    };
    const double GRID_SIZE = 0.005;
    const size_t N_PARTICLES = 100000000;

    std::vector<double> ION_TEMPERATURES(ION_RELATIVE_DENSITIES.size(),
        ION_TEMPERATURE);
    SimplePlasmaModel plasmamodel(
        "electron_density_hiisi.csv",
        ELECTRON_DENSITY, ION_RELATIVE_DENSITIES, ELECTRON_TEMPERATURE,
        ION_TEMPERATURES, ELEM_DATA);

    simthreadresources *thread_res = Util::allocateThreadResources(13122016);
    SimulationModel simModel;

    Surface radial_wall("model/radial_wall.stl",
        0.0, ROOM_TEMPERATURE_EV, false, 0.1);
    Surface end1("model/end1.stl",
        1.0, ROOM_TEMPERATURE_EV, true, 0.1);
    Surface end2("model/end2.stl",
        1.0, ROOM_TEMPERATURE_EV, false, 0.1);

    simModel.addSurface(&radial_wall);
    simModel.addSurface(&end1);
    simModel.addSurface(&end2);
    SurfaceEmission gasFeed(&radial_wall, 1e-12);
    simModel.addSource(&gasFeed);

    CollisionGenerator generator(simModel.getGrid(GRID_SIZE));
    plasmamodel.populateCollisionReactions(generator, thread_res);

    simModel.runSimulation(generator, N_PARTICLES, "test", true, GRID_SIZE, 0.1, 2.0, 4);

    Util::deallocateThreadResources(thread_res);

    return 0;
}

