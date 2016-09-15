#include <iostream>

#include "cgal_and_typedefs.h"
#include "surface.h"
#include "simulationmodel.h"

int main() {
    SimulationModel model;

    Surface* surf_walls = new Surface("models/cylinder/cylinder_walls.stl",
        0.0, 273.0, 0.0, 0.1);
    Surface* surf_source = new Surface("models/cylinder/cylinder_cap1.stl",
        1.0, 273.0, 1e-12, 0.1);
    Surface* surf_sink = new Surface("models/cylinder/cylinder_cap2.stl",
        1.0, 273.0, 0.0, 0.1);

    model.addSurface(surf_walls);
    model.addSurface(surf_source);
    model.addSurface(surf_sink);

    unsigned long nParticles = 1000000;
    model.runSimulation(nParticles, 1e-1);

    std::cout << "% of particles pumped to source: " <<
        (100 * surf_source->getPumpedParticles() / nParticles) << std::endl;
    std::cout << "% of particles pumped to sink: " <<
        (100 * surf_sink->getPumpedParticles() / nParticles) << std::endl;

    delete surf_walls;
    delete surf_source;
    delete surf_sink;

    return 0;
}

