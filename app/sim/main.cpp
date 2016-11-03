#include <iostream>

#include "cgal_and_typedefs.h"
#include "surface.h"
#include "simulationmodel.h"
#include "surfaceemission.h"

int main() {
    SimulationModel model;

    Surface surf_walls("models/elbow/elbow_3_4_body.stl", 0.0, ROOM_TEMPERATURE_EV, true, 0.1);
    Surface surf_source("models/elbow/elbow_3_4_cap1.stl", 1.0, ROOM_TEMPERATURE_EV, false, 0.1);
    Surface surf_sink("models/elbow/elbow_3_4_cap2.stl", 1.0, ROOM_TEMPERATURE_EV, false, 0.1);

    model.addSurface(&surf_walls);
    model.addSurface(&surf_source);
    model.addSurface(&surf_sink);

    SurfaceEmission gasFeed(&surf_source, 1e-12);
    model.addSource(&gasFeed);

    unsigned long nParticles = 100000;
    model.runSimulation(nParticles, "test", true, 0.1, 5);

    std::cout << "% of particles pumped to source: " <<
        (100 * surf_source.getPumpedParticles() / nParticles) << std::endl;
    std::cout << "% of particles pumped to sink: " <<
        (100 * surf_sink.getPumpedParticles() / nParticles) << std::endl;

    return 0;
}

