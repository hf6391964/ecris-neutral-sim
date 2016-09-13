#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>

#include "randutils.hpp"
#include "cgal_and_typedefs.h"
#include "surface.h"
#include "simulationmodel.h"

Rng rng;
int main() {
    SimulationModel model;

    Surface* surf_walls = new Surface("models/cylinder/cylinder_walls.stl",
        0.0, 273.0, 0.0, 0.1);
    Surface* surf_source = new Surface("models/cylinder/cylinder_cap1.stl",
        0.0, 273.0, 1e-12, 0.1);
    Surface* surf_sink = new Surface("models/cylinder/cylinder_cap2.stl",
        1.0, 273.0, 0.0, 0.1);

    model.addSurface(surf_walls);
    model.addSurface(surf_source);
    model.addSurface(surf_sink);

    model.runSimulation(100, 1e-5);

    delete surf_walls;
    delete surf_source;
    delete surf_sink;

    return 0;
}

