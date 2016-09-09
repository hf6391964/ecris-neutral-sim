#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>

#include "randutils.hpp"
#include "cgal_and_typedefs.h"
#include "surface.h"

Rng rng;
int main() {
    Surface surf_walls("models/cylinder/cylinder_walls.stl", 0.1);
    Surface surf_walls1("models/cylinder/cylinder_walls.stl");

    for (int i = 0; i < 10; i++) {
        surf_walls.getRandomPoint(rng);
    }

    surf_walls.computeIntersection(Ray(Point(2.0, 0.1, 0.05), Point(0.0, 0.0, 0.0)));
    std::cout << std::endl;
    surf_walls1.computeIntersection(Ray(Point(2.0, 0.1, 0.05), Point(0.0, 0.0, 0.0)));

    return 0;
}

