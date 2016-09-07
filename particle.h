#ifndef PARTICLE_H
#define PARTICLE_H

#include "cgal_end_typedefs.h"

class Particle {
    Point position;
    Direction direction;
    double velocity;

    public:
        enum State { Free, OnWall, Colliding };
};

#endif
