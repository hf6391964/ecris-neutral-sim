#ifndef PARTICLE_H
#define PARTICLE_H

#include "cgal_and_typedefs.h"

class Particle {
    Point position;
    Direction direction;
    double speed;

    public:
        enum State { Free, OnWall, Colliding };

        void updateVelocity(double speed, Direction direction);
};

#endif
