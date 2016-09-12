#ifndef PARTICLE_H
#define PARTICLE_H

#include "cgal_and_typedefs.h"
#include "surface.h"

class Particle {
    public:
        Particle() {};

        enum State { Free };

        void setVelocity(double speed, Direction direction);

    private:
        Point position_;
        Direction direction_;
        double speed_ = 0.0;
        double temperature_ = 0.0;
        State state_ = Free;
};

#endif
