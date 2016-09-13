#ifndef PARTICLE_H
#define PARTICLE_H

#include "cgal_and_typedefs.h"
#include "surface.h"

class Particle {
    public:
        Particle() {};

        enum State { Free };

        void setVelocity(double speed, Direction direction);

        void setNextIntersection(IntersectionPoint& ip);

        const Point& getPosition() const {
            return position_;
        }

        const Direction& getDirection() const {
            return direction_;
        }

        const Ray getRay() const {
            return Ray(position_, direction_);
        }

        double getMolarMass() const {
            return molarMass_;
        }

        double getTemperature() const {
            return temperature_;
        };

    private:
        Point position_;
        Direction direction_;
        double speed_ = 0.0;  // m/s
        double temperature_ = 0.0;  // K
        double molarMass_ = 1.0;  // g/mol
        State state_ = Free;
        IntersectionPoint nextIntersection_;
};

#endif
