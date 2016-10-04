#ifndef PARTICLE_H
#define PARTICLE_H

#include "cgal_and_typedefs.h"
#include "surface.h"
#include "util.h"

/* #define DEBUG */

class Particle {
    public:
        Particle() {};

        enum State { Free, Pumped };

        const Point getPosition() const {
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
        }

        State getState() const {
            return state_;
        }

        double getSpeed() const {
            return speed_;
        }

        void setPosition(Point position) {
            position_ = position;
        }

        void setVelocity(double speed, Direction direction) {
            speed_ = speed;
            direction_ = direction;
        }

        Vector getVelocity() const {
            return speed_ * direction_.vector();
        }

        bool hasNextIntersection() {
            return nextIntersection_.pSurface != NULL;
        }

        double distanceToIntersection() {
            return std::sqrt(CGAL::squared_distance(position_,
                nextIntersection_.point));
        }

        void goToIntersection(Rng& rng);

        void goForward(double dt) {
            double stepLength = dt * speed_;
            position_ = position_ + stepLength * direction_.vector();
        }

        bool findNextIntersection(
            std::vector<Surface*>::const_iterator itSurface,
            std::vector<Surface*>::const_iterator itEnd);

    private:
        Point position_;
        Direction direction_;
        double speed_ = 0.0;  // m/s
        double temperature_ = 273.0;  // K
        double molarMass_ = 1.0;  // g/mol
        State state_ = Free;
        IntersectionPoint nextIntersection_;
};

#endif
