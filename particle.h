#ifndef PARTICLE_H
#define PARTICLE_H

#include "cgal_and_typedefs.h"
#include "surface.h"
#include "util.h"

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

        bool hasNextIntersection() {
            return nextIntersection_.pSurface != NULL;
        }

        double distanceToIntersection() {
            return std::sqrt(CGAL::squared_distance(position_,
                nextIntersection_.point));
        }

        template<typename InputIterator>
        bool findNextIntersection(InputIterator itSurface,
            InputIterator itEnd) {
            Ray r = getRay();

            bool found = false;
            double nearestDistance = 0.0;
            Point* nearestPoint;
            face_descriptor nearestFaceId;
            Point* p = NULL;
            Surface* pSurface = NULL;
            nextIntersection_.pSurface = NULL;

            while (itSurface != itEnd) {
                std::list<Ray_intersection> intersections;
                (*itSurface)->computeIntersections(r,
                    std::back_inserter(intersections));

                for (Ray_intersection intersection : intersections) {
                    if (intersection &&
                        (p = boost::get<Point>(&(intersection->first)))) {
                        double distance = CGAL::squared_distance(position_, *p);
                        if (!found || distance < nearestDistance) {
                            found = true;
                            nearestPoint = p;
#ifdef DEBUG
                            std::cout << "Next intersection: ";
                            printPoint(*p);
                            std::cout << std::endl;
#endif
                            nearestDistance = distance;
                            nearestFaceId = intersection->second;
                            pSurface = *itSurface;
                        }
                    }
                }

                itSurface++;
            }

            if (found) {
                nextIntersection_.pSurface = pSurface;
                nextIntersection_.faceId = nearestFaceId;
                nextIntersection_.point = Point(
                    nearestPoint->x(),
                    nearestPoint->y(),
                    nearestPoint->z()
                );
            }

            return found;
        }

        void goToIntersection(Rng& rng) {
            if (nextIntersection_.pSurface == NULL) return;

            position_ = nextIntersection_.point;
            // TODO accommodation!
            speed_ = Util::getMBSpeed(rng,
                nextIntersection_.pSurface->getTemperature(), molarMass_);
            direction_ = nextIntersection_.pSurface->
                generateCosineLawDirection(nextIntersection_.faceId, rng);

            if (nextIntersection_.pSurface->checkIfPumped(rng)) {
                nextIntersection_.pSurface->addPumpedParticle();
                state_ = Pumped;
            }
        }

        void goForward(double dt) {
            double stepLength = dt * speed_;
            position_ = position_ + stepLength * direction_.vector();
        }

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
