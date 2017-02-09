#pragma once

#include "cgal_and_typedefs.h"
#include "surfacecollection.h"
#include "util.h"
#include "constants.h"
#include "element_data.h"


class Particle {
    public:
        enum State { Free, Pumped, DestroyedInReaction };

        Particle(Element element, double time);

        Element getElement() const;
        const Point getPosition() const;
        const Direction getDirection() const;
        const Ray getRay() const;
        double getMass_eV() const;
        State getState() const;
        void setState(State state);
        double getSpeed() const;
        double getTime() const;
        void setTime(double time);
        void setPosition(const Point &position);
        void setVelocity(double speed, const Direction &direction);
        void setVelocity(const Vector &vel);
        Vector getVelocity() const;
        bool hasNextIntersection() const;
        double distanceToIntersection() const;
        void goToIntersection(Rng& rng);
        void goForward(double dt);
        bool findNextIntersection(SurfaceCollection &surfaces);
        void setNextIntersection(const IntersectionPoint &ip);
        IntersectionPoint getNextIntersection() const;

    private:
        Point position_;
        Direction direction_;
        double speed_ = 0.0;  // m/s
        Element element_;
        const ElementData *elementData_;
        double mass_eV_;  // eV
        double time_ = 0.0;
        State state_ = Free;
        IntersectionPoint nextIntersection_;
};
