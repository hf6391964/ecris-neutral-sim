#pragma once

#include "cgal_and_typedefs.h"
#include "surfacecollection.h"
#include "util.h"
#include "constants.h"
#include "element_data.h"

/* #define DEBUG */

class Particle {
    public:
        enum State { Free, Pumped };

        Particle(Element element, double time);

        Element getElement() const;
        const Point getPosition() const;
        const Direction getDirection() const;
        const Ray getRay() const;
        void setMass_eV(double mass_eV);
        double getMass_eV() const;
        State getState() const;
        double getSpeed() const;
        double getTime() const;
        void setPosition(const Point &position);
        void setVelocity(double speed, const Direction &direction);
        void setVelocity(const Vector &vel);
        Vector getVelocity() const;
        bool hasNextIntersection();
        double distanceToIntersection();
        void goToIntersection(Rng& rng);
        void goForward(double dt);
        bool findNextIntersection(SurfaceCollection &surfaces);
        void setNextIntersection(const IntersectionPoint &ip);

    private:
        Point position_;
        Direction direction_;
        double speed_ = 0.0;  // m/s
        double mass_eV_ = 1.0;  // eV
        Element element_;
        double time_ = 0.0;
        const ElementData *elementData_;
        State state_ = Free;
        IntersectionPoint nextIntersection_;
};

