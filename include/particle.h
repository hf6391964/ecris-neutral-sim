#pragma once

#include "cgal_and_typedefs.h"
#include "surface.h"
#include "util.h"
#include "constants.h"
#include "element_data.h"

/* #define DEBUG */

class Particle {
    public:
        enum State { Free, Pumped };

        Particle(Element element);

        Element getElement() const;
        const Point getPosition() const;
        const Direction getDirection() const;
        const Ray getRay() const;
        void setMass_eV(double mass_eV);
        double getMass_eV() const;
        State getState() const;
        double getSpeed() const;
        void setPosition(const Point &position);
        void setVelocity(double speed, const Direction &direction);
        void setVelocity(const Vector &vel);
        Vector getVelocity() const;
        bool hasNextIntersection();
        double distanceToIntersection();
        void goToIntersection(Rng& rng);
        void goForward(double dt);
        bool findNextIntersection(
            std::vector<Surface*>::const_iterator itSurface,
            std::vector<Surface*>::const_iterator itEnd);

    private:
        Point position_;
        Direction direction_;
        double speed_ = 0.0;  // m/s
        double mass_eV_ = 1.0;  // eV
        Element element_;
        const ElementData *elementData_;
        State state_ = Free;
        IntersectionPoint nextIntersection_;
};
