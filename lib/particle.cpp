#include "particle.h"

#define USE_ALL_INTERSECTIONS

Particle::Particle(Element element) : element_(element) {
    elementData_ = ELEMENT_DATA.at(element_);
}

Element Particle::getElement() const {
    return element_;
}

const Point Particle::getPosition() const {
    return position_;
}

const Direction Particle::getDirection() const {
    return direction_;
}

const Ray Particle::getRay() const {
    return Ray(position_, direction_);
}

void Particle::setMass_eV(double mass_eV) {
    mass_eV_ = mass_eV;
}

double Particle::getMass_eV() const {
    return mass_eV_;
}

Particle::State Particle::getState() const {
    return state_;
}

double Particle::getSpeed() const {
    return speed_;
}

void Particle::setPosition(const Point &position) {
    position_ = position;
}

void Particle::setVelocity(double speed, const Direction &direction) {
    speed_ = speed;
    direction_ = direction;
}

void Particle::setVelocity(const Vector &vel) {
    direction_ = Direction(vel);
    speed_ = std::sqrt(vel.squared_length());
}

Vector Particle::getVelocity() const {
    return speed_ * direction_.vector();
}

bool Particle::hasNextIntersection() {
    return nextIntersection_.pSurface != NULL;
}

double Particle::distanceToIntersection() {
    return std::sqrt(CGAL::squared_distance(position_,
        nextIntersection_.point));
}

void Particle::goForward(double dt) {
    double stepLength = dt * speed_;
    position_ = position_ + stepLength * direction_.vector();
}

void Particle::goToIntersection(Rng& rng) {
    if (nextIntersection_.pSurface == NULL) return;

    position_ = nextIntersection_.point;

    // Take thermal accommodation in account

    // 1. Calculate a new particle speed in the case of full thermalization
    // with surface
    double thermalizedSpeed = Util::getMBSpeed(rng,
        nextIntersection_.pSurface->getTemperature(), mass_eV_);
    // Energy is quadratic in speed
    // Scaling by 1/2 * m is ignored here
    double thermalizedEnergy = thermalizedSpeed*thermalizedSpeed;
    double currentEnergy = speed_*speed_;
    // Query accommodation coefficient from the surface
    double alpha =
        nextIntersection_.pSurface->getAccommodationCoefficient(element_);
    // Accommodation coefficient is defined according to
    // http://prod.sandia.gov/techlib/access-control.cgi/2005/056084.pdf
    //
    // alpha = (E_in - E_re) / (E_in - E_w)
    // where E_in is incident energy flux, E_re reflected energy flux and
    // E_w energy flux in case of full thermalization with surface. Thus
    //
    // E_re = E_in - alpha * (E_in - E_w)
    double newEnergy = currentEnergy -
        alpha * (currentEnergy - thermalizedEnergy);
    // Speed is the square root of energy
    speed_ = std::sqrt(newEnergy);

    // Cosine law emission
    direction_ = nextIntersection_.pSurface->
        generateCosineLawDirection(nextIntersection_.faceId, rng);

    if (nextIntersection_.pSurface->checkIfPumped(rng)) {
        nextIntersection_.pSurface->addPumpedParticle();
        state_ = Pumped;
    }
}

bool Particle::findNextIntersection(
    std::vector<Surface*>::const_iterator itSurface,
    std::vector<Surface*>::const_iterator itEnd) {
    Ray r(position_ + 1e-12 * direction_.vector(), direction_);

    bool found = false;
    double nearestDistance = 0.0;
    Point* p = NULL;
    nextIntersection_.pSurface = NULL;
    Ray_intersection intersection;

    while (itSurface != itEnd) {
#ifdef USE_ALL_INTERSECTIONS
        std::vector<Ray_intersection> intersections;
        (*itSurface)->computeAllIntersections(r,
            std::back_inserter(intersections));
        for (Ray_intersection intersection : intersections) {
#else
        (*itSurface)->computeFirstIntersection(r, intersection);
#endif
            if (intersection &&
                (p = boost::get<Point>(&(intersection->first)))) {
                double distance = CGAL::squared_distance(position_, *p);
                if (!found || distance < nearestDistance) {
                    found = true;
                    nextIntersection_.point = *p;
                    nextIntersection_.faceId = intersection->second;
                    nextIntersection_.pSurface = *itSurface;
                    nearestDistance = distance;
                }
            }
#ifdef USE_ALL_INTERSECTIONS
        }
#endif

        itSurface++;
    }

    return found;
}
