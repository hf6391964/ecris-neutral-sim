#include "particle.h"

#define USE_ALL_INTERSECTIONS

Particle::Particle(Element element, double time)
    : element_(element), elementData_(ELEMENT_DATA.at(element_)),
      mass_eV_(elementData_->mass * ATOMIC_MASS_TO_EV),
      time_(time) {}

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

double Particle::getMass_eV() const {
    return mass_eV_;
}

Particle::State Particle::getState() const {
    return state_;
}

double Particle::getSpeed() const {
    return speed_;
}

double Particle::getTime() const {
    return time_;
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

bool Particle::hasNextIntersection() const {
    return nextIntersection_.pSurface != NULL;
}

double Particle::distanceToIntersection() const {
    return std::sqrt(CGAL::squared_distance(position_,
        nextIntersection_.point));
}

void Particle::goForward(double dt) {
    time_ += dt;
    double stepLength = dt * speed_;
    position_ = position_ + stepLength * direction_.vector();
}

void Particle::goToIntersection(Rng& rng) {
    if (!hasNextIntersection()) return;

    double distance =
        std::sqrt(CGAL::squared_distance(position_, nextIntersection_.point));
    double dt = distance / speed_;
    time_ += dt;

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
    nextIntersection_.pSurface->incrementCollisionCounter();

    if (nextIntersection_.pSurface->checkIfPumped(rng)) {
        nextIntersection_.pSurface->addPumpedParticle();
        state_ = Pumped;
    }
}

bool Particle::findNextIntersection(SurfaceCollection &surfaces) {
    Ray r(position_ + 1e-12 * direction_.vector(), direction_);

    return surfaces.findClosestIntersection(r, nextIntersection_);
}

void Particle::setNextIntersection(const IntersectionPoint &ip) {
    nextIntersection_ = ip;
}

IntersectionPoint Particle::getNextIntersection() const {
    return nextIntersection_;
}

void Particle::setState(Particle::State state) {
    state_ = state;
}
