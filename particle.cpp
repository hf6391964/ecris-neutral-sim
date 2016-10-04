#include "particle.h"

void Particle::goToIntersection(Rng& rng) {
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

bool Particle::findNextIntersection(
    std::vector<Surface*>::const_iterator itSurface,
    std::vector<Surface*>::const_iterator itEnd) {
    Direction direction = getDirection();
    Point position = getPosition() + 1e-12 * direction.vector();
    Ray r(position, direction);

    bool found = false;
    double nearestDistance = 0.0;
    Point* p = NULL;
    nextIntersection_.pSurface = NULL;
    Ray_intersection intersection;

    while (itSurface != itEnd) {
        (*itSurface)->computeFirstIntersection(r, intersection);

        if (intersection &&
            (p = boost::get<Point>(&(intersection->first)))) {
            double distance = CGAL::squared_distance(position_, *p);
            if (!found || distance < nearestDistance) {
                found = true;
                nextIntersection_.point = *p;
                nextIntersection_.faceId = intersection->second;
                nextIntersection_.pSurface = *itSurface;
                nearestDistance = distance;
#ifdef DEBUG
                std::cout << "Next intersection: ";
                Util::printPoint(*p);
                std::cout << std::endl;
#endif
            }
        }

        itSurface++;
    }

    return found;
}

