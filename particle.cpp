#include "particle.h"

#define USE_ALL_INTERSECTIONS

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

