#include "surfacecollection.h"

//#define USE_ALL_INTERSECTIONS

void SurfaceCollection::addSurface(std::shared_ptr<Surface> surface) {
    surfaces_.push_back(surface);
    bbox_ += surface->bbox();
}

Bbox SurfaceCollection::bbox() const {
    return bbox_;
}

bool SurfaceCollection::findClosestIntersection(Ray &r,
    IntersectionPoint &ip, bool skipCurrentFace) const {
    bool found = false;
    double nearestDistance = 0.0;
    Point* p = NULL;
    Ray_intersection intersection;
    Point position = r.point(0);
    skipCurrentFace = skipCurrentFace && ip.pSurface != NULL;
    ip.pSurface = NULL;
    Skip skip(ip.faceId);

    for (auto pSurface : surfaces_) {
#ifdef USE_ALL_INTERSECTIONS
        std::vector<Ray_intersection> intersections;
        pSurface->computeAllIntersections(r,
            std::back_inserter(intersections));
        for (Ray_intersection intersection : intersections) {
#else
        if (skipCurrentFace) {
            pSurface->computeFirstIntersection(r, intersection, true, skip);
        } else {
            pSurface->computeFirstIntersection(r, intersection, false, skip);
        }
#endif
            if (intersection &&
                (p = boost::get<Point>(&(intersection->first)))) {
                double distance = CGAL::squared_distance(position, *p);
                if (!found || distance < nearestDistance) {
                    found = true;
                    ip.point = *p;
                    ip.faceId = intersection->second;
                    ip.pSurface = pSurface.get();
                    nearestDistance = distance;
                }
            }
#ifdef USE_ALL_INTERSECTIONS
        }
#endif
    }

    return found;
}

void SurfaceCollection::writeStatistics(Logger &log) const {
    log << "Surface pumping and collision statistics:\n";
    for (auto pSurface : surfaces_) {
        log << pSurface->getLabel() << ": " << pSurface->getCollisionCount() <<
            " collisions, " << pSurface->getPumpedParticles() <<
            " particles pumped\n";
    }
    log << '\n';
}

