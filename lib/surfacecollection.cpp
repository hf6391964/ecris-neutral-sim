#include "surfacecollection.h"

#define USE_ALL_INTERSECTIONS

void SurfaceCollection::addSurface(Surface *surface) {
    surfaces_.push_back(surface);
    bbox_ += surface->bbox();
}

Bbox SurfaceCollection::bbox() const {
    return bbox_;
}

bool SurfaceCollection::findClosestIntersection(Ray &r,
    IntersectionPoint &ip) const {
    bool found = false;
    double nearestDistance = 0.0;
    Point* p = NULL;
    Ray_intersection intersection;
    Point position = r.point(0);

    for (Surface *pSurface : surfaces_) {
#ifdef USE_ALL_INTERSECTIONS
        std::vector<Ray_intersection> intersections;
        pSurface->computeAllIntersections(r,
            std::back_inserter(intersections));
        for (Ray_intersection intersection : intersections) {
#else
        pSurface->computeFirstIntersection(r, intersection);
#endif
            if (intersection &&
                (p = boost::get<Point>(&(intersection->first)))) {
                double distance = CGAL::squared_distance(position, *p);
                if (!found || distance < nearestDistance) {
                    found = true;
                    ip.point = *p;
                    ip.faceId = intersection->second;
                    ip.pSurface = pSurface;
                    nearestDistance = distance;
                }
            }
#ifdef USE_ALL_INTERSECTIONS
        }
#endif
    }

    return found;
}

