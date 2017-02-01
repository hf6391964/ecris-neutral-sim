#pragma once

#include <vector>

#include "cgal_and_typedefs.h"
#include "surface.h"
#include "logger.h"

class SurfaceCollection {
    std::vector<Surface*> surfaces_;
    Bbox bbox_;

    public:
        SurfaceCollection() {};

        Bbox bbox() const;

        void addSurface(Surface *surface);

        bool findClosestIntersection(Ray &r, IntersectionPoint &p) const;

        void writeStatistics(Logger& log) const;
};

