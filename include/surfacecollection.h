#pragma once

#include <vector>

#include "cgal_and_typedefs.h"
#include "surface.h"
#include "logger.h"

class SurfaceCollection {
    std::vector<std::shared_ptr<Surface>> surfaces_;
    Bbox bbox_;

    public:
        SurfaceCollection() {};

        Bbox bbox() const;

        void addSurface(const std::shared_ptr<Surface> &surface);

        bool findClosestIntersection(Ray &r, IntersectionPoint &p,
            bool skipCurrentFace = false) const;

        void writeStatistics(Logger& log) const;
};

