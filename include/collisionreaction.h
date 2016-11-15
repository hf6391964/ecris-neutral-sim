#pragma once

#include "cgal_and_typedefs.h"
#include "plasmadensities.h"

// Base class for a general collision reaction
class CollisionReaction {
    public:
        virtual double getMeanReactionRate(const Point &p) const = 0;
        virtual double getMajorantReactionRate(const Point &p) const = 0;

        virtual double getCrossSection(const Point &p, double speed)
            const = 0;
};

