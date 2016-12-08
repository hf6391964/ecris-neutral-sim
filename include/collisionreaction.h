#pragma once

#include "cgal_and_typedefs.h"
#include "plasmadensities.h"
#include "particle.h"

// Base class for a general collision reaction
class CollisionReaction {
    public:
        virtual double getMeanReactionRate(const Point &p,
            double relativeSpeed) const = 0;
        virtual double getMajorantReactionRate(const Point &p,
            double relativeSpeed) const = 0;

        virtual std::vector<Particle> computeReactionProducts(
            const Point &p, const Particle &target) const = 0;

        virtual double getCrossSection(double relativeSpeed) const = 0;
};

