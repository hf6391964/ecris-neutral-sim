#pragma once

#include "cgal_and_typedefs.h"
#include "particlepopulation.h"
#include "particle.h"

typedef std::pair<std::vector<Particle>, unsigned int> CollisionProducts;

// Base class for a general collision reaction
class CollisionReaction {
    protected:
        std::shared_ptr<ParticlePopulation> population_;

    public:
        CollisionReaction(std::shared_ptr<ParticlePopulation> population);

        std::shared_ptr<ParticlePopulation> getPopulation() const;

        virtual double getReactionRate(const Point &p, double particleSpeed,
            simthreadresources *thread_res) const = 0;

        virtual double getRateCoefficient(double particleSpeed,
            simthreadresources *thread_res) const = 0;

        // computeReactionProducts returns a vector of neutral reaction
        // products and the number of ionized reaction products
        virtual CollisionProducts computeReactionProducts(Rng &rng,
            const Point &p, const Particle &target) const = 0;

        virtual double getCrossSection(double relativeSpeed) const = 0;
};
