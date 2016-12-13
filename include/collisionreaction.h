#pragma once

#include "cgal_and_typedefs.h"
#include "particlepopulation.h"
#include "particle.h"

// Base class for a general collision reaction
class CollisionReaction {
    protected:
        std::shared_ptr<ParticlePopulation> population_;

    public:
        CollisionReaction(std::shared_ptr<ParticlePopulation> population)
            : population_(population) {}

        std::shared_ptr<ParticlePopulation> getPopulation() const {
            return population_;
        }

        virtual double getReactionRate(const Point &p, double particleSpeed,
            simthreadresources *thread_res) const = 0;

        // TODO count energy and particles "lost" in the plasma
        virtual std::vector<Particle> computeReactionProducts(Rng &rng,
            const Point &p, const Particle &target) const = 0;

        virtual double getCrossSection(double relativeSpeed) const = 0;
};

