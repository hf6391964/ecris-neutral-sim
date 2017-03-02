#pragma once

#include <atomic>

#include "cgal_and_typedefs.h"
#include "particlepopulation.h"
#include "particle.h"

typedef std::pair<std::vector<Particle>, unsigned int> CollisionProducts;

struct mc_integrate_resources {
    monte_state *ms;
    gsl_integration_workspace *ws;
    gsl_rng *gslrng;

    mc_integrate_resources(uint_least32_t seed) {
        ms = gsl_monte_vegas_alloc(3);
        ws = gsl_integration_workspace_alloc(RATE_COEFF_WORKSPACE_SIZE);
        gslrng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(gslrng, seed);
    }

    ~mc_integrate_resources() {
        gsl_monte_vegas_free(ms);
        gsl_integration_workspace_free(ws);
        gsl_rng_free(gslrng);
    }
};

// Base class for a general collision reaction
class CollisionReaction {
    private:
        std::atomic_ulong reactionCounter_;

    protected:
        std::shared_ptr<ParticlePopulation> population_;
        std::string label_ = "";

    public:
        CollisionReaction(std::shared_ptr<ParticlePopulation> population);

        std::shared_ptr<ParticlePopulation> getPopulation() const;

        std::string getLabel() const;

        virtual double getReactionRate(const Point &p, double particleSpeed,
            mc_integrate_resources &mc_res) const = 0;

        virtual double getRateCoefficient(double particleSpeed,
            mc_integrate_resources &mc_res) const = 0;

        // computeReactionProducts returns a vector of neutral reaction
        // products and the number of ionized reaction products
        virtual CollisionProducts computeReactionProducts(Rng &rng,
            const Point &p, const Particle &target) const = 0;

        virtual double getCrossSection(double relativeSpeed) const = 0;

        void incrementReactionCounter();

        unsigned long getReactionCount() const;
};
