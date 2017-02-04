#pragma once

#include "cgal_and_typedefs.h"
#include "spatialdistribution.h"
#include "constants.h"
#include "element_data.h"

class ParticlePopulation {
    protected:
        Element element_ = ELEMENT_NONE;
        const double particleMass_eV_;
        int chargeState_;
        DensityDistribution densityDistribution_;
        std::string label_ = "";

    public:
        ParticlePopulation(Element element, int chargeState,
            DensityDistribution densityDistribution)
            : element_(element),
              particleMass_eV_(ELEMENT_DATA.at(element_)->mass * ATOMIC_MASS_TO_EV),
              chargeState_(chargeState), densityDistribution_(densityDistribution) {
            if (chargeState_ == -1) {
                label_ = "electrons";
            } else {
                label_ = std::to_string(chargeState_) + "+ ions";
            }
        }

        ParticlePopulation(double particleMass_eV, int chargeState,
            DensityDistribution densityDistribution)
            : particleMass_eV_(particleMass_eV), chargeState_(chargeState),
              densityDistribution_(densityDistribution) {
            if (chargeState_ == -1) {
                label_ = "electrons";
            } else {
                label_ = std::to_string(chargeState_) + "+ ions";
            }
        }

        virtual ~ParticlePopulation() {};

        Element getElement() const {
            return element_;
        }

        std::string getLabel() const {
            return label_;
        }

        virtual Vector getRandomParticleVelocity(Rng &rng) const = 0;
        double getRandomParticleSpeed(Rng &rng) const;
        Point getRandomParticlePosition(Rng &rng) const;

        double getParticleMass_eV() const { return particleMass_eV_; }
        int getChargeState() const { return chargeState_; }

        double getDensityAt(const Point &p) const;

        virtual double getRelativeSpeed(double particleSpeed, monte_state *ms,
            gsl_rng *rng) const = 0;

        virtual double calculateRateCoefficient(double particleSpeed,
            monte_state *ms, gsl_rng *rng,
            double (*fCrossSection)(double, void *), void *fArgs) const = 0;

        void setCoordinateTransformation(const Aff_transformation
            &transformation);
        void removeCoordinateTransformation();
};

