#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <limits>

#include "cgal_and_typedefs.h"
#include "maxwellianpopulation.h"
#include "collisiongenerator.h"
#include "element_data.h"
#include "neutralizationgenerator.h"

class SimplePlasmaModel {
    std::vector<std::shared_ptr<MaxwellianPopulation>> particlePopulations_;
    unsigned int maxChargeState_ = 0;
    const ElementData *elementData_;
    const Element element_;
    double ionMassEv_ = 0.0;
    double maxIonTemp_ = 0.0;

    public:
        // Electron weight is the statistical weight of an electron in the
        // electron density distribution i.e.
        // # of real particles / test particle.
        //
        // Ion densities are relative to real electron density,
        // ionRelativeDensities[i] corresponds to charge state i+1
        // Quasineutrality should be satisfied:
        // sum(ionRelativeDensities[i] * (i+1)) = 1
        //
        // electron temperature is eV
        //
        // ion temperatures (eV) at indices corresponding to charge state
        //
        // ion mass should be in eV
        SimplePlasmaModel(std::string electronDensityFilename,
            double electronWeight, std::vector<double> ionRelativeDensities,
            double electronTemperature,
            std::vector<double> ionTemperatures, Element element);

        ~SimplePlasmaModel();

        void populateCollisionReactions(CollisionGenerator &generator,
            simthreadresources *thread_res, double speedStepSize = 0.01) const;
        void populateNeutralizationReactions(
            NeutralizationGenerator &generator,
            double confinementTime,
            const std::string &wallFilename,
            const std::string &end1Filename,
            const std::string &end2Filename,
            const std::string &flychkFilename,
            const SurfaceCollection &surfaces,
            double averageElectronDensity) const;

        void setCoordinateTransformation(const Aff_transformation &tf);

        // Here ion densities are indexed such that i = q.
        // A special case is i = 0 which corresponds to electron density.
        double getIonDensityAt(const Point &p, unsigned int chargeState) const;
        double getElectronDensityAt(const Point &p) const;

        double getIonMass() const {
            return ionMassEv_;
        }

        Vector getIsotropicIonVelocity(Rng &rng, unsigned int chargeState)
            const;
        Vector getIsotropicElectronVelocity(Rng &rng) const;
};

