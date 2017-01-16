#include "electronmodel.h"
#include "util.h"

ElectronModel::ElectronModel(const double B0, const double r0, const double dt,
    const double z1, const double z2, const double gridSize,
    const double confinementTime, const double a[]) : B0_(B0), r0_(r0), dt_(dt),
    z1_(z1), z2_(z2), gridSize_(gridSize), confinementTime_(confinementTime),
    a_(a) {
    grid_ = Grid(Bbox(-r0, -r0, z1, r0, r0, z2), gridSize_);
    resetCounters();
}

ElectronModel::~ElectronModel() {
    if (particleCount_ != NULL) {
        delete[] particleCount_;
    }
}

void ElectronModel::resetCounters() {
    if (particleCount_ != NULL) {
        delete[] particleCount_;
    }

    size_t arraySize = grid_.arraySize();
    particleCount_ = new unsigned long[arraySize];
    for (size_t i = 0; i < arraySize; i++) {
        particleCount_[i] = 0;
    }

    z1CollisionPoints_.clear();
    z2CollisionPoints_.clear();
    cylinderCollisionPoints_.clear();
    finalEnergies_.clear();
    lostVelocities_.clear();
    nonLostVelocities_.clear();
}

Vector ElectronModel::solenoidBfield(const double x, const double y,
    const double z, const double a[]) {
    double r2 = x*x + y*y;
    // Horner scheme w.r.t. z
    double solenoid_Bz =
           a[0] + r2*(0.5*a[2] + r2*((3.0/8.0)*a[4] - r2*(5.0/16.0)*a[6])) +
        z*(a[1] + r2*(-1.5*a[3] + r2*(15.0/8.0)*a[5]) +
        z*(a[2] + r2*(-3.0*a[4] + r2*(45.0/8.0)*a[6]) +
        z*(a[3] - r2*5.0*a[5] +
        z*(a[4] - r2*(15.0/2.0)*a[6] +
        z*(a[5] + z*a[6])))));
    double solenoid_Br_divided_by_r =
           -0.5*a[1] + r2*((3.0/8.0) - r2*(5.0/16.0)*a[5]) +
        z*(-a[2] + r2*(1.5*a[4] - r2*(15.0/8.0)*a[6]) +
        z*(-1.5*a[3] + r2*(15.0/4.0)*a[5] +
        z*(-2.0*a[4] + r2*(15.0/2.0)*a[6] +
        z*(-2.5*a[5] + z*(-3.0*a[6])))));

    return Vector(
        solenoid_Br_divided_by_r * x,
        solenoid_Br_divided_by_r * y,
        solenoid_Bz
    );
}

Vector ElectronModel::hexapoleBfield(const double x, const double y,
    const double B0, const double r0) {
    double C = -2.0*B0 / (r0*r0);
    return Vector(C*x*y, C*(x*x - y*y) / 2.0, 0.0);
}

Vector ElectronModel::totalBfield(const Vector vx, const double B0,
    const double r0, const double a[]) {
    double x = vx.x();
    double y = vx.y();
    double z = vx.z();
    return solenoidBfield(x, y, z, a) + hexapoleBfield(x, y, B0, r0);
}

bool ElectronModel::moveParticle() {
    Vector B = totalBfield(particlePosition_, B0_, r0_, a_);

    // Boris integrator with zero E field

    double Bsqnorm = B.squared_length();
    double Bnorm = std::sqrt(Bsqnorm);
    double f1 = std::tan(-ELECTRON_CHARGE_MASS_RATIO * dt_ * Bnorm / 2.0) /
        Bnorm;
    double f2 = 2.0 * f1 / (1.0 + f1*f1 * Bsqnorm);

    Vector v = particleVelocity_ + f1*CGAL::cross_product(particleVelocity_, B);
    particleVelocity_ = particleVelocity_ + f2*CGAL::cross_product(v, B);
    particlePosition_ = particlePosition_ + dt_ * particleVelocity_;

    t_ += dt_;

    double x = particlePosition_.x(),
           y = particlePosition_.y(),
           z = particlePosition_.z();
    if (z < z1_) {
        z1CollisionPoints_.push_back(particlePosition_);
    } else if (z > z2_) {
        z2CollisionPoints_.push_back(particlePosition_);
    } else if (x*x + y*y > r0_*r0_) {
        cylinderCollisionPoints_.push_back(particlePosition_);
    } else {
        size_t i;
        if (grid_.arrayIndex(x, y, z, i)) {
            particleCount_[i] += 1;
        }
        return false;
    }

    return true;
}

void ElectronModel::runSimulation(const unsigned long nParticles,
    std::function<double(Rng&)> getSpeed, const double Becr, Rng& rng) {
    for (unsigned long n = 0; n < nParticles; n++) {
        newParticle(getSpeed(rng),
            [Becr](Vector B) {
                return B.squared_length() < Becr*Becr;
            }, rng);

        Vector v = particleVelocity_;
        Vector B = totalBfield(particlePosition_, B0_, r0_, a_);

        double vRadial = v * B / std::sqrt(B.squared_length());
        Vector vt = v - vRadial * B / std::sqrt(B.squared_length());
        double vTransverse = std::sqrt(vt.squared_length());

        while (t_ < confinementTime_) {
            if (moveParticle()) break;
        }

        if (t_ < confinementTime_) {
            lostVelocities_.push_back(
                std::tuple<double, double>(vRadial, vTransverse));
        } else {
            nonLostVelocities_.push_back(
                std::tuple<double, double>(vRadial, vTransverse));
        }

        finalEnergies_.push_back(energy());
    }
}

void ElectronModel::runMonoenergeticSimulation(const unsigned long nParticles,
    const double energy, const double Becr, Rng& rng) {
    runSimulation(nParticles, [energy](Rng &) {
        return std::sqrt(2.0 * energy / ELECTRON_MASS_EV) * SPEED_OF_LIGHT;
    }, Becr, rng);
}

void ElectronModel::runMaxwellSimulation(const unsigned long nParticles,
    const double T_eV, const double Becr, Rng& rng) {
    runSimulation(nParticles, [T_eV](Rng &rng_) {
        return Util::getMBSpeed(rng_, T_eV, ELECTRON_MASS_EV);
    }, Becr, rng);
}

double ElectronModel::meanEnergy() const {
    return std::accumulate(finalEnergies_.begin(), finalEnergies_.end(),
        (double)0.0) / finalEnergies_.size();
}

double ElectronModel::energyStdDev() const {
    double mean = meanEnergy();
    return std::sqrt(std::accumulate(finalEnergies_.begin(),
        finalEnergies_.end(), (double)0.0, [mean](double acc, double e) {
            double de = e - mean;
            return acc + de*de;
        }
    ) / finalEnergies_.size());
}

void ElectronModel::writeDensityToFile(const std::string &filename,
    bool normalize) const {
    std::ofstream fout(filename);
    grid_.writeDimensions(fout);
    size_t N = grid_.arraySize();
    if (normalize) {
        unsigned long maxCount = 0;
        for (size_t i = 0; i < N; ++i) {
            maxCount = std::max(maxCount, particleCount_[i]);
        }

        unsigned long long nParticles = 0;
        unsigned long nPlasmaCells = 0;
        for (size_t i = 0; i < N; ++i) {
            if (particleCount_[i] * 4 >= maxCount) {
                nParticles += particleCount_[i];
                nPlasmaCells += 1;
            }
        }

        double averageDensity = (double)nParticles / nPlasmaCells,
               normalizationFactor = 1.0 / averageDensity;
        for (size_t i = 0; i < N; ++i) {
            double value = particleCount_[i] * normalizationFactor;
            fout << value << std::endl;
        }
    } else {
        for (size_t i = 0; i < N; ++i) {
            fout << particleCount_[i] << std::endl;
        }
    }
    fout.close();
}

void ElectronModel::writeElectronEndpoints(const std::string &z1Filename,
    const std::string &z2Filename, const std::string &radialFilename) const {
    std::ofstream fout;
    fout.open(z1Filename);
    for (Vector pos : z1CollisionPoints()) {
        fout << pos.x() << CSV_SEP << pos.y() << CSV_SEP << pos.z() << std::endl;
    }
    fout.close();
    fout.open(z2Filename);
    for (Vector pos : z2CollisionPoints()) {
        fout << pos.x() << CSV_SEP << pos.y() << CSV_SEP << pos.z() << std::endl;
    }
    fout.close();
    fout.open(radialFilename);
    for (Vector pos : cylinderCollisionPoints()) {
        fout << pos.x() << CSV_SEP << pos.y() << CSV_SEP << pos.z() << std::endl;
    }
    fout.close();
}

