#include "surfaceemission.h"
#include "argon_data.h"

SurfaceEmission::SurfaceEmission(Surface* pSurface, double emissionRate)
    : pSurface_(pSurface), emissionRate_(emissionRate) {}

void SurfaceEmission::generateParticle(Particle& particle, Rng& rng) const {
    Point p;
    face_descriptor fd;

    std::tie(p, fd) = pSurface_->getRandomPoint(rng);
    Direction d = pSurface_->generateCosineLawDirection(fd, rng);
    particle.setMass_eV(ARGON_DATA.mass * ATOMIC_MASS_TO_EV);
    double v = Util::getMBSpeed(rng, pSurface_->getTemperature(),
        particle.getMass_eV());
    particle.setPosition(p);
    particle.setVelocity(v, d);
}

