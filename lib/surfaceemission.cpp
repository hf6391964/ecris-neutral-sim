#include "surfaceemission.h"

SurfaceEmission::SurfaceEmission(Surface* pSurface, double emissionRate)
    : pSurface_(pSurface), emissionRate_(emissionRate) {}

void SurfaceEmission::generateParticle(Particle& particle, Rng& rng) const {
    Point p;
    face_descriptor fd;

    std::tie(p, fd) = pSurface_->getRandomPoint(rng);
    Direction d = pSurface_->generateCosineLawDirection(fd, rng);
    double v = Util::getMBSpeed(rng, particle.getTemperature(),
        particle.getMass_eV());
    particle.setPosition(p);
    particle.setVelocity(v, d);
}

