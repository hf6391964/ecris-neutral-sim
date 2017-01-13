#include "surfaceemission.h"

SurfaceEmission::SurfaceEmission(Surface* pSurface, double emissionRate,
    Element element)
    : pSurface_(pSurface), emissionRate_(emissionRate), element_(element) {
    elementData_ = ELEMENT_DATA.at(element_);
}

void SurfaceEmission::generateParticle(Particle& particle, Rng& rng) const {
    Point p;
    face_descriptor fd;

    std::tie(p, fd) = pSurface_->getRandomPoint(rng);
    Direction d = pSurface_->generateCosineLawDirection(fd, rng);
    particle.setMass_eV(elementData_->mass * ATOMIC_MASS_TO_EV);
    double v = Util::getMBSpeed(rng, pSurface_->getTemperature(),
        particle.getMass_eV());
    particle.setPosition(p);
    particle.setVelocity(v, d);
}

