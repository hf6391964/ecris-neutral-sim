#include "surfaceemission.h"

SurfaceEmission::SurfaceEmission(Surface* pSurface, double emissionRate,
    Element element)
    : pSurface_(pSurface), emissionRate_(emissionRate), element_(element) {
    elementData_ = ELEMENT_DATA.at(element_);
}

Particle SurfaceEmission::generateParticle(Rng& rng, double time) const {
    Point p;
    face_descriptor fd;
    Particle particle(element_, time);

    std::tie(p, fd) = pSurface_->getRandomPoint(rng);
    Direction d = pSurface_->generateCosineLawDirection(fd, rng);
    double v = Util::getMBSpeed(rng, pSurface_->getTemperature(),
        particle.getMass_eV());
    particle.setPosition(p);
    particle.setVelocity(v, d);
    return particle;
}

