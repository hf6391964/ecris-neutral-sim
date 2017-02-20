#include "surfaceemission.h"

SurfaceEmission::SurfaceEmission(std::shared_ptr<Surface> pSurface,
    double emissionRate, Element element, std::string label)
    : NeutralSource(label),
      pSurface_(pSurface), emissionRate_(emissionRate), element_(element) {
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
    IntersectionPoint ip;
    ip.pSurface = pSurface_.get();
    ip.faceId = fd;
    ip.point = p;
    particle.setNextIntersection(ip);
    particle.setPosition(p);
    particle.setVelocity(v, d);
    return particle;
}

