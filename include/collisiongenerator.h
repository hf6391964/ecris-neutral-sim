#pragma once

#include <vector>

#include "collisionreaction.h"

class CollisionGenerator {
    std::vector<CollisionReaction *> collisionReactions_;

    public:
        void addCollisionReaction(CollisionReaction *);
};

