#ifndef SURFACE_H
#define SURFACE_H

#include <string>

#include "cgal_and_typedefs.h"

class Surface {
    Polyhedron mesh_;

    public:
        static Surface loadFromSTL(std::string filename,
            double avgTriangleArea = -1.0);

        bool isLoaded();
};

#endif
