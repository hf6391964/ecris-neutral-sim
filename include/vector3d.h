#pragma once

#include <vector>

#include "grid.h"

template<typename T>
class Vector3D {
    public:
        Vector3D(const Grid &grid, T value) {
            unsigned int ix, iy, iz;
            std::tie(ix, iy, iz) = grid.dimensions();
            vec_ = std::vector<std::vector<std::vector<T>>>(iz,
                std::vector<std::vector<T>>(iy, std::vector<T>(ix, value)));
        }

        T &at(unsigned int ix, unsigned int iy, unsigned int iz) {
            return vec_.at(iz).at(iy).at(ix);
        }

        T &at(Index3D &idx) {
            unsigned int ix, iy, iz;
            std::tie(ix, iy, iz) = idx;
            return at(ix, iy, iz);
        }

        const T &at(unsigned int ix, unsigned int iy, unsigned int iz) const {
            return vec_.at(iz).at(iy).at(ix);
        }

        const T &at(Index3D &idx) const {
            unsigned int ix, iy, iz;
            std::tie(ix, iy, iz) = idx;
            return at(ix, iy, iz);
        }

    private:
        std::vector<std::vector<std::vector<T>>> vec_;
};

