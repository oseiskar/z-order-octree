#include "octree.hpp"

#include <array>
#include <vector>
#include <random>

int main() {
    using Coordinate = double;
    using Vertex = std::array<Coordinate, 3>;

    std::default_random_engine generator;
    std::normal_distribution<Coordinate> gauss;
    std::vector<Vertex> vertices;

    constexpr size_t N_WALKS = 10;
    constexpr size_t N_VERT_PER_WALK = 10;
    constexpr Coordinate WALK_START_STDEV = 5.0;
    constexpr Coordinate WALK_STEP_STDEV = 0.1;

    for (size_t iWalk = 0; iWalk < N_WALKS; ++iWalk) {
        Vertex v;
        for (int i = 0; i < 3; ++i) v[i] = gauss(generator) * WALK_START_STDEV;
        for (size_t iStep = 0; iStep < N_VERT_PER_WALK; ++iStep) {
            vertices.push_back(v);
            // log_debug("%g, %g, %g", v[0], v[1], v[2]);
            for (int i = 0; i < 3; ++i) v[i] += gauss(generator) * WALK_START_STDEV;
        }
    }

    ZOrderOctree<Vertex, Coordinate> octree(
        vertices.data(),
        vertices.data() + vertices.size(),
        { .leafSize = 0.1 });

    for (const auto *point : octree.lookup(Vertex { 0, 0, 0 }, 9)) {
        log_debug("%g\t%g\t%g", (*point)[0], (*point)[1], (*point)[2]);
    }

    return 0;
}
