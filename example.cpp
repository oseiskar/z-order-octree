#include "octree.hpp"

#include <array>
#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

namespace {
struct Timer {
    std::chrono::steady_clock::time_point t0;
    const char *name;
    Timer(const char *name) : name(name) {
        t0 = std::chrono::steady_clock::now();
    }
    void stop() {
        if (!name) return;
        const double tMs = (std::chrono::steady_clock::now() - t0).count() / 1e6;
        printf("TIMER %s: %dms\n", name, int(tMs + 0.5));
        name = nullptr;
    }
    ~Timer() {
        stop();
    }
};

template <class Vertex> void generateData(std::vector<Vertex> &vertices) {
    Timer timer("data generation");
    std::default_random_engine generator;
    std::normal_distribution<double> gauss;

    constexpr size_t N_WALKS = 1000;
    constexpr size_t N_VERT_PER_WALK = 1000;
    constexpr double WALK_START_STDEV = 5.0;
    constexpr double WALK_STEP_STDEV = 0.1;

    for (size_t iWalk = 0; iWalk < N_WALKS; ++iWalk) {
        Vertex v;
        for (int i = 0; i < 3; ++i) v[i] = gauss(generator) * WALK_START_STDEV;
        for (size_t iStep = 0; iStep < N_VERT_PER_WALK; ++iStep) {
            vertices.push_back(v);
            // log_debug("%g, %g, %g", v[0], v[1], v[2]);
            for (int i = 0; i < 3; ++i) v[i] += gauss(generator) * WALK_START_STDEV;
        }
    }
}

template <class Node> void traverse(const Node &node) {
    if (node.elements().size() > 130000) {
      const auto c = node.center();
      printf("node (%g, %g, %g) with %zu elements @ level %d, size %g\n",
        c[0], c[1], c[2],
        node.elements().size(),
        node.getLevel(),
        node.sideLength());
    }

    /*log_debug("%s%s%s",
      node.isLastSibling() ? " last sibling" : "",
      node.isLastAtThisLevel() ? " last at level" : "",
      node.isLeaf() ? " leaf" : "");*/

    if (!node.isLeaf()) for (auto child : node.children()) {
        // log_debug("child of node at level %d", node.getLevel());
        traverse(child);
    }
}
}

int main() {
    using Coordinate = double;
    using Vertex = std::array<Coordinate, 3>;

    std::vector<Vertex> vertices;
    generateData(vertices);
    printf("generated %zu points\n", vertices.size());

    Timer timer("octree build");
    ZOrderOctree<Vertex, Coordinate> octree(
        vertices.data(),
        vertices.data() + vertices.size(),
        { .leafSize = 0.1 });
    timer.stop();

    {
        Timer timer("lookup");
        for (const auto *point : octree.lookup(Vertex { 0, 0, 0 }, 4).elements()) {
            printf("%g\t%g\t%g\n", (*point)[0], (*point)[1], (*point)[2]);
        }
    }

    printf("%zu elements in a probably empty node\n",
        octree.lookup(Vertex { 10, 10, 10 }, 0).elements().size());

    {
        Timer timer("full traversal");
        traverse(octree.root());
    }

    {
        Timer timer("level traversal");
        for (const auto node : octree.nodesAtLevel(16)) {
            printf("lev %d: non-empty node with %zu elements\n", node.getLevel(), node.elements().size());
        }
    }

    return 0;
}
