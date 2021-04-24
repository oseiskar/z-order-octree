#include "octree.hpp"

#include <array>
#include <vector>
#include <random>

#include <cstdio>

namespace {
template <class Node> void traverse(const Node &node) {
    if (node.elements().size() > 5) {
      printf("node with %zu elements @ level %d\n",
        node.elements().size(),
        node.getLevel());
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

template <class Vertex> void generateData(std::vector<Vertex> &vertices) {
    std::default_random_engine generator;
    std::normal_distribution<double> gauss;

    constexpr size_t N_WALKS = 10;
    constexpr size_t N_VERT_PER_WALK = 10;
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
}

int main() {
    using Coordinate = double;
    using Vertex = std::array<Coordinate, 3>;

    std::vector<Vertex> vertices;
    generateData(vertices);

    ZOrderOctree<Vertex, Coordinate> octree(
        vertices.data(),
        vertices.data() + vertices.size(),
        { .leafSize = 0.1, .rootLevel = 10 });

    for (const auto *point : octree.lookup(Vertex { 0, 0, 0 }, 9).elements()) {
        printf("%g\t%g\t%g\n", (*point)[0], (*point)[1], (*point)[2]);
    }

    printf("%zu elements in a probably empty node\n",
        octree.lookup(Vertex { 10, 10, 10 }, 0).elements().size());

    traverse(octree.root());

    for (const auto node : octree.nodesAtLevel(8)) {
        printf("lev %d: non-empty node with %zu elements\n", node.getLevel(), node.elements().size());
    }

    return 0;
}
