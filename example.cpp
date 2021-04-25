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

template <class Node> size_t traverse(const Node &node) {
    const size_t nEls = node.elements().size();
    if (nEls > 130000) {
      const auto c = node.center();
      printf("node (%g, %g, %g) with %zu elements @ level %d, size %g\n",
        c[0], c[1], c[2],
        nEls,
        node.getLevel(),
        node.sideLength());

      auto corner0 = node.minCorner();
      auto corner1 = node.maxCorner();

      for (const auto *point : node.elements()) {
          for (int c = 0; c < 3; ++c) {
              float coord = (*point)[c];
              constexpr double MARGIN = 0.01;
              // printf("c%d: %g, [%g, %g]\n", c, coord, corner0[c], corner1[c]);
              assert(coord > corner0[c] - MARGIN && coord < corner1[c] + MARGIN);
          }
      }
    }

    if (!node.isLeaf()) {
      size_t childSizes = 0;
      for (auto child : node.children()) {
        // log_debug("child of node at level %d", node.getLevel());
        childSizes += traverse(child);
      }
      assert(childSizes == nEls);
    }

    return nEls;
}
}

int main() {
    using Coordinate = double;
    using Vertex = std::array<Coordinate, 3>;

    std::vector<Vertex> vertices;
    generateData(vertices);
    printf("generated %zu points\n", vertices.size());

    ZOrderOctree<Vertex, Coordinate> octree({ .leafSize = 0.1 });

    {
      Timer timer("octree build 1/2");
      octree.addData(
          vertices.data(),
          vertices.size() / 2);
    }
    {
      Timer timer("octree build 1/2");
      octree.addData(
          vertices.data() + vertices.size() / 2,
          vertices.size() - vertices.size() / 2);
    }
    {
      Timer timer("octree data removal");
      size_t nRemoved = 0;
      octree.removeData([&nRemoved](const Vertex &v) -> bool {
          if (v[0] >= 10.0 && v[1] > 20) {
            nRemoved++;
            return true;
          }
          return false;
      });
      printf("removed %zu point(s)\n", nRemoved);
    }

    octree.clearWorkspace();

    {
        Timer timer("lookup");
        auto node = octree.lookup(Vertex { 0, 0, 0 }, 4);
        for (const auto *point : node.elements()) {
            printf("%g\t%g\t%g\n", (*point)[0], (*point)[1], (*point)[2]);
        }
    }

    printf("%zu elements in a probably empty node\n",
        octree.lookup(Vertex { 10, 10, 10 }, 0).elements().size());

    {
        Timer timer("full traversal");
        traverse(octree.root());
    }

    size_t totalPoints = octree.root().elements().size();
    {
        for (int level = 0; level < octree.root().getLevel(); ++level) {
          Timer timer("level traversal");
          size_t levelPoints = 0;
          for (const auto node : octree.nodesAtLevel(16)) {
              levelPoints += node.elements().size();
              //printf("lev %d: non-empty node with %zu elements\n", node.getLevel(), node.elements().size());
          }
          printf("level %d: %zu points\n", level, levelPoints);
        }
    }

    return 0;
}
