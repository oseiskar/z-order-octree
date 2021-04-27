#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#ifdef Z_ORDER_OCTREE_DEBUG
#include <cstdio>
#define log_debug(...) do { std::printf(__VA_ARGS__); std::printf("\n"); } while (0)
#else
#define log_debug(...)
#endif

template <class Element, class Float = float, class ZIndex = std::uint64_t>
class ZOrderOctree {
public:
    // Note: "should" also work with 1 more bit but something is broken
    static constexpr size_t MAX_ROOT_LEVEL = std::numeric_limits<ZIndex>::digits / 3 - 1;
    using Vector3 = std::array<Float, 3>;
    struct Parameters {
        Vector3 origin = { 0, 0, 0 };
        Float leafSize = 1.0;
        bool stableSort = false;
        size_t rootLevel = MAX_ROOT_LEVEL;
    };

    ZOrderOctree(const Parameters &params) :
        params(params)
    {
        assert(params.rootLevel <= MAX_ROOT_LEVEL);
    }

    void clear() {
        zindices.clear();
        elements.clear();
    }

    void addData(const Element* elementsBegin, size_t nElements) {
        // log_debug("minCorner %g,%g,%g", minCorner[0], minCorner[1], minCorner[2]);
        tmp.clear();
        tmp.reserve(nElements + zindices.size());

        assert(zindices.size() == elements.size());

        // add existing data (note: it could be faster to do a custom merge-sort)
        for (size_t i = 0; i < zindices.size(); ++i) {
            tmp.zindices.push_back(zindices.at(i));
            tmp.elements.push_back(elements.at(i));
            tmp.order.push_back(tmp.order.size());
        }

        zindices.clear();
        elements.clear();

        for (const auto *itr = elementsBegin; itr != (elementsBegin + nElements); ++itr) {
            tmp.elements.push_back(&*itr);
            tmp.zindices.push_back(getZIndex(*itr));
            tmp.order.push_back(tmp.order.size());
        }

        assert(tmp.order.size() == tmp.zindices.size());
        const auto cmp = [this](size_t a, size_t b) -> int {
            return tmp.zindices.at(a) < tmp.zindices.at(b);
        };
        if (params.stableSort) {
            std::stable_sort(tmp.order.begin(), tmp.order.end(), cmp);
        } else {
            std::sort(tmp.order.begin(), tmp.order.end(), cmp);
        }

        elements.reserve(tmp.elements.size());
        zindices.reserve(tmp.zindices.size());
        for (size_t idx : tmp.order) {
            const auto zidx = tmp.zindices.at(idx);
            if (zidx == INVALID_COORD) break;
            else assert(zidx < INVALID_COORD);
            // log_debug("order: %zu (zindex %lx)", idx, zidx);
            zindices.push_back(tmp.zindices.at(idx));
            elements.push_back(tmp.elements.at(idx));
        }
    }

    template <class Predicate> void removeData(const Predicate &func) {
        auto zIdxIn = zindices.begin();
        auto elemIn = elements.begin();
        auto zIdxOut = zIdxIn;
        auto elemOut = elemIn;
        while (elemIn != elements.end()) {
            if (!func(**elemIn)) {
              *elemOut++ = *elemIn;
              *zIdxOut++ = *zIdxIn;
            }
            elemIn++;
            zIdxIn++;
        }
        const size_t nLeft = elemOut - elements.begin();
        elements.resize(nLeft);
        zindices.resize(nLeft);
    }

    void clearWorkspace() {
        tmp = {}; // should deallocate
    }

    class ElementRange {
    public:
        using Iterator = const Element * const *;
        Iterator begin() const { return b; }
        Iterator end() const { return e; }

        size_t size() const {
            return e - b;
        }

        ElementRange(Iterator b, Iterator e) : b(b), e(e) {}

    private:
        Iterator b = nullptr, e = nullptr;
    };

    class NodeRange;

    struct Node {
        Node(const ZOrderOctree &t, ZIndex zidx, int level) :
            tree(&t),
            level(level),
            zindex(zidx)
        {
            ZIndex mask = levelMask(level);
            zindex = zindex & mask;
            elementsBegin = tree->findRange(zindex, mask, true);
            elementsEnd = tree->findRange(zindex, mask, false);
            // log_debug("node %lx at level %d, elements %zu to %zu", zindex, level, elementsBegin, elementsEnd);
        }

        // "end node" marker
        Node() :
          tree(nullptr),
          level(-1),
          zindex(0),
          elementsBegin(0),
          elementsEnd(0)
        {}

        int getLevel() const {
            return level;
        }

        bool isLeaf() const {
            return level == 0 || isEndNode() || empty();
        }

        Node firstChild() const {
            assert(!isLeaf());
            return Node(*tree, tree->zindices.at(elementsBegin), level - 1);
        }

        bool isLastAtThisLevel() const {
            if (isEndNode()) return true;
            return elementsEnd == tree->zindices.size();
        }

        bool isLastSibling() const {
            if (isLastAtThisLevel()) return true;
            ZIndex parentMask = levelMask(level + 1);
            return (zindex & parentMask) != (tree->zindices.at(elementsEnd) & parentMask);
        }

        Node nextAtThisLevel() const {
            assert(!isLastAtThisLevel());
            return Node(*tree, tree->zindices.at(elementsEnd), level);
        }

        Node nextSibling() const {
            assert(!isLastSibling());
            return nextAtThisLevel();
        }

        bool isEndNode() const {
            return tree == nullptr;
        }

        bool operator==(const Node &other) const {
            if (isEndNode()) return other.isEndNode();
            return tree == other.tree &&
                level == other.level &&
                zindex == other.zindex;
        }

        NodeRange children() const {
            assert(!isLeaf());
            Node child = firstChild();
            assert(!child.empty());
            return NodeRange(child, true);
        }

        ElementRange elements() const {
            assert(!isEndNode());
            return tree->buildRange(elementsBegin, elementsEnd);
        }

        bool empty() const {
            return elementsEnd == elementsBegin;
        }

        Vector3 minCorner() const {
            return tree->zIndexToPoint(zindex, level, 0);
        }

        Vector3 maxCorner() const {
            return tree->zIndexToPoint(zindex, level, 1);
        }

        Vector3 center() const {
            return tree->zIndexToPoint(zindex, level, 0.5);
        }

        Float sideLength() const {
            return tree->params.leafSize * (1 << level);
        }

        bool isRoot() const {
            return level == tree->params.rootLevel;
        }

        // internal, do not use directly
        const ZOrderOctree *tree;
        int level;
        ZIndex zindex;
        size_t elementsBegin, elementsEnd;
    };

    class LevelIterator {
    public:
        LevelIterator(const Node &node, bool siblingsOnly) : node(node), siblings(siblingsOnly) {}

        const Node &operator*() const {
            assert(!node.isEndNode());
            return node;
        }

        LevelIterator &operator++() { // prefix OP, ++itr
            if ((siblings && node.isLastSibling()) ||
                (!siblings && node.isLastAtThisLevel())) {
                node = Node();
            } else {
                node = node.nextAtThisLevel();
            }
            return *this;
        }

        bool operator==(const LevelIterator &other) const {
            return node == other.node && siblings == other.siblings;
        }

        bool operator!=(const LevelIterator &other) const {
            return !(*this == other);
        }

    private:
        Node node;
        const bool siblings;
    };

    class NodeRange {
    private:
        LevelIterator b, e;

    public:
        NodeRange(Node beginNode, bool siblingsOnly) :
          b(beginNode, siblingsOnly),
          e(Node(), siblingsOnly)
        {}
        LevelIterator begin() const { return b; }
        LevelIterator end() const { return e; }
        bool empty() const { return b == e; }
    };

    template<class Point> Node lookup(const Point &point, int level) const {
        assert(level >= 0 && level < params.rootLevel);
        ZIndex zindex = getZIndex(point);
        if (zindex == INVALID_COORD) return Node();

        return Node(*this, zindex, level);
    }

    Node root() const {
        return Node(*this, zindices.empty() ? 0 : *zindices.begin(), params.rootLevel);
    }

    NodeRange nodesAtLevel(int level) const {
        assert(level >= 0 && level < params.rootLevel);
        return NodeRange(Node(*this, zindices.empty() ? 0 : *zindices.begin(), level), false);
    }

    class RadiusSearchIterator {
    public:
        const Element *operator*() const {
            assert(currentNodeIdx < nodeCount);
            return nodes[currentNodeIdx].tree->elements.at(currentElementIdx);
        }

        RadiusSearchIterator &operator++() { // prefix OP, ++itr
            assert(currentNodeIdx < nodeCount && currentNodeIdx < 8);
            findNext(false);
            return *this;
        }

        // note: do not compare these from two unrelated ranges
        bool operator==(const RadiusSearchIterator &other) const {
            return currentNodeIdx == other.currentNodeIdx &&
                currentElementIdx == other.currentElementIdx;
        }

        bool operator!=(const RadiusSearchIterator &other) const {
            return !(*this == other);
        }

        RadiusSearchIterator(
            size_t nodeCount,
            bool end = true,
            Float searchRadius = 0,
            const Vector3 *searchCenter = nullptr,
            Node *nodesPtr = nullptr)
        :
            nodeCount(nodeCount),
            currentNodeIdx(nodeCount),
            currentElementIdx(0)
        {
            assert(nodeCount <= 8);
            if (!end) {
                searchRadiusSquared = searchRadius*searchRadius;
                this->searchCenter = *searchCenter;
                for (size_t i = 0; i < nodeCount; ++i) nodes[i] = nodesPtr[i];
                currentNodeIdx = 0;
                if (nodeCount > 0) {
                    currentElementIdx = nodes[currentNodeIdx].elementsBegin;
                    findNext(true);
                }
            }
            // log_debug("iterator %zu/%zu, %zu", currentNodeIdx, nodeCount, currentElementIdx);
        }

        size_t getNodeCount() const {
            return nodeCount;
        }

    private:
        void findNext(bool checkFirst) {
            if (checkFirst && (
                currentNodeIdx >= nodeCount ||
                currentElementIdx >= nodes[currentNodeIdx].elementsEnd)) return;

            bool checkCurrent = checkFirst;
            while (true) {
                // log_debug("search %zu/%zu, %zu in %zu - %zu", currentNodeIdx, nodeCount, currentElementIdx,
                //    nodes[currentNodeIdx].elementsBegin,  nodes[currentNodeIdx].elementsEnd);
                if (checkCurrent) {
                    const Element &el = *nodes[currentNodeIdx].tree->elements.at(currentElementIdx);
                    Float r2 = 0;
                    for (int c = 0; c < 3; ++c) {
                        Float d = el[c] - searchCenter[c];
                        r2 += d*d;
                    }
                    // log_debug("%g/%g (%g, %g, %g)", std::sqrt(r2), std::sqrt(searchRadiusSquared), el[0], el[1], el[2]);
                    if (r2 < searchRadiusSquared) break;
                }
                checkCurrent = true;

                if (++currentElementIdx >= nodes[currentNodeIdx].elementsEnd) {
                    if (++currentNodeIdx < nodeCount) {
                        currentElementIdx = nodes[currentNodeIdx].elementsBegin;
                    } else {
                        // end node
                        currentElementIdx = 0;
                        break;
                    }
                }
            }
        }

        size_t currentNodeIdx, currentElementIdx;
        size_t nodeCount;
        Node nodes[8];
        Vector3 searchCenter;
        Float searchRadiusSquared;
    };

    class RadiusSearchRange {
    private:
        RadiusSearchIterator b, e;

    public:
        RadiusSearchRange(RadiusSearchIterator begin) :
          b(begin), e(begin.getNodeCount())
        {}

        RadiusSearchIterator begin() const { return b; }
        RadiusSearchIterator end() const { return e; }
        bool empty() const { return b == e; }
    };

    template<class Point> RadiusSearchRange searchWithRadius(const Point &point, Float radius) const {
        // TODO: breaks / overflows with very large radii like 1e10

        int searchLevel = 0;
        const Vector3 searchCenter = { point[0], point[1], point[2] };
        while (nodeCountWithinBox(searchCenter, radius, searchLevel) > 8) {
            assert(searchLevel < params.rootLevel);
            searchLevel++;
        }

        Node nodes[8];
        size_t nodeCount = 0;
        for (int dx = -1; dx <= 1; dx += 2) {
            for (int dy = -1; dy <= 1; dy += 2) {
                for (int dz = -1; dz <= 1; dz += 2) {
                    Vector3 corner = searchCenter;
                    corner[0] += dx * radius;
                    corner[1] += dy * radius;
                    corner[2] += dz * radius;
                    ZIndex zindex = getZIndex(corner, true) & levelMask(searchLevel);
                    // log_debug("corner (%g, %g, %g) zindex %lx", corner[0], corner[1], corner[2], zindex);
                    for (size_t j = 0; j < nodeCount; ++j) {
                        if (zindex == nodes[j].zindex) {
                            // log_debug("equal to corner %zu", j);
                            goto NEXT_CORNER;
                        }
                    }

                    assert(nodeCount < 8);
                    nodes[nodeCount++] = Node(*this, zindex, searchLevel);

                    NEXT_CORNER: (void)0;
                }
            }
        }
        // assert(nodeCount == nodeCountWithinBox(searchCenter, radius, searchLevel));
        return RadiusSearchRange(RadiusSearchIterator(nodeCount, false, radius, &searchCenter, nodes));
    }

private:
    struct Workspace {
        std::vector<size_t> order;
        std::vector<ZIndex> zindices;
        std::vector<const Element*> elements;

        void clear() {
            order.clear();
            zindices.clear();
            elements.clear();
        }

        void reserve(size_t n) {
            order.reserve(n);
            zindices.reserve(n);
            elements.reserve(n);
        }
    } tmp;

    ElementRange buildRange(size_t elementsBegin, size_t elementsEnd) const {
        return ElementRange(
          elements.data() + elementsBegin,
          elements.data() + elementsEnd);
    }

    static ZIndex levelMask(int level) {
        ZIndex mask = 0;
        for (int l = 0; l < level; ++l) {
            mask = (mask << 3l) | 0x7l;
        }
        return std::numeric_limits<ZIndex>::max() ^ mask;
    }

    static Vector3 saxpy(Float alpha, const Vector3 &x, const Vector3 &y) {
        Vector3 r;
        for (size_t i = 0; i < 3; ++i) r[i] = alpha * x[i] + y[i];
        return r;
    }

    static constexpr ZIndex INVALID_COORD = std::numeric_limits<ZIndex>::max();

    size_t findRange(ZIndex target, ZIndex mask, bool findBegin) const {
        // binary search
        size_t begin = 0;
        size_t end = zindices.size();
        while (end > begin) {
            size_t mid = begin + (end - begin) / 2;
            ZIndex cur = zindices[mid] & mask;
            // log_debug("findRange [%zu, %zu, %zu] -> [?, %lx, ?]", begin, mid, end, cur);
            if ((findBegin && (cur < target)) || (!findBegin && (cur <= target))) {
                if (begin == mid) mid++;
                begin = mid;
            } else {
                end = mid;
            }
        }
        return end;
    }

    inline int maxCoord() const {
        return 1 << params.rootLevel;
    }

    inline int halfMaxCoord() const {
        return (1 << params.rootLevel) / 2;
    }

    inline int floatToCoord(Float c, int dim) const {
        return std::floor((c - params.origin[dim]) / params.leafSize) + halfMaxCoord();
    }

    Vector3 zIndexToPoint(ZIndex zindex, int level, Float cellOffset) const {
        int coords[3] = { 0, 0, 0 };
        for (int l = params.rootLevel; l >= level; --l) {
            for (int d = 0; d < 3; ++d) {
                int bit = (zindex >> (3*l + d)) & 0x1;
                if (bit) coords[d] += 1 << l;
            }
        }
        Vector3 v;
        const Float offs = params.leafSize * (1 << level) * cellOffset;
        for (int d = 0; d < 3; ++d) {
            v[d] = (coords[d] - halfMaxCoord()) * params.leafSize + params.origin[d] + offs;
        }
        return v;
    }

    template <class Point> ZIndex getZIndex(const Point &xyz, bool capped = false) const {
        ZIndex zindex = 0;
        for (int d = 0; d < 3; ++d) {
            int coord = floatToCoord(xyz[d], d);
            // log_debug("getZIndex, coord. %d: %g -> %d", d, xyz[d], coord);
            if (coord < 0) {
                if (capped) coord = 0;
                return INVALID_COORD;
            }
            else if (coord >= maxCoord()) {
                if (capped) coord = maxCoord() - 1;
                return INVALID_COORD;
            }
            zindex |= interleaveBits(coord) << d;
        }
        // log_debug("getZIndex -> %lx", zindex);
        return zindex;
    }

    size_t nodeCountWithinBox(const Vector3 &center, Float radius, int level) const {
        size_t count = 1;
        for (int d = 0; d < 3; ++d) {
            int range = 1;
            for (int sign = -1; sign <= 1; sign += 2) {
                int coord = floatToCoord(center[d] + sign*radius, d);
                if (coord < 0) coord = 0;
                if (coord >= maxCoord()) coord = maxCoord() - 1;
                coord = coord >> level;
                range += coord * sign;
            }
            assert(range >= 0);
            count *= range;
        }
        // log_debug("nodeCountWithinBox %g, %d -> %zu", radius, level, count);
        return count;
    }

    static ZIndex interleaveBits(ZIndex coord) {
        auto x = static_cast<std::uint64_t>(coord);
        // https://stackoverflow.com/a/18528775/1426569
        x &= 0x1fffff;
        x = (x | x << 32) & 0x1f00000000ffffll;
        x = (x | x << 16) & 0x1f0000ff0000ffll;
        x = (x | x << 8) & 0x100f00f00f00f00fll;
        x = (x | x << 4) & 0x10c30c30c30c30c3ll;
        x = (x | x << 2) & 0x1249249249249249ll;
        return static_cast<ZIndex>(x);
    }

    const Parameters params;
    std::vector<size_t> zindices;
    std::vector<const Element*> elements;
};
