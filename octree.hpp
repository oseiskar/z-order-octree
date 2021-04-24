#pragma once

#include <algorithm>
#include <cassert>
#include <array>
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
    using Vector3 = std::array<Float, 3>;
    struct Parameters {
        Vector3 origin = { 0, 0, 0 };
        float leafSize = 1.0;
        bool stableSort = false;
        size_t rootLevel = std::numeric_limits<ZIndex>::digits / 3;
    };

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
    };

    ZOrderOctree(const Parameters &params)
    :
        params(params),
        minCorner(saxpy(-params.leafSize * (1 << (params.rootLevel - 1)), { 1, 1, 1 }, params.origin))
    {
        assert(params.rootLevel <= std::numeric_limits<ZIndex>::digits / 3);
    }

    void clear() {
        zindices.clear();
        elements.clear();
    }

    void addData(const Element* elementsBegin, size_t nElements, Workspace *workspace = nullptr) {
        // log_debug("minCorner %g,%g,%g", minCorner[0], minCorner[1], minCorner[2]);
        Workspace workNew;
        Workspace *tmp = workspace;
        if (!tmp) tmp = &workNew;
        tmp->clear();
        tmp->reserve(nElements + zindices.size());

        assert(zindices.size() == elements.size());

        // add existing data (note: it could be faster to do a custom merge-sort)
        for (size_t i = 0; i < zindices.size(); ++i) {
            tmp->zindices.push_back(zindices.at(i));
            tmp->elements.push_back(elements.at(i));
            tmp->order.push_back(tmp->order.size());
        }

        zindices.clear();
        elements.clear();

        for (const auto *itr = elementsBegin; itr != (elementsBegin + nElements); ++itr) {
            tmp->elements.push_back(&*itr);
            tmp->zindices.push_back(getZIndex(*itr));
            tmp->order.push_back(tmp->order.size());
        }

        assert(tmp->order.size() == tmp->zindices.size());
        const auto cmp = [&tmp](size_t a, size_t b) -> int {
            return tmp->zindices.at(a) < tmp->zindices.at(b);
        };
        if (params.stableSort) {
            std::stable_sort(tmp->order.begin(), tmp->order.end(), cmp);
        } else {
            std::sort(tmp->order.begin(), tmp->order.end(), cmp);
        }

        elements.reserve(tmp->elements.size());
        zindices.reserve(tmp->zindices.size());
        for (size_t idx : tmp->order) {
            const auto zidx = tmp->zindices.at(idx);
            if (zidx == INVALID_COORD) break;
            else assert(zidx < INVALID_COORD);
            // log_debug("order: %zu (zindex %lx)", idx, zidx);
            zindices.push_back(tmp->zindices.at(idx));
            elements.push_back(tmp->elements.at(idx));
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

    Workspace buildWorkspace() const {
        return {};
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

    class Node {
    public:
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
          zindex(0),
          level(-1),
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

    private:
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

private:
    ElementRange buildRange(size_t elementsBegin, size_t elementsEnd) const {
        return ElementRange(
          elements.data() + elementsBegin,
          elements.data() + elementsEnd);
    }

    Vector3 zIndexToPoint(ZIndex zindex, int level, Float cellOffset) const {
        int coords[3] = { 0, 0, 0 };
        for (int l = params.rootLevel; l >= level; --l) {
            for (int d = 0; d < 3; ++d) {
                int bit = (zindex >> (3*l + d)) & 0x1;
                if (bit) coords[d] += 1 << level;
            }
        }
        Vector3 v;
        const Float offs = params.leafSize * (1 << level) * cellOffset;
        for (int d = 0; d < 3; ++d) {
            v[d] = coords[d] * params.leafSize + minCorner[d] + offs;
        }
        return v;
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

    template <class Point> ZIndex getZIndex(const Point &xyz) const {
        ZIndex zindex = 0;
        const int maxCoord = 1 << params.rootLevel;
        for (int d = 0; d < 3; ++d) {
            int coord = (xyz[d] - minCorner[d]) / params.leafSize;
            // log_debug("getZIndex, coord. %d: %g -> %d", d, xyz[d], coord);
            if (coord < 0 || coord >= maxCoord) return INVALID_COORD;
            zindex |= interleaveBits(coord) << d;
        }
        // log_debug("getZIndex -> %lx", zindex);
        return zindex;
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
    const Vector3 minCorner;
    std::vector<size_t> zindices;
    std::vector<const Element*> elements;
};
