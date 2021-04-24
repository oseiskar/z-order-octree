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
    };

    struct Workspace {
        std::vector<size_t> order;
        std::vector<ZIndex> zindices;
        std::vector<const Element*> elements;
    };

    ZOrderOctree(const Element* elementsBegin, const Element *elementsEnd, const Parameters &params, Workspace *workspace = nullptr)
    :
        minCorner(saxpy(-params.leafSize * HALF_MAX_COORD, { 1, 1, 1 }, params.origin)),
        leafSize(params.leafSize)
    {
        // log_debug("minCorner %g,%g,%g", minCorner[0], minCorner[1], minCorner[2]);
        Workspace workNew;
        Workspace *tmp = workspace;
        if (!tmp) tmp = &workNew;
        tmp->elements.clear();
        tmp->zindices.clear();
        tmp->order.clear();

        for (const auto *itr = elementsBegin; itr != elementsEnd; ++itr) {
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
        for (size_t idx : tmp->order) {
            const auto zidx = tmp->zindices.at(idx);
            if (zidx == INVALID_COORD) break;
            else assert(zidx < INVALID_COORD);
            // log_debug("order: %zu (zindex %lx)", idx, zidx);
            zindices.push_back(tmp->zindices.at(idx));
            elements.push_back(tmp->elements.at(idx));
        }
    }

    struct Range {
        using Iterator = const Element**;
        Iterator b = nullptr, e = nullptr;
        Iterator begin() const { return b; }
        Iterator end() const { return e; }
    };

    template <class Point> Range lookup(const Point &point, int level) {
        ZIndex zindex = getZIndex(point);
        if (zindex == INVALID_COORD) return {};
        ZIndex mask = 0;
        for (int l = 0; l < level; ++l) {
            mask = (mask << 3l) | 0x7l;
        }
        mask = std::numeric_limits<ZIndex>::max() ^ mask;
        zindex = zindex & mask;

        // log_debug("lookup %g,%g,%g -> zindex %lx, mask %lx, sz %g^3", point[0], point[1], point[2], zindex, mask, leafSize * (1 << level));

        Range r;
        r.b = elements.data() + findRange(zindex & mask, mask, true);
        r.e = elements.data() + findRange(zindex & mask, mask, false);
        return r;
    }

private:
    static Vector3 saxpy(Float alpha, const Vector3 &x, const Vector3 &y) {
        Vector3 r;
        for (size_t i = 0; i < 3; ++i) r[i] = alpha * x[i] + y[i];
        return r;
    }

    static constexpr ZIndex BITS_PER_COORD = std::numeric_limits<ZIndex>::digits / 3;
    static constexpr ZIndex HALF_MAX_COORD = 1 << (BITS_PER_COORD - 1);
    static constexpr ZIndex MAX_COORD = 2 * HALF_MAX_COORD;
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
        for (int d = 0; d < 3; ++d) {
            int coord = (xyz[d] - minCorner[d]) / leafSize;
            // log_debug("getZIndex, coord. %d: %g -> %d", d, xyz[d], coord);
            if (coord < 0 || coord >= MAX_COORD) return INVALID_COORD;
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

    const Vector3 minCorner;
    const float leafSize;
    std::vector<size_t> zindices;
    std::vector<const Element*> elements;
};
