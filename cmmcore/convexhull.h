//
// Created by Andrew Astakhov on 08.10.24.
//

#ifndef CONVEXHULL_H
#define CONVEXHULL_H
#include "_convexhull/convexhull2d.h"
#include "_convexhull/quickhull.hpp"
#include "cmmcore/plane.h"
#include <cstring>
namespace cmmcore {
    enum class ConvexHullResult {
        VOLUME, FLAT, LINE, POINT, EMPTY
    };

    inline ConvexHullResult convex_hull(const std::vector<vec3> &pts1, std::vector<vec3>& pts1hull) {
        if (pts1.empty()) {
            return ConvexHullResult::EMPTY;
        }
        const size_t pts_length=pts1.size();
        if (isPtsCoplanar(pts1)) {
            if (pts_length> 3) {
                std::vector<vec2> _pt1(pts_length);
                for (size_t i = 0; i < pts_length; i++) {
                    _pt1[i].set(pts1[i].x, pts1[i].y);
                }
                auto hull12d = convex_hull2d(_pt1);
                pts1hull.reserve(hull12d.size());
                for (auto &i: hull12d) {
                    auto res = std::find(pts1.begin(), pts1.end(), i);
                    pts1hull.emplace_back(*res);
                }
            } else {
                pts1hull.resize(pts_length);
                   memcpy(pts1hull.data(), pts1.data(),  pts_length * sizeof(vec3));

            }
        } else {
            convex_hull3d(pts1, pts1hull);
        };
        switch (pts1hull.size()) {
            case 1: return ConvexHullResult::POINT;
            case 2: return ConvexHullResult::LINE;
            case 3: return ConvexHullResult::FLAT;
            default: return ConvexHullResult::VOLUME;

        }
    }

    inline ConvexHullResult convex_hull(const std::vector<vec2> &pts1,std::vector<vec2>& pts1hull) {
        if (pts1.empty()) {
            return ConvexHullResult::EMPTY;
        }
        std::vector<vec2> pts = pts1;
        pts1hull = convex_hull2d(pts);

        switch (pts1hull.size()) {
            case 1: return ConvexHullResult::POINT;
            case 2: return ConvexHullResult::LINE;
            default: return ConvexHullResult::FLAT;


        }
    }
}
#endif //CONVEXHULL_H
