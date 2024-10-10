//
// Created by Andrew Astakhov on 08.10.24.
//

#ifndef POLYGON_H
#define POLYGON_H
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

#ifdef CYTHON_ABI
#include "vec.h"
#else
#include "cmmcore/vec.h"
#endif

namespace cmmcore{


enum class PolygonRelationship2D {
    INTERSECT,
    TOUCH,
    DISTANCE
};



bool isvec2Inside2d(const std::vector<vec2>& polygon, const vec2& p) {
    int n = polygon.size();
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        if (((polygon[i].y <= p.y && p.y < polygon[j].y) ||
             (polygon[j].y <= p.y && p.y < polygon[i].y)) &&
            (p.x < (polygon[j].x - polygon[i].x) * (p.y - polygon[i].y) /
                   (polygon[j].y - polygon[i].y) + polygon[i].x)) {
            inside = !inside;
        }
    }
    return inside;
}

PolygonRelationship2D polygonRelationship2D(const std::vector<vec2>& poly1, const std::vector<vec2>& poly2) {
    auto projectPolygon = [](const std::vector<vec2>& poly, const vec2& axis) {
        double min = std::numeric_limits<double>::max();
        double max = std::numeric_limits<double>::lowest();
        for (const auto& p : poly) {
            double proj = dot(vec2(p.x, p.y), axis);
            min = std::min(min, proj);
            max = std::max(max, proj);
        }
        return std::make_pair(min, max);
    };

    auto checkOverlap = [](std::pair<double, double> proj1, std::pair<double, double> proj2) {
        return proj1.first <= proj2.second && proj2.first <= proj1.second;
    };

    std::vector<vec2> axes;
    for (int i = 0; i < poly1.size(); i++) {
        vec2 edge(poly1[i], poly1[(i + 1) % poly1.size()]);
        axes.push_back(vec2(-edge.y, edge.x));
    }
    for (int i = 0; i < poly2.size(); i++) {
        vec2 edge(poly2[i], poly2[(i + 1) % poly2.size()]);
        axes.push_back(vec2(-edge.y, edge.x));
    }

    bool touching = false;
    for (const auto& axis : axes) {
        auto proj1 = projectPolygon(poly1, axis);
        auto proj2 = projectPolygon(poly2, axis);

        if (!checkOverlap(proj1, proj2)) {
            return PolygonRelationship2D::DISTANCE;
        }

        if (proj1.first == proj2.second || proj1.second == proj2.first) {
            touching = true;
        }
    }

    // Check if one polygon is inside the other
    bool poly1InsidePoly2 = std::all_of(poly1.begin(), poly1.end(),
        [&poly2](const vec2& p) { return isvec2Inside2d(poly2, p); });
    bool poly2InsidePoly1 = std::all_of(poly2.begin(), poly2.end(),
        [&poly1](const vec2& p) { return isvec2Inside2d(poly1, p); });

    if (poly1InsidePoly2 || poly2InsidePoly1) {
        return PolygonRelationship2D::INTERSECT;
    }

    return touching ? PolygonRelationship2D::TOUCH : PolygonRelationship2D::INTERSECT;
}
}


#endif //POLYGON_H
