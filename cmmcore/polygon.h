//
// Created by Andrew Astakhov on 08.10.24.
//

#ifndef POLYGON_H
#define POLYGON_H
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <unordered_set>
#ifdef CYTHON_ABI
#include "vec.h"
#else
#include "cmmcore/vec.h"
#endif

namespace cmmcore{
    using polygon2=std::vector<vec2>;
    struct Normal2Hash {
        bool operator()(const vec2& a) const {
            auto _a=a.unit();

            auto _x=quantize(_a.x,std::numeric_limits<double>::digits10);
            auto _y=quantize(_a.y,std::numeric_limits<double>::digits10);
            auto _x1=quantize(-_a.x,std::numeric_limits<double>::digits10);
            auto _y1=quantize(-_a.y,std::numeric_limits<double>::digits10);

            return   std::max(_x ^ (_y << 1),    _x1 ^ (_y1 << 1));

        };
    };
    template<typename T>
    struct NormalEqual {
        bool operator()(const T& a, const T& b) const {
            return a.collinear(b);

        };
    };

    using Normal2Equal=NormalEqual<vec2>;
    using UnorderedNormal2Set=std::unordered_set<vec2,Normal2Hash,Normal2Equal>;


    inline bool contains( const UnorderedNormal2Set& s, const vec2& val) {
        return std::find(s.begin(), s.end(), val) != s.end();
    }

    enum class PolygonRelationship2D {
        INTERSECT,
        TOUCH,
        DISTANCE
    };



inline bool pointInPolygon2D(const polygon2& polygon, const vec2& p) {
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

inline PolygonRelationship2D polygonRelationship2D(const polygon2& poly1, const polygon2& poly2) {
    auto projectPolygon = [](const polygon2& poly, const vec2& axis) {
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
        [&poly2](const vec2& p) { return pointInPolygon2D(poly2, p); });
    bool poly2InsidePoly1 = std::all_of(poly2.begin(), poly2.end(),
        [&poly1](const vec2& p) { return pointInPolygon2D(poly1, p); });

    if (poly1InsidePoly2 || poly2InsidePoly1) {
        return PolygonRelationship2D::INTERSECT;
    }

    return touching ? PolygonRelationship2D::TOUCH : PolygonRelationship2D::INTERSECT;
}

    inline void polygonNormals( const polygon2& polygon, std::vector<vec2>& normals, const bool normalize = false ) {
        normals.reserve( polygon.size() );
        for( size_t i = 0; i < polygon.size(); i++ ) {
            size_t  j = ( i + 1 ) % polygon.size();
            normals.emplace_back( polygon[j].y - polygon[i].y, polygon[i].x - polygon[j].x );
            if (normalize) {
                normals[i].unitize();
            }
        }
    };


    inline void polygonUniqueNormals( const polygon2& polygon, UnorderedNormal2Set& normals, const bool normalize = false ) {

        for( size_t i = 0; i < polygon.size(); i++ ) {
            size_t  j = ( i + 1 ) % polygon.size();
            vec2 n={polygon[j].y - polygon[i].y, polygon[i].x - polygon[j].x};

            if (contains(normals,n)) {
                continue;
            } else {
                if (normalize) {
                    n.unitize();
                }
                normals.insert( n);

            }
        }
    };



}



#endif //POLYGON_H
