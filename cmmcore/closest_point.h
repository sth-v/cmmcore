//
// Created by Andrew Astakhov on 11.10.24.
//

#ifndef CLOSEST_POINT_H
#define CLOSEST_POINT_H


#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#ifdef CYTHON_ABI
#include "vec.h"
#else
#include "cmmcore/vec.h"
#endif
namespace cmmcore{

    inline vec3 closestPointOnLine(const vec3& a, const vec3& b) {
        vec3 ab = b - a;
        float t = -(a.dot(ab)) / ab.sqLength();
        t = std::clamp<double>(t, 0.0, 1.0);
        return a + ab * t;
    }


    inline vec3 closestPointOnTriangle(const vec3& a, const vec3& b, const vec3& c) {
        // Barycentric coordinates method
        vec3 ab = b - a;
        vec3 ac = c - a;
        vec3 ao = -a;

        float d1 = ab.dot(ao);
        float d2 = ac.dot(ao);
        if (d1 <= 0 && d2 <= 0) return a; // Vertex region A

        vec3 bo = -b;
        float d3 = ab.dot(bo);
        float d4 = ac.dot(bo);
        if (d3 >= 0 && d4 <= d3) return b; // Vertex region B

        float vc = d1 * d4 - d3 * d2;
        if (vc <= 0 && d1 >= 0 && d3 <= 0) {
            float v = d1 / (d1 - d3);
            return a + ab * v; // Edge region AB
        }

        vec3 co = -c;
        float d5 = ab.dot(co);
        float d6 = ac.dot(co);
        if (d6 >= 0 && d5 <= d6) return c; // Vertex region C

        float vb = d5 * d2 - d1 * d6;
        if (vb <= 0 && d2 >= 0 && d6 <= 0) {
            float w = d2 / (d2 - d6);
            return a + ac * w; // Edge region AC
        }

        float va = d3 * d6 - d5 * d4;
        if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
            float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return b + (c - b) * w; // Edge region BC
        }

        // Inside face region
        float denom = 1.0f / (va + vb + vc);
        float v = vb * denom;
        float w = vc * denom;
        return a + ab * v + ac * w;
    }
    // closestPoint protocol support
    inline vec3 closestPoint(const vec3& a, const vec3& b){
        return closestPointOnLine(a,  b);
    }
    inline vec3 closestPoint(const vec3& a, const vec3& b, const vec3& c){
        return closestPointOnTriangle(a,  b,  c);
    }

}
#endif //CLOSEST_POINT_H
