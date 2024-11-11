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
    inline void closestPointOnRay(const vec3& a, const vec3& start,const vec3& rayDirection,vec3& result )
    {
        vectorProjection( a-start,rayDirection,result);
        result+=start;
    }
    inline vec3 closestPointOnRay(const vec3& a, const vec3& start,const vec3& rayDirection)
    {
        vec3 result;
        closestPointOnRay(a,start,rayDirection,result);
        return result;

    }
    inline void closestPointOnLine(const vec3& a, const vec3& start,const vec3& end,vec3& result)
    {


        vectorProjection(a-start,  end-start,result);
        result+=start ;
    }
    inline vec3 closestPointOnLine(const vec3& a, const vec3& start,const vec3& end)
    {

        vec3 result;
        closestPointOnLine(a,start,end,result);
        return result;
    }
    inline void closestPointOnSegm(vec3& a, vec3& start, vec3& end,vec3& result)
    {
        vec3 d=end-start;
        vec3 o=a-start;
        double t=d.dot(o);

        if (t<0.)
        {
            result.set(start);

        } else if (t>1.)
        {
            result.set(end);
        } else
        {    vectorProjection(o,  d,result);
            result+=start ;

        }

    }
    inline vec3 closestPointOnSegm(vec3& a, vec3& start, vec3& end)
    {
        vec3 res;
        closestPointOnSegm(a,start,end,res);
        return res;
    }
    // Closest point on the segment (a,b) to the point (0.,0.,0.)
    inline vec3 closestPointOnSegm(const vec3& a, const vec3& b) {
        const vec3 ab = b - a;
        double t = -(a.dot(ab)) / ab.sqLength();
        t = std::clamp<double>(t, 0.0, 1.0);
        return a + ab * t;
    }

    // Closest point on the triangle (a, b, c) to the point (0.,0.,0.)
    inline vec3 closestPointOnTriangle(const vec3& a, const vec3& b, const vec3& c) {
        // Barycentric coordinates method
        vec3 ab = b - a;
        vec3 ac = c - a;
        vec3 ao = -a;

        double d1 = ab.dot(ao);
        double d2 = ac.dot(ao);
        if (d1 <= 0 && d2 <= 0) return a; // Vertex region A

        vec3 bo = -b;
        double d3 = ab.dot(bo);
        double d4 = ac.dot(bo);
        if (d3 >= 0 && d4 <= d3) return b; // Vertex region B

        double vc = d1 * d4 - d3 * d2;
        if (vc <= 0 && d1 >= 0 && d3 <= 0) {
            double v = d1 / (d1 - d3);
            return a + ab * v; // Edge region AB
        }

        vec3 co = -c;
        double d5 = ab.dot(co);
        double d6 = ac.dot(co);
        if (d6 >= 0 && d5 <= d6) return c; // Vertex region C

        double vb = d5 * d2 - d1 * d6;
        if (vb <= 0 && d2 >= 0 && d6 <= 0) {
            double w = d2 / (d2 - d6);
            return a + ac * w; // Edge region AC
        }

        double va = d3 * d6 - d5 * d4;
        if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
            double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return b + (c - b) * w; // Edge region BC
        }

        // Inside face region
        double denom = 1.0f / (va + vb + vc);
        double v = vb * denom;
        double w = vc * denom;
        return a + ab * v + ac * w;
    }
    inline void closestPointOnTriangle(const vec3& pt,const vec3& a, const vec3& b, const vec3& c,vec3& result)
    {
        result.set(closestPointOnTriangle(a-pt,b-pt,c-pt));

    }
    inline vec3 closestPointOnTriangle(const vec3& pt,const vec3& a, const vec3& b, const vec3& c)
    {
        vec3 res;
        closestPointOnTriangle(pt,a,b,c,res);
        return res;
    }



}
#endif //CLOSEST_POINT_H
