//
// Created by Andrew Astakhov on 11.10.24.
//

#ifndef CMMCORE_NUMERIC_UTILS_H
#define CMMCORE_NUMERIC_UTILS_H
#ifdef CYTHON_ABI
#include "vec.h"

#else
#include "cmmcore/vec.h"

#endif
#include <cstring>
namespace cmmcore {
inline double segmentsDistance(const vec3& p1,const vec3& p2,const vec3& q1,const vec3& q2,
                      double& t,double& u, vec3& cpt1, vec3& cpt2 ){
  auto d1=p2-p1;
  auto d2=q2-q1;
  auto r=p1-q1;
  double a,b,c,d,e,denom;
  a=d1.dot(d1);
  b = d1.dot( d2)   ;
  c = d2.dot( d2)   ;
  d = d1.dot( r)    ;
  e = d2.dot( r)    ;
    denom = a * c - b * b;
   if (denom != 0){
            t = (b * e - c * d) / denom;
            u = (a * e - b * d) / denom;
} else {
  t = u = 0;
}
    p1.add(t * d1,cpt1);
    q1.add(u * d2,cpt2);
    auto diffvec=cpt2-cpt1;
    return diffvec.length();
}
    inline double segmentsDistance(const vec3& p1,const vec3& p2,const vec3& q1,const vec3& q2
                       ){
    auto d1=p2-p1;
    auto d2=q2-q1;
    auto r=p1-q1;
    double a,b,c,d,e,denom,t,u;
    a=d1.dot(d1);
    b = d1.dot( d2)   ;
    c = d2.dot( d2)   ;
    d = d1.dot( r)    ;
    e = d2.dot( r)    ;
    denom = a * c - b * b;
    if (denom != 0){
        t = (b * e - c * d) / denom;
        u = (a * e - b * d) / denom;
    } else {
        t = u = 0;
    }
    auto cpt1=p1+t * d1;
    auto cpt2=q1+u * d2;
    return (cpt2-cpt1).length();
}

    inline void centerProjectionToXY(const vec3& origin, const vec3& point, vec2& projected) {
    if (origin[2] == point[2]) {
        // If the z-coordinates are the same, no projection is needed
        projected.set(point.x,point.y);

        return;

    }

    // Calculate the scaling factor
    const double t = -origin[2] / (point[2] - origin[2]);

    // Calculate the projected point

    projected[0] = origin[0] + t * (point[0] - origin[0]);
    projected[1] = origin[1] + t * (point[1] - origin[1]);



}
}
#endif //CMMCORE_NUMERIC_UTILS_H
