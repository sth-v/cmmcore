//
// Created by Andrew Astakhov on 09.10.24.
//

#ifndef SSX_H
#define SSX_H
#include "cmmcore/nurbs.h"
#include "cmmcore/gauss_map.h"
namespace cmmcore {
    struct PointIntersection {
        vec2 uv1;
        vec2 uv2;
        vec3 xyz;
        PointIntersection()=default;
        PointIntersection(const vec2& _uv1, const vec2& _uv2, const vec3& _xyz):uv1(_uv1),uv2(_uv2),xyz(_xyz) {

        }
    };
    struct PatchIntersection {
        NURBSSurface s1;
        NURBSSurface s2;
        PatchIntersection()=default;
        PatchIntersection(NURBSSurface& _s1, NURBSSurface& _s2): s1(std::move(_s1)), s2(std::move(_s2)) {}

    };
    struct Intersection {
        std::vector<PatchIntersection> patches;
        std::vector<PointIntersection> points;
        Intersection():patches(),points() {

        }
    };
    void detectIntersections(NURBSSurface& surface1,NURBSSurface& surface2, GaussMap& gaussMap1,  GaussMap& gaussMap2, double tol, Intersection& intersection) {
        auto& bb1=surface1.bbox();
        auto& bb2=surface2.bbox();

        if (!(bb1.intersects(bb2))) {

            return;
        }
        auto bb1d=bb1.max-bb1.min;
        auto bb2d=bb2.max-bb2.min;
        if (bb1d.x<=tol&&bb1d.y<=tol&&bb1d.z<=tol&&bb2d.x<=tol&&bb2d.y<=tol&&bb2d.z<=tol) {
            intersection.points.push_back(
               {{ 0.5*(surface1._interval[0][1]+ surface1._interval[0][0]),
                    0.5*(surface1._interval[1][1]+ surface1._interval[1][0])},
                {0.5*(surface2._interval[0][1]+ surface2._interval[0][0]),
                    0.5*(surface2._interval[1][1]+ surface2._interval[1][0])},
                 (bb1.min+bb1.max+bb2.min+bb2.max)*0.25
        }
        );

        }

        AABB aabb1;
        bb1.intersection(bb2,aabb1);
        auto _dims =(aabb1.max-aabb1.min);

        if ((aabb1.volume()==0) || ( _dims.x<=tol)||(_dims.y<=tol)||(_dims.z<=tol)){
            return;
        }

        if (gaussMap1.separabilityTest(gaussMap2)) {
            intersection.patches.emplace_back(surface1,surface2);
            return;
        }
        std::array<GaussMap,4> g1;// g11,g12,g13,g14
        std::array<GaussMap,4> g2;//g21,g22,g23,g24;
         std::array<NURBSSurface ,4>s1;//s11,s12,s13,s14;
        std::array<NURBSSurface ,4>s2;//s21,s22,s23,s24;
        surface1.subdivide( s1[0],s1[1],s1[2],s1[3]);
        surface2.subdivide( s2[0],s2[1],s2[2],s2[3]);
        gaussMap1.subdivide(  g1[0],g1[1],g1[2],g1[3]);
        gaussMap2.subdivide(  g2[0],g2[1],g2[2],g2[3]);
        for (int i=0;i<4;i++) {
            for (int j=0;j<4;j++) {
                detectIntersections(s1[i],s2[j],g1[i],g2[j],tol,intersection);
            }
        }


    }
}
#endif //SSX_H
