//
// Created by Andrew Astakhov on 26.10.24.
//

#ifndef CSX_H
#define CSX_H
#include <vector>
#include <array>
#include <string>
#include <tuple>
#include <optional>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "separability.h"
#include "cmmcore/nurbs.h"
#include "cmmcore/vec.h"
#include "cmmcore/newthon.h"
namespace cmmcore {

class NURBSCurveSurfaceIntersector {
public:
    NURBSCurveSurfaceIntersector(NURBSCurve& curve, NURBSSurface& surface, double tolerance = 1e-3, double ptol = 1e-7)
        : curve(curve), surface(surface), tolerance(tolerance), ptol(ptol) {}

    std::vector<std::tuple<std::string, vec3, std::array<double, 3>>> intersect() {
        intersections.clear();
        _curve_surface_intersect(curve, surface);
        return intersections;
    }

private:
    NURBSCurve& curve;
    NURBSSurface& surface;
    double tolerance;
    double ptol;
    std::vector<std::tuple<std::string, vec3, std::array<double, 3>>> intersections;

    int _find_new_intersection(const NURBSCurve& curve, const NURBSSurface& surface,vec3& pt, std::vector<double>&tuv) {
        // Implement the logic to find a new intersection

        auto function = [&curve,&surface](const std::vector<double>& v) {
            vec3 temp1,temp2;
            curve.evaluate(v[0],temp1);
            surface.evaluate(v[1],v[2],temp2);


            return (temp1-temp2).sqLength();
        };
        vec3 tmp;
        tuv.resize(3);
        tuv[0]=(curve._interval[1]+curve._interval[0])*0.5;
        tuv[1]=(surface._interval[0][0]+surface._interval[0][1])*0.5;
        tuv[2]= (surface._interval[1][0]+surface._interval[1][1])*0.5;


        auto res=newtonsMethod(function, tuv, tolerance);
        if (res==0)
        {
            surface.evaluate(tuv[1],tuv[2],pt);
        }


        return res;



    }

    void _curve_surface_intersect(const NURBSCurve& curve, const NURBSSurface& surface) {
        if (_no_new_intersections(curve, surface)) return;
        std::pair<vec3, std::vector<double>> new_point;
        int res=_find_new_intersection(curve, surface,new_point.first,new_point.second);
        auto [u0, u1] = surface._interval[0];
        auto [v0, v1] = surface._interval[1];

        if (!(res==0)) {
            auto [t0, t1] = curve._interval;
            auto [curve1, curve2] = curve.split((t0 + t1) * 0.5);
            std::array<NURBSSurface,4> surfaces;
            surface.subdivide(surfaces[0],surfaces[1],surfaces[2],surfaces[3]);

         for (auto& s : {surfaces[0],surfaces[1], surfaces[2], surfaces[3]}) {
                for (auto& c : {curve1, curve2}) {
                    _curve_surface_intersect(c, s);
                }
            }
        } else {
            auto [point, params] = *new_point;
            auto [t, u, v] = params;

            const std::string intersection_type = _is_degenerate(params, curve, surface) ? "degenerate" : "transversal";
            intersections.emplace_back(intersection_type, point, params);

            if (std::abs(u - u0) < ptol || std::abs(u - u1) < ptol || std::abs(v - v0) < ptol || std::abs(v - v1) < ptol)
                return;


            auto [curve1, curve2] = curve.split(t);
            std::array<NURBSSurface,4> surfaces;
            surface.subdivide(surfaces[0],surfaces[1],surfaces[2],surfaces[3],u, v);
            for (size_t i = 0; i < surfaces.size(); i++)
            {

            }

            for (auto& s : {surfaces[0],surfaces[1], surfaces[2], surfaces[3]}) {
                std::vector<vec3> scpts(s.control_points.size()-1);

                auto scpts2=s.control_points_flat3d();
                size_t jc=0;
                for (size_t i = 0; i < s.control_points.size(); i++)
                {
                    scpts2[i].sub(point, scpts[jc]);
                    if (!(scpts[jc].x<=std::numeric_limits<double>::epsilon()
                    &&scpts[jc].y<=std::numeric_limits<double>::epsilon()
                    &&scpts[jc].z<=std::numeric_limits<double>::epsilon()))
                    {
                        jc++;
                    }





                }
                for (auto& c : {curve1, curve2}) {
                    std::vector<vec3> ccpts(c.control_points.size()-1);

                    std::vector<vec3> ccpts2;
                    c.get_control_points3d(ccpts2);
                    size_t ic = 0;

                    for (size_t i = 0; i < c.control_points.size(); i++)
                    {
                       ccpts2[i].sub(point, ccpts[ic]);
                        if (!(ccpts[ic].x<=std::numeric_limits<double>::epsilon()
                        &&ccpts[ic].y<=std::numeric_limits<double>::epsilon()
                        &&ccpts[ic].z<=std::numeric_limits<double>::epsilon()))
                        {
                            ic++;
                        }


                    }

                    if (!SphericalSeparabilityTest(ccpts2, scpts2))
                    {
                        _curve_surface_intersect(c, s);
                    };


                }
            }
        }
    }

    bool _no_new_intersections(const NURBSCurve& curve, const NURBSSurface& surface) const {
        std::vector<vec3> cpts;
        return SAT3D(curve.get_control_points3d(cpts), surface.control_points_flat3d(), tolerance);
    }

    bool _is_degenerate(const std::array<double, 3>& params, const NURBSCurve& curve, const NURBSSurface& surface) const {
        double t = params[0];
        double u = params[1];
        double v = params[2];

        vec3 curve_tangent;
        curve.derivative(t, curve_tangent);

        vec3 surface_normal;
        surface.normal(u, v, surface_normal);

        return std::abs(curve_tangent.dot( surface_normal)) < tolerance;
    }
};

} // namespace cmmcore
#endif //CSX_H
