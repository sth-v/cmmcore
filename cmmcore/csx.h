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

#include "separability.h"
#include "cmmcore/nurbs.h"
#include "cmmcore/vec.h"
#include "cmmcore/newthon2.h"


namespace cmmcore {




    class NURBSCurveSurfaceIntersector {
    public:
        NURBSCurveSurfaceIntersector(const NURBSCurve &curve, const NURBSSurface &surface, double tolerance = 1e-3, double ptol = 1e-7)
            : curve(curve), surface(surface), tolerance(tolerance), ptol(ptol) {}

         std::vector<std::tuple<std::string, vec3, Vector<3>>> intersect() {
            intersections.clear();
            _curve_surface_intersect(curve, surface);
            return intersections;
        }

    private:
        const NURBSCurve &curve;
        const NURBSSurface &surface;
        double tolerance;
        double ptol;
        std::vector<std::tuple<std::string, vec3, Vector<3>>> intersections;

        bool _is_valid_point(vec3 point, const NURBSCurve& curve, const NURBSSurface& surface,
                              const Vector<3>& params) {
            vec3 curve_pt;

                curve.evaluate(params[0], curve_pt);
                vec3 surface_pt;
                surface.evaluate(params[1], params[2], surface_pt);
                return (curve_pt - surface_pt).length() <= tolerance;



            return true;
        }

        void _curve_surface_intersect(const NURBSCurve &curve, const NURBSSurface &surface) {
            // Check if we're already at a too small subdivision
            if (curve._interval[1] - curve._interval[0] < ptol ||
                surface._interval[0][1] - surface._interval[0][0] < ptol ||
                surface._interval[1][1] - surface._interval[1][0] < ptol) {
                return;
            }

            if (_no_new_intersections(curve, surface)) {
                return;
            }

            auto new_point = _find_new_intersection(curve, surface);
            auto [u0, u1] = surface._interval[0];
            auto [v0, v1] = surface._interval[1];
            double u = (u0 + u1) * 0.5;
            double v = (v0 + v1) * 0.5;

            if (std::abs(u - u0) < ptol || std::abs(u - u1) < ptol ||
                std::abs(v - v0) < ptol || std::abs(v - v1) < ptol) {
                return;
            }

            if (!new_point.has_value()) {
                auto [curve1, curve2] = curve.split(0.5 * (curve._interval[0] + curve._interval[1]), false); // normalize_knots=false
                auto [surface1, surface2, surface3, surface4] = surface.subdivide(u, v); // normalize_knots=false

                _curve_surface_intersect(curve1, surface1);
                _curve_surface_intersect(curve1, surface2);
                _curve_surface_intersect(curve1, surface3);
                _curve_surface_intersect(curve1, surface4);
                _curve_surface_intersect(curve2, surface1);
                _curve_surface_intersect(curve2, surface2);
                _curve_surface_intersect(curve2, surface3);
                _curve_surface_intersect(curve2, surface4);
            } else {
                auto [point, params] = new_point.value();
                double t = params[0], u = params[1], v = params[2];

                if (std::abs(u - u0) < ptol || std::abs(u - u1) < ptol ||
                    std::abs(v - v0) < ptol || std::abs(v - v1) < ptol) {
                    return;
                }

                if (_is_degenerate(params, curve, surface)) {
                    intersections.emplace_back("degenerate", point, params);
                } else {
                    intersections.emplace_back("transversal", point, params);
                }

                // Add spherical separability check


                // Split curve and surface at intersection point
                auto [curve1, curve2] = curve.split(t, false); // normalize_knots=false
                auto surfaces = surface.subdivide(u, v); // normalize_knots=false

                // Recursively check each combination
                for (const auto& s : surfaces) {
                    for (const auto& nurbs_curve : {curve1, curve2})
                    {
                        std::vector<vec3> curve_points = nurbs_curve.get_control_points3d();
                        std::vector<vec3> surface_points = s.control_points_flat3d();

                        if ( SphericalCentralProjectionTest2(point,curve_points, surface_points)) {
                            continue;

                        }else{
                            _curve_surface_intersect(nurbs_curve, s);
                        }
                    }


                }
            }
        }

        bool _no_new_intersections(const NURBSCurve &curve, const NURBSSurface &surface) {
            // Implement separability test
            auto e1=curve.get_control_points3d();
              auto e12=surface.control_points_flat3d();
            if (AABBTest(e1,e12))
            {
                return true;
            } else
            {
                std::vector<vec3> simplex;
                vec3 cpt;
                GJK(e1,e12,simplex,cpt,1e-5);

                return cpt.length()>tolerance;
            }

        }

        std::optional<std::pair<vec3, Vector<3>>> _find_new_intersection(const NURBSCurve &curve, const NURBSSurface &surface) {
            auto equation = [&](const Vector<3> &params) -> double {
                // First validate parameters are in range
                    vec3 curve_pt, surface_pt;
                    curve.evaluate(params[0], curve_pt);
                    surface.evaluate(params[1], params[2], surface_pt);
                    const auto d = (curve_pt - surface_pt);
                    return d.sqLength();

            };

            Vector<3> initial_guess = {(curve._interval[0] + curve._interval[1]) * 0.5,
                                                   (surface._interval[0][0] + surface._interval[0][1]) * 0.5,
                                                   (surface._interval[1][0] + surface._interval[1][1]) * 0.5};

            auto result = newtonsMethod<3>(equation, initial_guess,ptol, 15);

            if (_is_valid_parameter(result, curve._interval, surface._interval)) {
                vec3 curve_pt, surface_pt;
                curve.evaluate(result[0], curve_pt);
                surface.evaluate(result[1], result[2], surface_pt);
                if ((curve_pt - surface_pt).length() <= tolerance) {
                    return std::make_pair(curve_pt, result);
                }

            }
            return std::nullopt;
        }

        bool _is_valid_parameter(const Vector<3> &params, const std::array<double, 2> &curve_interval,
                                 const std::array<std::array<double, 2>, 2> &surface_interval) {
            return (params[0] >= curve_interval[0] && params[0] <= curve_interval[1]) &&
                   (params[1] >= surface_interval[0][0] && params[1] <= surface_interval[0][1]) &&
                   (params[2] >= surface_interval[1][0] && params[2] <= surface_interval[1][1]);
        }

        bool _is_degenerate(const Vector<3> &params, const NURBSCurve &curve, const NURBSSurface &surface) {
            vec3 curve_tangent;
            curve.derivative(params[0], curve_tangent);
            vec3 surface_normal;
            surface.normal(params[1], params[2], surface_normal);
            return std::abs(curve_tangent.dot(surface_normal)) < tolerance;
        }
    };
    inline double computePtol(const NURBSCurve &curve,double tolerance)
    {
        return tolerance*(curve._interval[1]-curve._interval[0])/curve.length() ;

    }
inline void csx( const NURBSCurve &curve, const NURBSSurface &surface,  std::vector<std::tuple<std::string, vec3, Vector<3>>> & result,double tol= 1e-3)
{
    double ptol_=computePtol(curve, tol);
    auto intersector=NURBSCurveSurfaceIntersector(curve,surface, tol,ptol_);
    result=std::move(intersector.intersect());


}

} // namespace cmmcore
#endif //CSX_H
