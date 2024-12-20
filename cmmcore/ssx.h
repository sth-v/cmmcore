//
// Created by Andrew Astakhov on 09.10.24.
//

#ifndef SSX_H
#define SSX_H

#include "cmmcore/nurbs.h"
#include "cmmcore/gauss_map.h"
#include "cmmcore/boundary_intersection.h"
#include "cmmcore/point_inversion.h"

namespace cmmcore
{
    namespace ssx_helpers
    {
        inline bool improve_uv(const vec3& du, const vec3& dv, const vec3& delta, vec2& result)
        {
            Matrix2x2 mat[3]{
                {
                    {du.x, dv.x},
                    {du.y, dv.y}
                },
                {
                    {du.x, dv.x},
                    {du.z, dv.z}
                },
                {

                    {du.y, dv.y},
                    {du.z, dv.z}
                }
            };
            double dett = std::numeric_limits<double>::min();
            double temp;
            size_t n = -1;

            double y[3][2] = {
                {delta.x, delta.y},
                {delta.x, delta.z},
                {delta.y, delta.z}
            };
            for (size_t i = 0; i < 3; ++i)
            {
                temp = mat[i].det();
                if (temp > dett)
                {
                    dett = temp;
                    n = i;
                }
            }
            if (n == -1)
            {
                return false;
            }
            if (std::abs(dett) <= CMMCORE_NUMERIC_EPSILON_DOUBLE)
            {
                return false;
            }
            result.x = (y[n][0] * mat[n][1][1] - mat[n][0][1] * y[n][1]) / dett;
            result.y = (mat[n][0][0] * y[n][1] - y[n][0] * mat[n][1][0]) / dett;
            return true;
        }
    }

    inline void improve_uv_robust(
        const NURBSSurface& surf,
        const vec2& uv_old,
        const vec3& du,
        const vec3& dv,
        const vec3& xyz_old,
        const vec3& xyz_new,
        vec2& uv_new,
        double tol = 1e-7, size_t max_iter = 100)
    {
        auto delta = xyz_new - xyz_old;
        if (ssx_helpers::improve_uv(du, dv, delta, uv_new))
        {
            uv_new += uv_old;
        }

        uv_new.set(point_inversion_surface(surf, xyz_new, uv_old.x, uv_old.y, tol, tol, max_iter));
    }

    static int counter = 0;
    using UnorderedSetPairDouble = std::unordered_set<
        std::pair<double, double>, _hashimpl::PairDoubleHash, _hashimpl::PairDoubleEqual>;


    struct PointIntersection
    {
        vec2 uv1;
        vec2 uv2;
        vec3 xyz;
        PointIntersection() = default;

        PointIntersection(const vec2& _uv1, const vec2& _uv2, const vec3& _xyz): uv1(_uv1), uv2(_uv2), xyz(_xyz)
        {
        }
    };

    struct PatchIntersection
    {
        //NURBSSurface s1, s2;
        IntersectionPoint start, end;
        PatchIntersection() = default;

        PatchIntersection(const IntersectionPoint& a,
                          const IntersectionPoint& b): start(a), end(b)
        {
        }
    };

    struct Intersection
    {
        std::vector<PatchIntersection> patches{};
        std::vector<PointIntersection> points{};
        Intersection() = default;
    };

    inline void processSimpleSurfaceIntersection(const NURBSSurface& surface1, const NURBSSurface& surface2,
                                                 const double tol, Intersection& intersection)
    {
        std::vector<IntersectionPoint> ixs;
        //printf("OOO\n");
        find_boundary_intersections(surface1, surface2, ixs, tol);
        auto n = ixs.size();
        if (n >= 2)
        {
            intersection.patches.emplace_back(ixs[0], ixs[n - 1]);
        }

        else
        {
#ifdef CMMCORE_DEBUG
            printf("\nError in processSimpleSurfaceIntersection %lu\n", ixs.size());
            printf("[[");
            for (auto& p :surface1.control_points_flat3d())
            {

                printf("[%f,%f,%f],", p.x,p.y,p.z);
            }
            printf("],\n");
            printf("[");
            for (auto& p :surface2.control_points_flat3d())
            {

                printf("[%f,%f,%f],", p.x,p.y,p.z);
            }
            printf("]]\n---\n");


#endif
        }
    }

    //
    inline void detectIntersections(NURBSSurface& surface1, NURBSSurface& surface2, NURBSSurface& gaussMap1,
                                    NURBSSurface& gaussMap2, const double tol, Intersection& intersection,
                                    int depth = 0)
    {
        counter++;

        const auto& bb1 = surface1.bbox();
        const auto& bb2 = surface2.bbox();

        if (!(bb1.intersects(bb2)))
        {
#ifdef CMMCORE_DEBUG
            printf("no intersection bb . depth: %d\n", depth);
            printf("[[%f,%f,%f],[%f,%f,%f]]\n", bb1.min.x, bb1.min.y, bb1.min.z,bb1.max.x, bb1.max.y, bb1.max.z);
            printf("[[%f,%f,%f],[%f,%f,%f]]\n", bb2.min.x, bb2.min.y, bb2.min.z,bb2.max.x, bb2.max.y, bb2.max.z);
#endif

            return;
        }
        auto bb1d = bb1.max - bb1.min;
        auto bb2d = bb2.max - bb2.min;

        if (bb1d.x <= tol && bb1d.y <= tol && bb1d.z <= tol && bb2d.x <= tol && bb2d.y <= tol && bb2d.z <= tol)
        {
#ifdef CMMCORE_DEBUG
            printf("both surfaces are too small. depth: %d\n", depth);
#endif

            intersection.points.push_back(
                {
                    {
                        0.5 * (surface1._interval[0][1] + surface1._interval[0][0]),
                        0.5 * (surface1._interval[1][1] + surface1._interval[1][0])
                    },
                    {
                        0.5 * (surface2._interval[0][1] + surface2._interval[0][0]),
                        0.5 * (surface2._interval[1][1] + surface2._interval[1][0])
                    },
                    (bb1.min + bb1.max + bb2.min + bb2.max) * 0.25

                }
            );
            return;
        }

        AABB aabb1;
        bb1.intersection(bb2, aabb1);
        auto _dims = (aabb1.max - aabb1.min);

        if ((aabb1.volume() < tol || std::abs(_dims.x) <= tol) || (std::abs(_dims.y) <= tol) || (std::abs(_dims.z) <=
            tol))
        {
#ifdef CMMCORE_DEBUG
            printf("a surface is planar. depth: %d\n", depth);
            printf("[[%f,%f,%f],[%f,%f,%f]],[[%f,%f,%f],[%f,%f,%f]]\n", bb1.min.x, bb1.min.y, bb1.min.z,bb1.max.x, bb1.max.y, bb1.max.z,bb2.min.x, bb2.min.y,bb2.min.z,bb2.max.x,bb2.max.y,bb2.max.z);
#endif
            return;
        }
        if (SAT3D(surface1.control_points_flat3d(), surface2.control_points_flat3d(), tol))
        {
#ifdef CMMCORE_DEBUG
            printf("SAT3D. depth: %d\n", depth);
#endif

            return;
        };

        if (gaussMapSeparability(gaussMap1, gaussMap2))
        {
            auto vv = surface1.control_points_flat3d();
            auto vv1 = surface2.control_points_flat3d();
#ifdef CMMCORE_DEBUG
            printf(("["+format_vec3vec( vv)+"],["+format_vec3vec(  vv1)+"]\n\n").c_str(), depth);
            printf("gauss maps are separable. depth: %d\n", depth);
#endif
            processSimpleSurfaceIntersection(surface1, surface2, tol, intersection);

            return;
        }

        std::array<NURBSSurface, 4> g1; // g11,g12,g13,g14
        std::array<NURBSSurface, 4> g2; //g21,g22,g23,g24;
        std::array<NURBSSurface, 4> s1; //s11,s12,s13,s14;
        std::array<NURBSSurface, 4> s2; //s21,s22,s23,s24;
        surface1.subdivide(s1[0], s1[1], s1[2], s1[3]);
        surface2.subdivide(s2[0], s2[1], s2[2], s2[3]);
        gaussMap1.subdivide(g1[0], g1[1], g1[2], g1[3]);
        gaussMap2.subdivide(g2[0], g2[1], g2[2], g2[3]);


        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                detectIntersections(s1[i], s2[j], g1[i], g2[j], tol, intersection, depth + 1);
            }
        }
    }

    // Finds at least one point on each intersection branch (usually much more).
    // Used primarily for obtaining the start points of the marching method
    inline void findStartIntersectionPoints(const NURBSSurface& ns1, const NURBSSurface& ns2, const double tol,
                                            Intersection& intersection)
    {
        std::vector<NURBSSurface> surfs;
        std::vector<NURBSSurface> surfs2;


        ns1.decompose(surfs);
        ns2.decompose( surfs2);

        std::vector<cmmcore::NURBSSurface> g1;
        std::vector<cmmcore::NURBSSurface> g2;


        for (size_t i = 0; i < surfs.size(); i++)
        {
            //surfs[i].bbox();
            g1.push_back(cmmcore::gaussMap(surfs[i]));
            //g1[i].bbox();
        }
        for (size_t i = 0; i < surfs2.size(); i++)
        {
            ///surfs2[i].bbox();
            g2.push_back(cmmcore::gaussMap(surfs2[i]));
            //g2[i].bbox();
        }

        for (size_t i = 0; i < surfs.size(); i++)
        {
            for (size_t j = 0; j < surfs2.size(); j++)
            {
                detectIntersections(surfs[i],
                                             surfs2[j],
                                             g1[i],
                                             g2[j],
                                             tol,
                                             intersection);
            }
        }
    }
}
#endif //SSX_H
