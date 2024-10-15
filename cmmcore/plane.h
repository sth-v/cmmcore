//
// Created by Andrew Astakhov on 11.10.24.
//

#ifndef PLANE_H
#define PLANE_H
#include "cmmcore/vec.h"
#include <tuple>

namespace cmmcore
{
    inline void plane_eq_from_pts(const vec3& pt0, const vec3& pt1, const vec3& pt2, vec4& equation)
    {
        const float x3 = pt0.y * pt1.z;
        const float x4 = pt1.y * pt2.z;
        const float x5 = pt2.y * pt0.z;
        equation.x = x3 + x4 + x5 - pt0.y * pt2.z - pt1.y * pt0.z - pt2.y * pt1.z;
        equation.y = -pt0.x * pt1.z + pt0.x * pt2.z + pt1.x * pt0.z - pt1.x * pt2.z - pt2.x * pt0.z + pt2.x * pt1.z;
        equation.z = pt0.x * pt1.y - pt0.x * pt2.y - pt1.x * pt0.y + pt1.x * pt2.y + pt2.x * pt0.y - pt2.x * pt1.y;
        equation.w = -pt0.x * x4 + pt0.x * pt2.y * pt1.z - pt1.x * x5 + pt1.x * pt0.y * pt2.z - pt2.x * x3 + pt2.x * pt1
            .y * pt0.z;
    }
    inline double planeSignedDistance(const vec3& pt, const vec4& plane)
    {
        return plane.x * pt.x + plane.y * pt.y + plane.z * pt.z + plane.w;
    };

    class plane
    {
        double denominator = 1;

    public:
        vec4 equation{0, 0, 1, 0};
        vec3 normal{0, 0, 1};
        vec3 origin{0, 0, 0};
        plane() = default;

        plane(double a, double b, double c, double d = 0) : equation(a, b, c, d)
        {
            computeDenominator();

            normal = vec3(equation.x / denominator, equation.y / denominator, equation.z / denominator);
            origin.set(normal * equation.w);
        }

        plane(const vec3& _normal) : denominator(_normal.length()), equation(_normal.x, _normal.y, _normal.z, 0),
                                     normal(_normal / denominator)
        {
        }

        plane(const vec3& _origin, const vec3& _normal) : denominator(_normal.length()),
                                                          equation(_normal.x, _normal.y, _normal.z,
                                                                   -(_normal.x * _origin.x, _normal.y * _origin.y,
                                                                       _normal.z * _origin.z)),
                                                          normal(_normal / denominator), origin(_origin)
        {
        }

        plane(const std::tuple<vec3, vec3, vec3>& _pts)
        {
            plane_eq_from_pts(std::get<0>(_pts), std::get<1>(_pts), std::get<2>(_pts), equation);
            computeDenominator();
            normal = vec3(equation.x / denominator, equation.y / denominator, equation.z / denominator);
            origin.set(normal * equation.w);
        };

        void computeDenominator()
        {
            denominator = std::sqrt(equation.x * equation.x + equation.y * equation.y + equation.z * equation.z);
        }

        double signedDistance(const vec3& pt) const
        {
            return (pt.x * equation.x + pt.y * equation.y + pt.z * equation.z + equation.w) / denominator;
        };

        void signedDistance(const std::vector<vec3>& pts, std::vector<double>& result) const
        {
            result.resize(pts.size());
            for (size_t i = 0; i < pts.size(); ++i)
            {
                result[i] = signedDistance(pts[i]);
            }
        }

        bool planeRayIntersection(const vec3& start, const vec3& direction, vec3& pt) const
        {
            constexpr double epsilon = 1e-8;
            double N = -(equation.x * start.x + equation.y * start.y + equation.z * start.z + equation.w);
            double D = equation.x * direction.x + equation.y * direction.y + equation.z * direction.z;

            if (std::abs(D) < epsilon)
            {
                // The direction is parallel to the plane
                if (std::abs(N) < epsilon)
                {
                    // The starting point is on the plane
                    pt = start;
                }
                else
                {
                    // No intersection
                    return false;
                }
            }
            else
            {
                double t = N / D;
                pt = start + direction * t;
            }
            return true;
        }

        bool planeLineIntersection(const vec3& start, const vec3& end, vec3& pt) const

        {
            constexpr double epsilon = std::numeric_limits<double>::epsilon();
            vec3 direction = end - start;
            double N = -(equation.x * start.x + equation.y * start.y + equation.z * start.z + equation.w);
            double D = equation.x * direction.x + equation.y * direction.y + equation.z * direction.z;

            if (std::abs(D) < epsilon)
            {
                // The segment is parallel to the plane
                if (std::abs(N) < epsilon)
                {
                    // The entire segment lies on the plane
                    pt = start;
                    return true;
                }
                else
                {
                    // No intersection
                    return false;
                }
            }
            else
            {
                double t = N / D;
                if (t < 0.0 || t > 1.0)
                {
                    // The intersection point is outside the segment
                    return false;
                }
                // The intersection point is within the segment
                pt = start + direction * t;
                return true;
            }
        }

        bool planePlaneIntersection(const plane& other, vec3& pt1, vec3& pt2) const
        {
            constexpr double epsilon = std::numeric_limits<double>::epsilon();
            vec3 dir = normal.cross(other.normal);

            if (dir.sqLength() < epsilon)
            {
                // The planes are parallel

                return false;
            }

            double a1 = equation.x;
            double b1 = equation.y;
            double c1 = equation.z;
            double d1 = equation.w;

            double a2 = other.equation.x;
            double b2 = other.equation.y;
            double c2 = other.equation.z;
            double d2 = other.equation.w;

            // Choose the variable to eliminate
            int maxIndex = 0;
            double maxAbs = std::abs(a1);

            if (std::abs(b1) > maxAbs)
            {
                maxIndex = 1;
                maxAbs = std::abs(b1);
            }
            if (std::abs(c1) > maxAbs)
            {
                maxIndex = 2;
            }

            vec3 point;

            if (maxIndex == 0)
            {
                // Eliminate x
                double det = b1 * c2 - b2 * c1;
                if (std::abs(det) < epsilon)
                {
                    return false;
                }
                double y = (c2 * -d1 - c1 * -d2) / det;
                double z = (b1 * -d2 - b2 * -d1) / det;
                point = vec3(0, y, z);
            }
            else if (maxIndex == 1)
            {
                // Eliminate y
                double det = a1 * c2 - a2 * c1;
                if (std::abs(det) < epsilon)
                {
                    return false;
                }
                double x = (c2 * -d1 - c1 * -d2) / det;
                double z = (a1 * -d2 - a2 * -d1) / det;
                point = vec3(x, 0, z);
            }
            else
            {
                // Eliminate z
                double det = a1 * b2 - a2 * b1;
                if (std::abs(det) < epsilon)
                {
                    return false;
                }
                double x = (b2 * -d1 - b1 * -d2) / det;
                double y = (a1 * -d2 - a2 * -d1) / det;
                point = vec3(x, y, 0);
            }

            pt1 = point;
            pt2 = point + dir;
            return true;
        }

        bool planePlanePlaneIntersection(const plane& other1, const plane& other2, vec3& pt) const
        {
            constexpr double epsilon = 1e-8;
            double a1 = equation.x;
            double b1 = equation.y;
            double c1 = equation.z;
            double d1 = equation.w;

            double a2 = other1.equation.x;
            double b2 = other1.equation.y;
            double c2 = other1.equation.z;
            double d2 = other1.equation.w;

            double a3 = other2.equation.x;
            double b3 = other2.equation.y;
            double c3 = other2.equation.z;
            double d3 = other2.equation.w;

            double Det = determinant3x3(a1, b1, c1,
                                        a2, b2, c2,
                                        a3, b3, c3);

            if (std::abs(Det) < epsilon)
            {
                return false;
            }

            double Dx = determinant3x3(-d1, b1, c1,
                                       -d2, b2, c2,
                                       -d3, b3, c3);

            double Dy = determinant3x3(a1, -d1, c1,
                                       a2, -d2, c2,
                                       a3, -d3, c3);

            double Dz = determinant3x3(a1, b1, -d1,
                                       a2, b2, -d2,
                                       a3, b3, -d3);

            double x = Dx / Det;
            double y = Dy / Det;
            double z = Dz / Det;

            pt = vec3(x, y, z);
            return true;
        }

    private:
        static double determinant3x3(double m00, double m01, double m02,
                                     double m10, double m11, double m12,
                                     double m20, double m21, double m22)
        {
            return m00 * (m11 * m22 - m12 * m21)
                - m01 * (m10 * m22 - m12 * m20)
                + m02 * (m10 * m21 - m11 * m20);
        }
    };

    inline bool isPtsCoplanar(const std::vector<vec3>& pts, const double eps = std::numeric_limits<double>::epsilon())
    {
        if (pts.size() <= 3)
        {
            return true;
        }
        plane p = {{pts[0], pts[1], pts[2]}};


        for (size_t i = 4; i < pts.size(); i++)
        {
            if (std::abs(p.signedDistance(pts[i])) >= eps)
            {
                return false;
            }
        }
        return true;
    }

    class cplane : plane
    {
    public:
        vec3 xaxis{1, 0, 0};
        vec3 yaxis{0, 1, 0};

        cplane() = default;

        cplane(const vec3& _origin, const vec3& xaxis, const vec3& yaxis) : plane(_origin, xaxis.cross(yaxis)),
                                                                            xaxis(xaxis), yaxis(yaxis)
        {
        }

        cplane(const std::tuple<vec3, vec3, vec3>& _pts): plane(_pts)
        {
            xaxis = std::get<0>(_pts) - origin;
            xaxis.unitize();
            yaxis = normal.cross(xaxis);
            yaxis.unitize();
        };
        double signedDistance(const vec3& pt) const { return plane::signedDistance(pt - origin); };

        void evaluate(double u, double v, vec3& result) const
        {
            result.set(origin);
            result += (xaxis * u);
            result += (yaxis * v);
        }

        void evaluateInv(const vec3& xyz, vec3& uvh) const
        {
            auto temp = xyz - origin;
            uvh.x = temp.dot(xaxis);
            uvh.y = temp.dot(yaxis);
            uvh.z = temp.dot(normal);
        }

        void evaluateInv(const std::vector<vec3>& xyz, std::vector<vec3>& uvh) const
        {
            uvh.resize(xyz.size());
            vec3 temp;
            for (size_t i = 0; i < xyz.size(); ++i)
            {
                xyz[i].sub(origin, temp);
                uvh[i].x = temp.dot(xaxis);
                uvh[i].y = temp.dot(yaxis);
                uvh[i].z = temp.dot(normal);
            }
        }

    };

};
#endif //PLANE_H
