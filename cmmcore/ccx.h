/*
Copyright (c) 2024 Andrew Astakhov <aa@contextmachine.ru>. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef CCX_H
#define CCX_H
#ifdef CYTHON_ABI
#include "nurbs.h"
#else
#include "cmmcore/nurbs.h"
#endif

namespace cmmcore {
    /**
    * @brief Recursive intersection detection algorithm for two NURBS curves.
    *
    * This function attempts to find intersections between two given NURBS curves by
    * recursively splitting and testing their bounding boxes (AABB - Axis-Aligned Bounding Box).
    * If the bounding boxes of the curves intersect, the function refines the test by splitting
    * the curves into smaller segments, continuing the process until the tolerance (`tol`) is met.
    * The results of intersection points are stored as pairs of parameter values from the
    * parameter spaces of the two curves.
    *
    * @param c1 NURBSCurve&
    *     The first NURBS curve involved in the intersection detection.
    *
    * @param c2 NURBSCurve&
    *     The second NURBS curve involved in the intersection detection.
    *
    * @param tol double
    *     The tolerance value that determines the precision of the intersection.
    *     If the size of the bounding box is smaller than this tolerance, the midpoint
    *     of the curve intervals is considered a potential intersection.
    *
    * @param result std::vector<std::pair<double, double>>&
    *     A vector that stores the results of the intersection. Each pair contains two `double` values,
    *     representing the parameter values from the respective parameter spaces of `c1` and `c2`
    *     where the intersection is found.
    *
    * @details
    * The function follows these steps:
    * 1. Rebuild the Axis-Aligned Bounding Box (AABB) for both curves.
    * 2. If the bounding boxes intersect and are smaller than the given tolerance,
    *    it computes the midpoint of both curve intervals and checks if it is close
    *    enough to any previous results in the `result` vector.
    * 3. If the bounding boxes are larger than the tolerance, the curves are recursively
    *    split in half, and the intersection test is applied to the four combinations
    *    of the split segments (`a1, b1` and `a2, b2`).
    *
    * @remarks
    * :
    * - The function assumes that the `rebuildAABB`, `_update_interval`, `aabb`, and `interval` methods are
    *   defined for the `NURBSCurve` class.
    * - This algorithm uses an adaptive refinement approach, where curves are split recursively
    *   until the bounding boxes are small enough (i.e., below the given tolerance).
    *
    * @example
    * Here is an example of how to use the `ccx` function:
    *
    * .. code-block:: cpp
    *
    *     NURBSCurve curve1, curve2;
    *     double tolerance = 0.001;
    *     std::vector<std::pair<double, double>> intersections;
    *
    *     // Find intersections between curve1 and curve2
    *     ccx(curve1, curve2, tolerance, intersections);
    *
    *     // Output the results
    *     for (const auto& [t1, t2] : intersections) {
    *         std::cout << "Intersection at parameter: (" << t1 << ", " << t2 << ")" << std::endl;
    *     }
    */
    inline void ccx(NURBSCurve &c1, NURBSCurve &c2, double tol, std::vector<std::pair<double, double> > &result) {

        auto &bb1 = c1.aabb();
        auto &bb2 = c2.aabb();
        auto iv1 = c1.interval();
        auto iv2 = c2.interval();
        if (bb1.intersects(bb2)) {
            if ((bb1.max - bb1.min).length() < tol && (bb2.max - bb2.min).length() < tol) {
                size_t l = result.size();
                double _first = 0.5 * (iv1[0] + iv1[1]);
                double _second = 0.5 * (iv2[0] + iv2[1]);
                if (l > 0) {
                    if ((std::abs(result[l - 1].first - _first) <= tol) && (
                            std::abs(result[l - 1].second - _second) <= tol)) {
                        return;
                    }
                }
                result.emplace_back(_first, _second);
                return;
            }

            auto [a1,b1] = c1.split((iv1[1] - iv1[0]) * 0.5 + iv1[0], false);
            auto [a2,b2] = c2.split((iv2[1] - iv2[0]) * 0.5 + iv2[0], false);
            ccx(a1, a2, tol, result);
            ccx(b1, b2, tol, result);
            ccx(a1, b2, tol, result);
            ccx(b1, a2, tol, result);
        }
    }




}
#endif //CCX_H
