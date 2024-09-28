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
#include "cmmcore/nurbs.h"

namespace cmmcore {
    inline void ccx(NURBSCurve &c1, NURBSCurve &c2, double tol, std::vector<std::pair<double, double> > &result) {
        c1.rebuildAABB();
        c2.rebuildAABB();
        c1._update_interval();
        c2._update_interval();
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
