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


#ifndef GJK_H
#define GJK_H
#ifdef CYTHON_ABI
#include "vec.h"
#else
#include "cmmcore/vec.h"
#endif
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <csignal>
#include <regex>


namespace cmmcore{
inline int support_vector(const std::vector<vec3> &vertices, const vec3 &d) {
  double highest = -std::numeric_limits<double>::max();
  int supportIndex = -1;

  for (size_t i = 0; i < vertices.size(); ++i) {
    double dotValue = vertices[i].dot( d);
    if (dotValue > highest) {
      highest = dotValue;
      supportIndex = static_cast<int>(i);
    }
  }

  return supportIndex;
}

inline bool handle_simplex2(std::vector<vec3>& simplex, vec3& d, double tol) {
    switch (simplex.size()) {
        case 1: {
            // Point simplex
            vec3 ao = -simplex[0];
            if (ao.length() < tol) {
                // Origin is at the point
                return true;
            }
            d = ao;
            break;
        }
        case 2: {
            const auto& a = simplex[1];
            const auto& b = simplex[0];
            vec3 ab = b - a;
            vec3 ao = -a;

            if (ab.dot(ao) > tol) {
                vec3 temp = ab.cross(ao).cross(ab);
                if (temp.length() < tol) {
                    // Origin is on the line segment
                    return true;
                }
                d = temp;
            } else {
                simplex.erase(simplex.begin());
                d = ao;
            }
            break;
        }
        case 3: {
            const auto& a = simplex[2];
            const auto& b = simplex[1];
            const auto& c = simplex[0];
            vec3 ab = b - a;
            vec3 ac = c - a;
            vec3 ao = -a;
            vec3 abc = ab.cross( ac);

            vec3 ab_abc = abc.cross(ab);
            vec3 abc_ac = ac.cross(abc);

            if (ab_abc.dot(ao) > tol) {
                simplex.erase(simplex.begin());
                d = ab_abc;
            } else if (abc_ac.dot(ao) > tol) {
                    simplex.erase(simplex.begin() + 1);
                d = abc_ac;
                } else {
                if (std::abs(abc.dot(ao)) < tol) {
                    // Origin is on the triangle
                    return true;
                }
                    if (abc.dot( ao) > tol) {
                        d = abc;
                    } else {
                        std::swap(simplex[0], simplex[1]);
                    d = -abc;
                }
            }
            break;
        }
        case 4: {
          const auto &a = simplex[3];
          const auto &b = simplex[2];
          const auto &c = simplex[1];
          const auto &dPoint = simplex[0];

          vec3 ab = b - a;
          vec3 ac = c - a;
          vec3 ad = dPoint - a;
          vec3 ao = {-a[0], -a[1], -a[2]};

          vec3 abc = ab.cross( ac);
          vec3 acd = ac.cross( ad);
          vec3 adb = ad.cross( ab);

          if (abc.dot( ao) > 0) {
            simplex.erase(simplex.begin());
            return handle_simplex2(simplex, d, tol);
          } if (acd.dot( ao) > 0) {
            simplex.erase(simplex.begin() + 1);
            return handle_simplex2(simplex, d, tol);
          }  if (acd.dot( ao) > 0) {
            simplex.erase(simplex.begin() + 2);
            return handle_simplex2(simplex, d, tol);
          } 
            return true;

        }
        default:
            throw std::invalid_argument("Invalid simplex size "+simplex.size());
        }
    return false;
}


inline bool gjk_collision_detection( const std::vector<vec3>& vertices1, const std::vector<vec3>& vertices2, const double tol, size_t maxIter = 0) {
    if (vertices1.empty() || vertices2.empty()) {
        // Handle empty sets
        return false;
    }

    if (maxIter == 0) {

        maxIter = std::min<size_t>(vertices1.size() * vertices2.size(), 1000);
    }

        const size_t rows = vertices1.size();
        const size_t cols = vertices2.size();
        std::vector<bool> visited(rows * cols, false);
        std::vector<vec3> cache(rows * cols);

        std::vector<vec3> simplex;
        simplex.reserve(4);  // Pre-allocate space for efficiency

        auto support = [&](const vec3& d) -> vec3 {
            int i = support_vector(vertices1, d);
            int j = support_vector(vertices2, {-d[0], -d[1], -d[2]});
            return vertices1[i] - vertices2[j];
        };

        vec3 d = {1, 0, 0};
        auto newPoint = support(d);
        size_t index = static_cast<size_t>(support_vector(vertices1, d)) * cols + static_cast<size_t>(support_vector(vertices2, {-d[0], -d[1], -d[2]}));
        visited[index] = true;
        cache[index] = newPoint;
        simplex.push_back(newPoint);
        d = {-newPoint[0], -newPoint[1], -newPoint[2]};

        for (size_t iter = 0; iter < maxIter; ++iter) {
            newPoint = support(d);


            if (newPoint.dot( d) < 0) {
                printf("newPoint.dot( d) {%f} < 0\n", newPoint.dot( d) );
                return false;
            }

            simplex.push_back(newPoint);

            if (handle_simplex2(simplex, d, tol)) {
                return true;
            }
        }
    printf("WARNING!!!! GJK failed to converge after %zu iterations\n", maxIter);

        return false;
    }

}


#endif //GJK_H
