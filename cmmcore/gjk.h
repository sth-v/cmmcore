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


#ifndef CMMCORE_GJK_H
#define CMMCORE_GJK_H



#ifdef CYTHON_ABI
#include "vec.h"
#include "closest_point.h"
#else
#include "cmmcore/vec.h"
#include "cmmcore/closest_point.h"

#endif
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <array>
#include <map>

namespace cmmcore{
static const size_t CMMCORE_SIMPLEX_FACES[4][3]= {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};

inline bool handleSimplex(std::vector<vec3>& simplex, vec3& direction, double tol) {
    switch (simplex.size()) {
        case 2: {
            // Line segment case
            const vec3& a = simplex[1];
            const vec3& b = simplex[0];
            vec3 ab = b - a;
            vec3 ao = -a;

            if (ab.dot(ao) > tol) {
                // Direction is perpendicular to ab towards the origin
                direction = ab.cross(ao).cross(ab);
            } else {
                // Remove point not contributing to the simplex
                simplex.erase(simplex.begin());
                direction = ao;
            }
            break;
        }
        case 3: {
            // Triangle case
            const vec3& a = simplex[2];
            const vec3& b = simplex[1];
            const vec3& c = simplex[0];
            vec3 ab = b - a;
            vec3 ac = c - a;
            vec3 ao = -a;
            vec3 abc = ab.cross(ac);

            vec3 abPerp = abc.cross(ac);
            vec3 acPerp = ab.cross(abc);

            if (abPerp.dot(ao) > tol) {
                // Region outside AB
                simplex.erase(simplex.begin());
                direction = abPerp;
            } else if (acPerp.dot(ao) > tol) {
                // Region outside AC
                simplex.erase(simplex.begin() + 1);
                direction = acPerp;
            } else {
                if (abc.dot(ao) > tol) {
                    // Origin is below the triangle
                    direction = abc;
                } else {
                    // Origin is above the triangle
                    std::swap(simplex[0], simplex[1]);
                    direction = -abc;
                }
            }
            break;
        }
        case 4: {
            // Tetrahedron case
            const vec3& a = simplex[3];
            const vec3& b = simplex[2];
            const vec3& c = simplex[1];
            const vec3& d = simplex[0];

            vec3 ab = b - a;
            vec3 ac = c - a;
            vec3 ad = d - a;
            vec3 ao = -a;

            vec3 abc = ab.cross(ac);
            vec3 acd = ac.cross(ad);
            vec3 adb = ad.cross(ab);

            if (abc.dot(ao) > 0) {
                // Remove point D
                simplex.erase(simplex.begin());
                direction = abc;
            } else if (acd.dot(ao) > 0) {
                // Remove point B
                simplex.erase(simplex.begin() + 2);
                direction = acd;
            } else if (adb.dot(ao) > 0) {
                // Remove point C
                simplex.erase(simplex.begin() + 1);
                direction = adb;
            } else {
                // Origin is inside the tetrahedron
                return true;
            }
            break;
        }
        default:
            throw std::invalid_argument("Invalid simplex size " + std::to_string(simplex.size()));
    }
    return false;
}

// Overload operators for vec3 to use in maps and comparisons
inline bool operator<(const vec3& lhs, const vec3& rhs) {
    if (lhs[0] < rhs[0]) return true;
    if (lhs[0] == rhs[0] && lhs[1] < rhs[1]) return true;
    if (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] < rhs[2]) return true;
    return false;
}



// Function to compute the closest point to the origin on the current simplex
inline vec3 closestPointOnSimplex(const std::vector<vec3>& simplex) {
#ifdef CMMCORE_DEBUG
    printf("GJK: closest point at  %d-simplex\n", simplex.size());
#endif
    if (simplex.size() == 1) {
        // Only one point in the simplex
        return simplex[0];
    } else if (simplex.size() == 2) {
        // Line segment case
        return closestPointOnLine(simplex[1], simplex[0]);
    } else if (simplex.size() == 3) {
        // Triangle case
        return closestPointOnTriangle(simplex[2], simplex[1], simplex[0]);
    } else {

        // Tetrahedron case (should not happen in non-colliding case)
        vec3 cpt;
        double mindist=std::numeric_limits<double>::max();
        for (size_t i = 0; i < simplex.size(); ++i) {

            auto pt=closestPointOnTriangle(simplex[CMMCORE_SIMPLEX_FACES[i][0]], simplex[ CMMCORE_SIMPLEX_FACES[i][1]], simplex[CMMCORE_SIMPLEX_FACES[i][2]]);
            double dst=pt.sqLength();
            if (dst<mindist) {
              cpt=pt;
            }
            
        }
        return cpt;
    }
}

    inline vec3 closestPoint(const std::vector<vec3>& simplex){
      return closestPointOnSimplex(simplex);
    }

inline int supportVector(const std::vector<vec3> &vertices, const vec3 &d) {
    double highest = -std::numeric_limits<double>::max();
    int supportIndex = -1;

    for (size_t i = 0; i < vertices.size(); ++i) {
      double dotValue = vertices[i].dot(d);
      if (dotValue > highest) {
        highest = dotValue;
        supportIndex = static_cast<int>(i);
      }
    }
    return supportIndex;

}
// Support function for Minkowski Difference
inline vec3 support(const std::vector<vec3>& verticesA, const std::vector<vec3>& verticesB, const vec3& d) {
    int indexA = supportVector(verticesA, d);
    int indexB = supportVector(verticesB, -d);
    return verticesA[indexA] - verticesB[indexB];
}

inline bool GJK(const std::vector<vec3>&verticesA,
    const std::vector<vec3>& verticesB,
    std::vector<vec3>& simplex,
    vec3& closestPointToOrigin,
    const double tolerance = std::numeric_limits<double>::epsilon(),
    const size_t maxIterations = 25) {

 // Swap A and B to make the algorithm easier to read




    vec3 direction = verticesB[0] - verticesA[0]; // Initial direction
    if (direction.sqLength() == 0) {
        direction = vec3(1, 0, 0); // Arbitrary direction
    }

    // Initial support point
    vec3 point = support(verticesA, verticesB, direction);
    simplex.push_back(point);

    // New direction towards the origin
    direction = -point;


    int iterations = 0;


    while (iterations++ < maxIterations) {
        point = support(verticesA, verticesB, direction);
        if (point.dot(direction) < 0) {
            // No collision
            // Compute the closest point to the origin on the simplex
            closestPointToOrigin = closestPointOnSimplex(simplex);
            #ifdef CMMCORE_DEBUG
            printf("GJK iterations: %d\n", iterations);
            #endif

            return false;
        }

        simplex.push_back(point);

        if (handleSimplex(simplex, direction, tolerance)) {
            // Collision detected
            #ifdef CMMCORE_DEBUG
            printf("GJK iterations: %d\n", iterations);
            #endif

            return true;
        }
    }
    closestPointToOrigin = closestPointOnSimplex(simplex);
    // Should not reach here
#ifdef CMMCORE_DEBUG
    printf("GJK iterations: %d\n", iterations);
    printf(("["+format_vec3vec( verticesA)+","+format_vec3vec( verticesB)+"]\n\n" ).c_str() );
#endif
    closestPointToOrigin = {0.0, 0.0, 0.0};
    printf("WARNING GJK did not converge");
    return true;
}
}
#endif //CMMCORE_GJK_H
