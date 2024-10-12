//
// Created by Andrew Astakhov on 08.10.24.
//

#ifndef __CONVEXHULL2D_H
#define __CONVEXHULL2D_H
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "cmmcore/vec.h"

namespace cmmcore{
// A utility function to find next to top in a stack
inline vec2 nextToTop(std::vector<vec2>& hull) {
    vec2 p = hull.back();
    hull.pop_back();
    vec2 res = hull.back();
    hull.push_back(p);
    return res;
}

// A utility function to swap two points
inline void swap(vec2& p1, vec2& p2) {
    vec2 temp = p1;
    p1 = p2;
    p2 = temp;
}

// A utility function to return square of distance between p1 and p2
inline double distSq(const vec2& p1, const vec2& p2) {
    return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
inline int orientation(const vec2& p, const vec2& q, const vec2& r) {
    double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    const double epsilon = 1e-9;
    if (std::abs(val) < epsilon) return 0;  // collinear
    return (val > 0) ? 1 : 2; // clock or counterclock wise
}

// Custom comparator for sorting points
class PointComparator {
private:
    vec2 p0;

public:
    PointComparator(const vec2& referencePoint) : p0(referencePoint) {}

    bool operator()(const vec2& p1, const vec2& p2) const {
        int o = orientation(p0, p1, p2);
        if (o == 0) {
            return distSq(p0, p1) < distSq(p0, p2);
        }
        return (o == 2);
    }
};

// Prints convex hull of a set of n points.
inline std::vector<vec2> convex_hull2d(std::vector<vec2>& points) {
    if (points.size() < 3) return points;  // Convex hull not possible with less than 3 points

    int n = points.size();
    
    // Find the bottommost point (and leftmost if there is a tie)
    int ymin = 0;
    for (int i = 1; i < n; i++) {
        if (points[i].y < points[ymin].y || 
           (std::abs(points[i].y - points[ymin].y) < 1e-9 && points[i].x < points[ymin].x))
            ymin = i;
    }

    // Place the bottom-most point at first position
    swap(points[0], points[ymin]);

    // Sort n-1 points with respect to the first point.
    // A point p1 comes before p2 in sorted output if p2
    // has larger polar angle (in counterclockwise direction)
    // than p1
    PointComparator comp(points[0]);
    std::sort(points.begin() + 1, points.end(), comp);

    // Keep only the farthest point when there are multiple points with the same angle
    int m = 1;
    for (int i = 1; i < n; i++) {
        while (i < n - 1 && orientation(points[0], points[i], points[i + 1]) == 0) {
            i++;
        }
        points[m] = points[i];
        m++;
    }

    // If modified array of points has less than 3 points, convex hull is not possible
    if (m < 3) return {};

    // Create an empty stack and push first three points to it.
    std::vector<vec2> hull;
    hull.push_back(points[0]);
    hull.push_back(points[1]);
    hull.push_back(points[2]);

    // Process remaining n-3 points
    for (int i = 3; i < m; i++) {
        while (hull.size() > 1 && orientation(nextToTop(hull), hull.back(), points[i]) != 2) {
            hull.pop_back();
        }
        hull.push_back(points[i]);
    }

    return hull;
}


}
#endif //__CONVEXHULL2D_H
