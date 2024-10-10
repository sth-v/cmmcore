//
// Created by Andrew Astakhov on 08.10.24.
//
#include "cmmcore/convexhull.h"
using namespace cmmcore;

int main() {
    std::vector<vec2> points = {{0.0, 3.0}, {2.0, 2.0}, {1.0, 1.0}, {2.0, 1.0}, {3.0, 0.0}, {0.0, 0.0}, {3.0, 3.0}};
    std::vector<vec2> hull ;
    hull=convex_hull2d(points);

    std::cout << "Convex Hull vec2s:" << std::endl;
    for (vec2 p : hull)
        std::cout << "(" << p.x << ", " << p.y << ")" << std::endl;

    return 0;
}