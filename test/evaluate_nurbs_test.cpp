//
// Created by Andrew Astakhov on 28.09.24.
//
#include <iostream>
#include "cmmcore/nurbs_utils.h"
using namespace cmmcore;
int main() {
    // Example knot vector
    std::vector<double> U = {0, 0, 0, 0.5, 1, 1, 1};
    // Example control points: Each control point has (x, y, z, weight)
    std::vector<std::vector<double>> P = {
        {0.0, 0.0, 0.0, 1.0},
        {0.5, 1.0, 0.0, 1.0},
        {1.0, 0.0, 0.0, 1.0},
        {1.5, -1.0, 0.0, 1.0},
        {2.0, 0.0, 0.0, 1.0}
    };

    int n = 4;          // Number of basis functions minus one
    int p = 2;          // Degree of the B-spline
    double u = 0.3;     // Parameter value
    bool is_periodic = false;

    std::vector<double> result;  // To store the computed point
    curve_point(n, p, U, P, u, result, is_periodic);

    std::cout << "Computed curve point: ("
              << result[0] << ", "
              << result[1] << ", "
              << result[2] << ")\n";

    return 0;
}