//
// Created by Andrew Astakhov on 28.09.24.
//
#include <iostream>
#include "cmmcore/nurbs_utils.h"
using namespace cmmcore;


int main() {
    // Example parameters (to be filled with actual data)
    int p = 3;  // Degree of B-spline curve
    std::vector<double> U = {0, 0, 0, 0.5, 1, 1, 1};  // Knot vector
    std::vector<vec4> P = {
        {0.0, 0.0, 0.0, 1.0},
        {0.5, 1.0, 0.0, 1.0},
        {1.0, 0.0, 0.0, 1.0},
        {1.5, -1.0, 0.0, 1.0},
        {2.0, 0.0, 0.0, 1.0}
    };  // Control points

    int d = 2;    // Number of derivatives to compute
    int r1 = 0;   // Start index
    int r2 = 3;   // End index

    std::vector<std::vector<vec4>> PK;  // Output derivative control points

    // Compute derivative control points
    curve_deriv_cpts(p, U, P, d, r1, r2, PK);

    // Output the derivative control points
    for (int k = 0; k <= d; ++k) {
        std::cout << "Derivative order " << k << ":\n";
        for (const auto& point : PK[k]) {

                std::cout <<"["<< point.x << ", " << point.y << ", " << point.z << "]," << std::endl;

            std::cout << "\n";
        }
        std::cout << "\n";
    }

    return 0;
}