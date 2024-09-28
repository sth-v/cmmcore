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
    double u = 0.5;     // Parameter value
    int d = 2;          // Number of derivatives to compute
    bool is_periodic = false;

    std::vector<std::vector<double>> CK;  // To store the computed derivatives

    curve_derivs_alg1(n, p, U, P, u, d, CK, is_periodic);

    // Output the computed derivatives
    for (int k = 0; k <= d; ++k) {
        std::cout << "Derivative " << k << ": (";
        for (size_t l = 0; l < P[0].size(); ++l) {
            std::cout << CK[k][l];
            if (l < P[0].size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << ")\n";
    }

    return 0;
}