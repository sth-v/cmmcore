//
// Created by Andrew Astakhov on 08.10.24.
//

#include "cmmcore/monomial.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace cmmcore;

void printMatrix(const Matrix& mat, const std::string& name) {
    std::cout << name << ":" << std::endl;
    for (size_t i = 0; i < mat.getRows(); ++i) {
        for (size_t j = 0; j < mat.getCols(); ++j) {
            std::cout << std::setw(10) << std::setprecision(4) << mat(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void printTensor3D(const Tensor3D& tensor, const std::string& name) {
    std::cout << name << ":" << std::endl;
    for (size_t d = 0; d < tensor.size(); ++d) {
        std::cout << "Dimension " << d << ":" << std::endl;
        printMatrix(tensor[d], "");
    }
}

int main() {
    // Create a simple 3x3 Bézier surface patch
    Tensor3D control_points(3, Matrix(3, 3));

    // Fill control points with some example values
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            control_points[i](j, 0) = std::sin(i * j * M_PI / 4);  // x-coordinate
            control_points[i](j, 1) = std::cos(i * j * M_PI / 4);  // y-coordinate
            control_points[i](j, 2) = (i + j) / 4.0;  // z-coordinate
        }
    }

    std::cout << "Original Bézier control points:" << std::endl;
    printTensor3D(control_points, "Control Points");

    // Convert Bézier to monomial form
    Tensor3D monomial_coeffs = bezier_to_monomial(control_points);

    std::cout << "Monomial coefficients:" << std::endl;
    printTensor3D(monomial_coeffs, "Monomial Coefficients");

    // Convert monomial back to Bézier form
    Tensor3D reconstructed_control_points = monomial_to_bezier(monomial_coeffs);

    std::cout << "Reconstructed Bézier control points:" << std::endl;
    printTensor3D(reconstructed_control_points, "Reconstructed Control Points");

    // Calculate and print the maximum difference between original and reconstructed control points
    double max_diff = 0.0;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            for (size_t d = 0; d < 3; ++d) {
                double diff = std::abs(control_points[i](j, d) - reconstructed_control_points[i](j, d));
                max_diff = std::max(max_diff, diff);
            }
        }
    }

    std::cout << "Maximum difference between original and reconstructed control points: "
              << std::scientific << std::setprecision(6) << max_diff << std::endl;
    assert(max_diff<std::numeric_limits<double>::epsilon());
    return 0;
}