//
// Created by Andrew Astakhov on 26.10.24.
//

#ifndef CMMCORE_NEWTHON
#define CMMCORE_NEWTHON
#include <numeric>
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <limits>

namespace cmmcore {
    // Utility function to compute the gradient of a scalar function
    inline std::vector<double> computeGradient(const std::function<double(const std::vector<double>&)>& f,
                                               const std::vector<double>& point, double h = 1e-5) {
        size_t n = point.size();
        std::vector<double> grad(n);
        for (size_t i = 0; i < n; ++i) {
            std::vector<double> pointForward = point;
            std::vector<double> pointBackward = point;
            pointForward[i] += h;
            pointBackward[i] -= h;
            grad[i] = (f(pointForward) - f(pointBackward)) / (2 * h);
        }
        return grad;
    }
    // Utility function to compute the Hessian matrix of a scalar function
    inline std::vector<std::vector<double>> computeHessian(const std::function<double(const std::vector<double>&)>& f,
                                                           const std::vector<double>& point, double h = 1e-5) {
        size_t n = point.size();
        std::vector<std::vector<double>> H(n, std::vector<double>(n, 0.0));
        double fp = f(point);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                if (i == j) {
                    std::vector<double> forward = point;
                    std::vector<double> backward = point;
                    forward[i] += h;
                    backward[i] -= h;
                    H[i][i] = (f(forward) - 2 * fp + f(backward)) / (h * h);
                } else {
                    std::vector<double> forward_i_j = point;
                    std::vector<double> forward_i_backward_j = point;
                    std::vector<double> backward_i_forward_j = point;
                    std::vector<double> backward_i_j = point;
                    forward_i_j[i] += h;
                    forward_i_j[j] += h;
                    forward_i_backward_j[i] += h;
                    forward_i_backward_j[j] -= h;
                    backward_i_forward_j[i] -= h;
                    backward_i_forward_j[j] += h;
                    backward_i_j[i] -= h;
                    backward_i_j[j] -= h;
                    H[i][j] = H[j][i] = (f(forward_i_j) - f(forward_i_backward_j) - f(backward_i_forward_j) + f(backward_i_j)) / (4 * h * h);
                }
            }
        }
        return H;
    }
    // Newton's method implementation to find the root of a function
    inline int newtonsMethod(const std::function<double(const std::vector<double>&)>& f,
                             std::vector<double> point,
                             double tol = 1e-5, int maxIter = 100) {
        const size_t n = point.size();
        for (int iter = 0; iter < maxIter; ++iter) {
            std::vector<double> grad = computeGradient(f, point);
            std::vector<std::vector<double>> H = computeHessian(f, point);
            double det = H[0][0] * H[1][1] - H[0][1] * H[1][0];
            if (std::abs(det) < std::numeric_limits<double>::epsilon()) {
                std::cerr << "Warning: Hessian is singular at iteration " << iter << "." << std::endl;
                return 1; // Return point as the best guess
            }
            std::vector<std::vector<double>> H_inv(2, std::vector<double>(2));
            H_inv[0][0] = H[1][1] / det;
            H_inv[0][1] = -H[0][1] / det;
            H_inv[1][0] = -H[1][0] / det;
            H_inv[1][1] = H[0][0] / det;
            std::vector<double> step(n, 0.0);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    step[i] += H_inv[i][j] * grad[j];
                }
            }
            for (size_t i = 0; i < n; ++i) {
                point[i] -= step[i];
            }
            double norm = std::sqrt(std::inner_product(step.begin(), step.end(), step.begin(), 0.0));
            if (norm < tol) {
                std::cout << "Converged in " << iter << " iterations." << std::endl;
                return 0;
            }
        }
        std::cerr << "Iteration limit reached without convergence." << std::endl;
        return 1;
    }
}

#endif //CMMCORE_NEWTHON
