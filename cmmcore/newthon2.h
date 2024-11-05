//
// Created by Andrew Astakhov on 04.11.24.
//

#ifndef NEWTHON2_H
#define NEWTHON2_H
#include <array>
#include <functional>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <iostream>

namespace cmmcore{


// Template parameter for dimension
template <size_t N>
using Vector = std::array<double, N>;

template <size_t N>
using matrix = std::array<std::array<double, N>, N>;

// Helper function to compute dot product
template <size_t N>
constexpr inline double dotProduct(const Vector<N>& a, const Vector<N>& b) {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

// Compute the gradient of f at point using central differences
template <size_t N>
constexpr Vector<N> computeGradient(const std::function<double(const Vector<N>&)>& f,
                          const Vector<N>& point, double h = 1e-7) {
    Vector<N> grad{};
    Vector<N> point_forward = point;
    Vector<N> point_backward = point;

    for (size_t i = 0; i < N; ++i) {
        double original = point[i];

        point_forward[i] = original + h;
        point_backward[i] = original - h;

        double f_forward = f(point_forward);
        double f_backward = f(point_backward);

        // Check for Inf or NaN values
        if (!std::isfinite(f_forward) || !std::isfinite(f_backward)) {
            std::cerr << "Warning: Non-finite function value encountered in gradient computation." << std::endl;
            std::cerr <<"( ";
            for (auto& v : point)
            {
                std::cerr << v<<", ";
            }
            std::cerr <<") f(point) = "<< f(point) <<std::endl;
            // Set gradient component to zero or a large finite value
            grad[i] = 0.0;
        } else {
            grad[i] = (f_forward - f_backward) / (2 * h);
        }

        // Restore original value
        point_forward[i] = original;
        point_backward[i] = original;
    }
    return grad;
}

// Compute the Hessian matrix of f at point using central differences
template <size_t N>
constexpr matrix<N> computeHessian(const std::function<double(const Vector<N>&)>& f,
                         const Vector<N>& point, double h = 1e-5) {
    matrix<N> H{};
    Vector<N> point_perturbed = point;
    double f0 = f(point);

    // Check for non-finite f0
    if (!std::isfinite(f0)) {
        std::cerr << "Warning: Non-finite function value at the base point in Hessian computation. "<< std::endl;
        std::cerr <<"( ";
        for (auto& v : point)
        {
            std::cerr << v<<", ";
        }
        std::cerr <<") f(point) = "<< f(point) <<std::endl;
        // Return identity matrix scaled by a large value
        for (size_t i = 0; i < N; ++i) {
            H[i][i] = 1e6;
        }
        return H;
    }

    // Precompute function evaluations to reuse
    std::array<double, N> f_forward{}, f_backward{};
    for (size_t i = 0; i < N; ++i) {
        point_perturbed[i] = point[i] + h;
        f_forward[i] = f(point_perturbed);

        point_perturbed[i] = point[i] - h;
        f_backward[i] = f(point_perturbed);

        // Check for non-finite values
        if (!std::isfinite(f_forward[i]) || !std::isfinite(f_backward[i])) {
            std::cerr << "Warning: Non-finite function value encountered in Hessian computation." << std::endl;
            std::cerr <<"( ";
            for (auto& v : point)
            {
                std::cerr << v<<", ";
            }
            std::cerr <<") f(point) = "<< f(point) <<std::endl;
            f_forward[i] = f_backward[i] = f0;
        }

        point_perturbed[i] = point[i]; // Restore
    }

    for (size_t i = 0; i < N; ++i) {
        // Diagonal elements
        double f_ii = f_forward[i];
        double f__ii = f_backward[i];
        H[i][i] = (f_ii - 2 * f0 + f__ii) / (h * h);

        for (size_t j = i + 1; j < N; ++j) {
            // Off-diagonal elements
            double orig_i = point[i];
            double orig_j = point[j];

            // f(x_i + h, x_j + h)
            point_perturbed[i] = orig_i + h;
            point_perturbed[j] = orig_j + h;
            double f_pp = f(point_perturbed);

            // f(x_i + h, x_j - h)
            point_perturbed[j] = orig_j - h;
            double f_pm = f(point_perturbed);

            // f(x_i - h, x_j + h)
            point_perturbed[i] = orig_i - h;
            point_perturbed[j] = orig_j + h;
            double f_mp = f(point_perturbed);

            // f(x_i - h, x_j - h)
            point_perturbed[j] = orig_j - h;
            double f_mm = f(point_perturbed);

            H[i][j] = H[j][i] = (f_pp - f_pm - f_mp + f_mm) / (4 * h * h);

            // Restore original values
            point_perturbed[i] = orig_i;
            point_perturbed[j] = orig_j;
        }
    }
    return H;
}

// Cholesky decomposition for solving H * x = b
// Returns true if successful, false if matrix is not positive definite
template <size_t N>
constexpr bool choleskyDecomposition(const matrix<N>& H, matrix<N>& L) {
    L = {}; // Initialize L to zero

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double sum = H[i][j];

            for (size_t k = 0; k < j; ++k)
                sum -= L[i][k] * L[j][k];

            if (i == j) {
                if (sum <= 0.0)
                    return false; // Not positive definite
                L[i][j] = std::sqrt(sum);
            } else {
                L[i][j] = sum / L[j][j];
            }
        }
    }
    return true;
}

// Solve L * y = b
template <size_t N>
constexpr void forwardSubstitution(const matrix<N>& L, const Vector<N>& b, Vector<N>& y) {
    for (size_t i = 0; i < N; ++i) {
        double sum = b[i];
        for (size_t k = 0; k < i; ++k)
            sum -= L[i][k] * y[k];
        y[i] = sum / L[i][i];
    }
}

// Solve L^T * x = y
template <size_t N>
constexpr void backwardSubstitution(const matrix<N>& L, const Vector<N>& y, Vector<N>& x) {
    for (int i = N - 1; i >= 0; --i) {
        double sum = y[i];
        for (size_t k = i + 1; k < N; ++k)
            sum -= L[k][i] * x[k];
        x[i] = sum / L[i][i];
    }
}

// Solve H * x = b using Cholesky decomposition
template <size_t N>
constexpr bool solveLinearSystem(matrix<N> H, const Vector<N>& b, Vector<N>& x) {
    matrix<N> L;
    if (!choleskyDecomposition(H, L)) {
        return false; // Not positive definite
    }

    Vector<N> y;
    forwardSubstitution(L, b, y);
    backwardSubstitution(L, y, x);

    return true;
}

// Backtracking line search with Armijo condition
template <size_t N>
constexpr double lineSearch(const std::function<double(const Vector<N>&)>& f,
                  const Vector<N>& point, const Vector<N>& direction,
                  const Vector<N>& grad, double alpha_init = 1.0,
                  double rho = 0.5, double c = 1e-4) {
    double alpha = alpha_init;
    double f0 = f(point);
    double grad_dot_dir = dotProduct(grad, direction);
    Vector<N> new_point;

    while (true) {
        for (size_t i = 0; i < N; ++i)
            new_point[i] = point[i] + alpha * direction[i];
        double f_new = f(new_point);
        if (f_new <= f0 + c * alpha * grad_dot_dir)
            break;
        alpha *= rho;
    }
    return alpha;
}

// Newton's method implementation to find the root (minimum) of a function
template <size_t N>
constexpr Vector<N> newtonsMethod(const std::function<double(const Vector<N>&)>& f,
                        Vector<N> point,
                        double tol = 1e-5, int maxIter = 100) {
    Vector<N> grad{};
    matrix<N> H{};
    Vector<N> step{};

    for (int iter = 0; iter < maxIter; ++iter) {
        grad = computeGradient<N>(f, point);

        // Check for non-finite gradient
        if (std::any_of(grad.begin(), grad.end(), [](double v){ return !std::isfinite(v); })) {
            std::cerr << "Error: Non-finite gradient encountered at iteration " << iter << "." << std::endl;
            std::cerr <<"( ";
            for (auto& v : point)
            {
                std::cerr << v<<", ";
            }
            std::cerr <<") f(point) = "<< f(point) <<std::endl;
            return point;
        }

        H = computeHessian<N>(f, point);

        // Check for non-finite Hessian
        bool hessian_finite = true;
        for (const auto& row : H) {
            for (const auto& val : row) {
                if (!std::isfinite(val)) {
                    hessian_finite = false;
                    break;
                }
            }
            if (!hessian_finite)
                break;
        }

        if (!hessian_finite) {
            std::cerr << "Error: Non-finite Hessian encountered at iteration " << iter << "." << std::endl;
            return point;
        }

        // Modify Hessian to be positive definite if necessary
        double mu = 1e-6;
        bool success = false;
        while (mu <= 1e6) {
            matrix<N> H_mod = H;
            for (size_t i = 0; i < N; ++i) {
                H_mod[i][i] += mu;
            }
            success = solveLinearSystem<N>(H_mod, grad, step);
            if (success)
                break;
            mu *= 10;
        }

        if (!success) {
            std::cerr << "Error: Unable to solve the linear system at iteration " << iter << "." << std::endl;
            return point;
        }

        // Step direction is negative of the solution
        for (size_t i = 0; i < N; ++i)
            step[i] = -step[i];

        // Line search to find acceptable step size
        double alpha = lineSearch<N>(f, point, step, grad);

        if (alpha == 0.0) {
            std::cerr << "Warning: Alpha is zero at iteration " << iter << ", stopping optimization." << std::endl;
            return point;
        }

        // Update the point: point = point + alpha * step
        bool converged = true;
        for (size_t i = 0; i < N; ++i) {
            double delta = alpha * step[i];
            point[i] += delta;
            if (std::abs(delta) >= tol * std::max(1.0, std::abs(point[i]))) {
                converged = false;
            }
        }

        if (converged) {
            // std::cout << "Converged in " << iter + 1 << " iterations." << std::endl;
            return point;
        }
    }

    // std::cerr << "Iteration limit reached without convergence." << std::endl;
    return point;
}


}
#endif //NEWTHON2_H
