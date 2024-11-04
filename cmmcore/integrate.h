//
// Created by Andrew Astakhov on 04.11.24.
//

#ifndef CMMCORE_INTEGRATE_H
#define CMMCORE_INTEGRATE_H

#include <cmath>
#include <functional>
#include <iostream>

#include <map>
#include <vector>
namespace cmmcore{


// Constants for adaptive Simpson's method
constexpr int MAX_DEPTH = 20;       // Maximum recursion depth
constexpr double TOLERANCE = 1e-6;  // Desired tolerance for error estimation

// Function prototype for the integration procedure

// Internal recursive function implementing adaptive Simpson's method
constexpr void adaptiveSimpson(const std::function<double(double)>& fun, const double a, const double b,
                                const double fa, const double fb,const  double fm,const  double S,
                                const double tol, double& result, double& error, const size_t depth, const size_t max_depth) {
        const double m = 0.5 * (a + b);
        const double h = b - a;
        const
        double m1 = 0.5 * (a + m);
        const double m2 = 0.5 * (m + b);
        const
        double f1 = fun(m1);
        const double f2 = fun(m2);
        const
        double S_left = (h / 12.0) * (fa + 4.0 * f1 + fm);
        const double S_right = (h / 12.0) * (fm + 4.0 * f2 + fb);
        const
        double S_total = S_left + S_right;
        const double error_estimate = (S_total - S) / 15.0;

        if (depth >= max_depth || std::fabs(error_estimate) < tol) {
            result += S_total;
            error += std::fabs(error_estimate);
        } else {
            adaptiveSimpson(fun, a, m, fa, fm, f1, S_left, tol / 2.0, result, error, depth + 1,max_depth);
            adaptiveSimpson(fun, m, b, fm, fb, f2, S_right, tol / 2.0, result, error, depth + 1,max_depth);
        }
    }

constexpr void integrate(const std::function<double(double)>& fun, const double t0, const double t1,
               double& result, double& error, const double tol=1e-8,const size_t max_depth=20) {
    const double fa = fun(t0);
    const double fb = fun(t1);
    const double m = 0.5 * (t0 + t1);
    const double fm = fun(m);
    const double h = t1 - t0;
    const double S = (h / 6.0) * (fa + 4.0 * fm + fb);  // Initial Simpson's rule estimate

    result = 0.0;
    error = 0.0;

    adaptiveSimpson(fun, t0, t1, fa, fb, fm, S, tol, result, error, 0, max_depth);
}

inline double find_t1_newton(double t0, double I, const std::function<double(double)>& fun,
                             double initial_guess,const double tol=1e-8, const size_t max_iter=20) {
    auto G = [&](double t) {
        double result = 0.0;
        double error = 0.0;
        integrate(fun, t0, t, result, error, tol,max_iter);
        return result - I;
    };

    auto G_prime = [&](double t) {
        return fun(t);
    };



    double t = initial_guess;

    for (int iter = 0; iter < max_iter; ++iter) {
        double Gt = G(t);
        double Gt_prime = G_prime(t);

        if (std::fabs(Gt_prime) < 1e-10) {
            std::cerr << "Derivative too small; Newton-Raphson may fail." << std::endl;
            return NAN;
        }

        double t_new = t - Gt / Gt_prime;

        if (std::fabs(t_new - t) < tol) {
            return t_new;
        }

        t = t_new;
    }

    std::cerr << "Warning: Newton-Raphson method did not converge within maximum iterations." << std::endl;
    return t;
}


}
#endif //CMMCORE_INTEGRATE_H
