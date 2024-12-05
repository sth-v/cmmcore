//
// Ported _minimize_scalar_bounded function from scipy/optimize/_optimize.py
//

#ifndef CMMCORE_MINIMIZE_SCALAR_H
#define CMMCORE_MINIMIZE_SCALAR_H
#include <cmath>
#include <limits>
#include <string>
#include <functional>
#include <ostream>
namespace cmmcore{


template <typename T>
struct OptimizeResult {
    T fun;      // Value of objective function at minimum
    T x;        // Location of minimum
    int status;      // Status of optimization (0=success)
    bool success;    // Whether optimization succeeded
    std::string message;  // Description of status
    int nfev;       // Number of function evaluations
    int nit;        // Number of iterations

    friend std::ostream& operator<<(std::ostream& os, const OptimizeResult& result) {
        os << "message: " << result.message << "\n"
           << "success: " << (result.success ? "True" : "False") << "\n"
           << "status: " << result.status << "\n"
           << "fun: " << result.fun << "\n"
           << "x: " << result.x << "\n"
           << "nit: " << result.nit << "\n"
           << "nfev: " << result.nfev;
        return os;
    }
};
template <typename T>
bool is_finite_scalar(T x) {
    return std::isfinite<T>(x);
}
template <typename T>
void minimize_scalar_bounded(
    const std::function<T(T)>& func,
    std::pair<T, T> bounds,
    OptimizeResult<T>& result,
    T xatol = 1e-5,
    int maxiter = 500,
    int disp = 0
) {
    int maxfun = maxiter;
    T x1 = bounds.first;
    T x2 = bounds.second;

    // Validate bounds
    if (!is_finite_scalar(x1) || !is_finite_scalar(x2)) {
        throw std::invalid_argument("Optimization bounds must be finite scalars.");
    }

    if (x1 > x2) {
        throw std::invalid_argument("The lower bound exceeds the upper bound.");
    }

    int flag = 0;
    const T sqrt_eps = std::sqrt(2.2e-16);
    const T golden_mean = 0.5 * (3.0 - std::sqrt(5.0));

    T a = x1;
    T b = x2;
    T fulc = a + golden_mean * (b - a);
    T nfc = fulc, xf = fulc;
    T rat = 0.0, e = 0.0;
    T x = xf;
    T fx = func(x);
    int num = 1;
    T fu = std::numeric_limits<T>::infinity();

    T ffulc = fx;
    T fnfc = fx;
    T xm = 0.5 * (a + b);
    T tol1 = sqrt_eps * std::abs(xf) + xatol / 3.0;
    T tol2 = 2.0 * tol1;

    if (disp > 2) {
        printf(" Func-count     x          f(x)          Procedure\n");
        printf("%5d   %12.6g %12.6g %s\n", 1, xf, fx, "initial");
    }

    // Main optimization loop
    while (std::abs(xf - xm) > (tol2 - 0.5 * (b - a))) {
        bool golden = true;

        // Check for parabolic fit
        if (std::abs(e) > tol1) {
            golden = false;
            T r = (xf - nfc) * (fx - ffulc);
            T q = (xf - fulc) * (fx - fnfc);
            T p = (xf - fulc) * q - (xf - nfc) * r;
            q = 2.0 * (q - r);

            if (q > 0.0) {
                p = -p;
            }
            q = std::abs(q);
            T r1 = e;
            e = rat;

            // Check for acceptability of parabola
            if ((std::abs(p) < std::abs(0.5 * q * r1)) &&
                (p > q * (a - xf)) &&
                (p < q * (b - xf))) {
                rat = p / q;
                x = xf + rat;

                if ((x - a) < tol2 || (b - x) < tol2) {
                    T si = std::copysign(1.0, xm - xf);
                    rat = tol1 * si;
                }
            } else {
                golden = true;
            }
        }

        // Golden section step
        if (golden) {
            if (xf >= xm) {
                e = a - xf;
            } else {
                e = b - xf;
            }
            rat = golden_mean * e;
        }

        T si = std::copysign(1.0, rat);
        x = xf + si * std::max(std::abs(rat), tol1);
        fu = func(x);
        num++;

        if (disp > 2) {
            printf("%5d   %12.6g %12.6g %s\n", num, x, fu,
                   golden ? "golden" : "parabolic");
        }

        if (fu <= fx) {
            if (x >= xf) {
                a = xf;
            } else {
                b = xf;
            }
            fulc = nfc;
            ffulc = fnfc;
            nfc = xf;
            fnfc = fx;
            xf = x;
            fx = fu;
        } else {
            if (x < xf) {
                a = x;
            } else {
                b = x;
            }
            if (fu <= fnfc || nfc == xf) {
                fulc = nfc;
                ffulc = fnfc;
                nfc = x;
                fnfc = fu;
            } else if (fu <= ffulc || fulc == xf || fulc == nfc) {
                fulc = x;
                ffulc = fu;
            }
        }

        xm = 0.5 * (a + b);
        tol1 = sqrt_eps * std::abs(xf) + xatol / 3.0;
        tol2 = 2.0 * tol1;

        if (num >= maxfun) {
            flag = 1;
            break;
        }
    }

    if (std::isnan(xf) || std::isnan(fx) || std::isnan(fu)) {
        flag = 2;
    }

    std::string message;
    switch (flag) {
        case 0:
            message = "Solution found.";
            break;
        case 1:
            message = "Maximum number of function calls reached.";
            break;
        case 2:
            message = "NaN result encountered.";
            break;
        default:
            message = "Unknown error.";
    }


    result.fun = fx;
    result.status = flag;
    result.success = (flag == 0);
    result.message = message;
    result.x = xf;
    result.nfev = num;
    result.nit = num;

    return ;
}
}
#endif //CMMCORE_MINIMIZE_SCALAR_H
