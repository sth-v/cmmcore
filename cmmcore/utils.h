//
// Created by Andrew Astakhov on 28.09.24.
//

#ifndef UTILS_H
#define UTILS_H
#include <cmath>
#include <limits>
namespace cmmcore {
#define EPS_DELTA 1e-12

/**
 * @brief Calculates an epsilon value for a given number.
 *
 * @param x The input number.
 * @return The calculated epsilon value.
 */
inline double calc_epsilon(const double x) noexcept {
    const double relative_epsilon = std::numeric_limits<double>::epsilon() * std::fabs(x);
    constexpr double absolute_epsilon = std::numeric_limits<double>::denorm_min();

    if ( std::fabs(x) < EPS_DELTA) {
        return absolute_epsilon;
    }
    return relative_epsilon;
}
}
#endif //UTILS_H
