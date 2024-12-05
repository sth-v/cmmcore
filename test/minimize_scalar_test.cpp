//
// Created by Andrew Astakhov on 05.12.24.
//
// Example usage
#include <iostream>
#include <cassert>
#include "cmmcore/minimize_scalar.h"
using namespace cmmcore;
double f(double x) {
    return (x - 2) * (x - 2);  // Minimum at x=2
}
float f2(float x) {
    return (x - 2) * (x - 2);  // Minimum at x=2
}

int main() {
    OptimizeResult<double> result;
    std::cout << "using double: " << "\n";
    minimize_scalar_bounded<double>(f, {0, 4}, result, 1e-5, 500, 0);
    std::cout << result << std::endl;
    assert(result.success);
    assert(result.x==2);
    assert(result.fun==0);
    if (result.success) {
        std::cout << "Minimum found at x = " << result.x << "\n";
        std::cout << "Minimum value = " << result.fun << "\n";
    }
    std::cout << "using float: " << "\n";
    OptimizeResult<float> result2;
    minimize_scalar_bounded<float>(f2, {0, 4}, result2, 1e-5, 500, 0);
    std::cout << result2 << std::endl;
    assert(result2.success);
    assert(result2.x==2);
    assert(result2.fun==0);
    if (result2.success) {
        std::cout << "Minimum found at x = " << result2.x << "\n";
        std::cout << "Minimum value = " << result2.fun << "\n";
    }
    return 0;
}