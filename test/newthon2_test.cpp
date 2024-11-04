//
// Created by Andrew Astakhov on 04.11.24.
//
// Example usage
#include <iostream>

#define CMMCORE_DEBUG
#include <cassert>
#include <cmmcore/nurbs.h>
#include <cmmcore/vec.h>
#include "cmmcore/newthon2.h"
using namespace cmmcore;

int main() {
    constexpr size_t N = 2;

    // Define a sample function: Rosenbrock function
    auto rosenbrock = [](const Vector<N>& x) -> double {
        constexpr double a = 1.0;
        constexpr double b = 100.0;
        double term1 = a - x[0];
        double term2 = x[1] - x[0] * x[0];
        return term1 * term1 + b * term2 * term2;
    };

    // Initial guess
    Vector<N> initial_point = { -1.2, 1.0 };

    // Perform Newton's method
    Vector<N> minimum = newtonsMethod<N>(rosenbrock, initial_point);

    // Output the result
    std::cout << "Minimum found at: ";
    for (const auto& val : minimum) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    std::vector<std::vector<cmmcore::vec4>> pts1={{{-25.0, -25.0, -10.0, 1.0}, {-25.0, -15.0, -5.0, 1.0}, {-25.0, -5.0, 0.0, 1.0}, {-25.0, 5.0, 0.0, 1.0}, {-25.0, 15.0, -5.0, 1.0}, {-25.0, 25.0, -10.0, 1.0}}, {{-15.0, -25.0, -8.0, 1.0}, {-15.0, -15.0, -4.0, 1.0}, {-15.0, -5.0, -4.0, 1.0}, {-15.0, 5.0, -4.0, 1.0}, {-15.0, 15.0, -4.0, 1.0}, {-15.0, 25.0, -8.0, 1.0}}, {{-5.0, -25.0, -5.0, 1.0}, {-5.0, -15.0, -3.0, 1.0}, {-5.0, -5.0, -8.0, 1.0}, {-5.0, 5.0, -8.0, 1.0}, {-5.0, 15.0, -3.0, 1.0}, {-5.0, 25.0, -5.0, 1.0}}, {{5.0, -25.0, -3.0, 1.0}, {5.0, -15.0, -2.0, 1.0}, {5.0, -5.0, -8.0, 1.0}, {5.0, 5.0, -8.0, 1.0}, {5.0, 15.0, -2.0, 1.0}, {5.0, 25.0, -3.0, 1.0}}, {{15.0, -25.0, -8.0, 1.0}, {15.0, -15.0, -4.0, 1.0}, {15.0, -5.0, -4.0, 1.0}, {15.0, 5.0, -4.0, 1.0}, {15.0, 15.0, -4.0, 1.0}, {15.0, 25.0, -8.0, 1.0}}, {{25.0, -25.0, -10.0, 1.0}, {25.0, -15.0, -5.0, 1.0}, {25.0, -5.0, 2.0, 1.0}, {25.0, 5.0, 2.0, 1.0}, {25.0, 15.0, -5.0, 1.0}, {25.0, 25.0, -10.0, 1.0}}};
    std::array<int,2> deg1={3,3};
    cmmcore::NURBSSurface ns1(pts1,deg1);
    auto pt=cmmcore::vec3(0.,1.,5);
    auto function = [&pt,&ns1](const Vector<N>& v) {
        cmmcore::vec3 vv;
        ns1.evaluate(v[0],v[1],vv);

        return (pt-vv).sqLength();
    };
    auto timer=Timer();
    Vector<N> initialPoint = {1.0, 2.0};
    timer.start();
    Vector<N> result = cmmcore::newtonsMethod<N>(function, initialPoint);
    timer.stop();
    timer.print("newton2 at: ");

    std::cout << "Root found at: (" << result[0] << ", " << result[1] << ")" << std::endl;
    assert(result[0]-1.51302<1e-5&& result[1]-1.79546<1e-5);
    return 0;
}