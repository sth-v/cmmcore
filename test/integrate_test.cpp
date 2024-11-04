//
// Created by Andrew Astakhov on 04.11.24.
//
#include "cmmcore/integrate.h"

#include <cassert>
#include <iostream>
#include <cmmcore/nurbs.h>


using namespace cmmcore;

int main()
{
    double result = 0.0;
    double error = 0.0;
    // Example function to integrate
    auto func = [](double t) -> double { return std::sin(t); };
    // Wrap the function in a std::function
    constexpr double tol=1e-8;

    integrate(func, 0.0, M_PI, result, error, tol);

    std::cout << "Integral = " << result << std::endl;
    std::cout << "Estimated error = " << error << std::endl;
    assert (std::abs(result - 2.0) < tol);



    double result2 = find_t1_newton(0.0, 2.,func, 2.5);
    std::cout << "Integral from t0 to t1 = " << result2 << std::endl;

    NURBSCurve curve={{{82.590441760000004, 93.433423149999996, 1.3664691899999999, 1.}, {74.593468410000000, 107.28458731000001, 11.677832530000000, 1.}, {67.794196490000004, 107.28458731000001, 11.677832530000000, 1.}, {66.055171930000000, 100.79445928000000, 18.428662989999999, 1.}, {79.915180860000007, 94.044880919999997, 8.4899551199999994, 1.}, {80.435717909999994, 90.849708669999998, 7.3044441400000002, 1.}, {77.964498780000000, 89.422949639999999, 8.9282437800000007, 1.}, {71.405593510000003, 93.209735370000004, 18.662592450000002, 1.}, {64.604393130000005, 97.136410229999996, 11.677832530000000, 1.}, {60.529016300000002, 90.077650500000004, 7.3565421100000004, 1.}, {61.852706509999997, 85.137571379999997, 11.677832530000000, 1.}, {65.517514250000005, 86.119553659999994, 11.677832530000000, 1.}, {81.279475120000001, 88.787741519999997, 16.404133649999999, 1.}, {82.966008479999999, 86.544302020000003, 11.677832530000000, 1.}, {74.432499269999994, 87.204450989999998, 11.677832530000000, 1.}, {72.507462320000002, 80.052285970000000, 11.677832530000000, 1.}, {71.159364179999997, 76.195238169999996, 8.4119728499999997, 1.}, {84.598869510000000, 74.761491860000007, 11.677832530000000, 1.}}, 3};
    auto timer=Timer();
    timer.start();
    double total_length=curve.length();
    timer.stop();
    timer.print("curve length at: ");
    std::cout << "Curve length = " << total_length << std::endl;
    timer.start();
    double param=curve.lengthAt(total_length/3, 1e-3);
    timer.stop();
    timer.print("curve lengthAt at: ");
    std::cout << "Param at:"<< (total_length/3) <<" = "  <<param<<std::endl;
    std::cout << "Param min/max: "<< (curve._interval[0]) <<" = "  <<curve._interval[1]<<std::endl;


}
