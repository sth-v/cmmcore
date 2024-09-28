//
// Created by Andrew Astakhov on 28.09.24.
//
#include <iostream>
#include "cmmcore/nurbs_utils.h"
using namespace cmmcore;
int main() {
    int degree = 3;
    std::vector<double> knotvector = {0, 0, 0, 0, 1, 1, 1, 1};
    std::vector<std::vector<double>> ctrlpts = {
        {0.0, 0.0},
        {1.0, 2.0},
        {2.0, 0.0},
        {3.0, 2.0}
    };
    double u = 0.5;
    int num = 1;
    int s = find_multiplicity(u, knotvector);
    int span = find_span(ctrlpts.size() - 1, degree, u, knotvector, false);

    auto new_ctrlpts = knot_insertion(degree, knotvector, ctrlpts, u, num, s, span);

    // new_ctrlpts now contains the updated control points after knot insertion
    printf("[");
    for (auto& cpt : new_ctrlpts) {
        printf("[%f,%f],", cpt[0], cpt[1]);
    }
    printf("]");
    return 0;
}