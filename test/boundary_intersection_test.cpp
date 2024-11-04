//
// Created by Andrew Astakhov on 04.11.24.
//
#include "cmmcore/boundary_intersection.h"

using namespace cmmcore;

#define create_test_surface1() NURBSSurface({                                           \
{{0., 0., 0.}, {0., 1., 0.}, {0., 2., 0.}}, \
{{1., 0., 0.}, {1., 1., 1.}, {1., 2., 0.}}, \
{{2., 0., 0.}, {2., 1., 0.}, {2., 2., 0.}}  \
}, {2, 2} )                                      \

#define create_test_surface2() NURBSSurface({                                           \
{{1., -1., -1.}, {1., 0., 1.}, {1., 1., -1.}},\
{{1., -1., 0.}, {1., 0., 2.}, {1., 1., 0.}},\
{{1., -1., 1.}, {1., 0., 1.}, {1., 1., 1.}}\
}, {2, 2} )                                      \

void test_extract_surface_boundaries(){
    auto surface = create_test_surface1();
    std::array<NURBSCurve,4> boundaries{};
    extract_surface_boundaries(surface,boundaries);
    vec3 res1,res2;
    boundaries[0].evaluate(0.,res1);
    boundaries[2].evaluate(0,res2);
    assert ((res1-res2)<1e-6);

    boundaries[0].evaluate(1.,res1);
    boundaries[3].evaluate(0,res2);
    assert ((res1-res2)<1e-6);

    boundaries[1].evaluate(0.,res1);
    boundaries[2].evaluate(1,res2);
    assert ((res1-res2)<1e-6);

    boundaries[1].evaluate(1.,res1);
    boundaries[3].evaluate(1,res2);
    assert ((res1-res2)<1e-6);


}

int main() {
  test_extract_surface_boundaries();
  return 0;


}