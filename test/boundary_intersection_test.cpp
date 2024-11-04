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
void test_find_boundary_intersections()
{
    auto surface1 = create_test_surface1();
    auto surface2 = create_test_surface2();
    std::vector<IntersectionPoint> intersections;
    find_boundary_intersections(surface1,surface2,intersections);
    assert(intersections.size()>0);
    vec3 res1,res2;
    surface1.evaluate(0.,0.,res1);
    surface2.evaluate(0.,0.,res2);


    for (const auto& intersection : intersections)
    {
        // Check that all points lie on both surfaces using surface1_params and surface2_params
        vec3 pt1, pt2;
        surface1.evaluate(intersection.surface1_params[0], intersection.surface1_params[1], pt1);
        surface2.evaluate(intersection.surface2_params[0], intersection.surface2_params[1], pt2);

        std::cout << "Point: [" << intersection.point[0] << ", " << intersection.point[1] << ", " << intersection.point[
            2] << "]\n";
        std::cout << "Surface1 evaluation at [" << intersection.surface1_params[0] << ", " << intersection.
            surface1_params[1] << "]: ["
            << pt1[0] << ", " << pt1[1] << ", " << pt1[2] << "]\n";
        std::cout << "Surface2 evaluation at [" << intersection.surface2_params[0] << ", " << intersection.
            surface2_params[1] << "]: ["
            << pt2[0] << ", " << pt2[1] << ", " << pt2[2] << "]\n";

        double tol = 1e-6;
        assert((pt1 - intersection.point).length() < tol && "Point should lie on first surface");
        assert((pt2 - intersection.point).length() < tol && "Point should lie on second surface");
        assert((pt1 - pt2).length() < tol && "Evaluations should match each other");

        // Verify that get_start_params still works correctly for backward compatibility
        auto start_params = intersection.get_start_params();
        vec3 start_pt1, start_pt2;
        surface1.evaluate(start_params.first[0], start_params.first[1], start_pt1);
        surface2.evaluate(start_params.second[0], start_params.second[1], start_pt2);

        assert(
            (start_pt1 - intersection.point).length() < tol &&
            "get_start_params should give correct parameters for first surface");
        assert(
            (start_pt2 - intersection.point).length() < tol &&
            "get_start_params should give correct parameters for second surface");
    }

}
int main() {
  test_extract_surface_boundaries();
  test_find_boundary_intersections();
  return 0;


}