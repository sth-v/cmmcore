//
// Created by Andrew Astakhov on 26.10.24.
//
#include "cmmcore/newthon2.h"
#include "cmmcore/nurbs.h"
int main() {
    std::vector<std::vector<cmmcore::vec4>> pts1={{{-25.0, -25.0, -10.0, 1.0}, {-25.0, -15.0, -5.0, 1.0}, {-25.0, -5.0, 0.0, 1.0}, {-25.0, 5.0, 0.0, 1.0}, {-25.0, 15.0, -5.0, 1.0}, {-25.0, 25.0, -10.0, 1.0}}, {{-15.0, -25.0, -8.0, 1.0}, {-15.0, -15.0, -4.0, 1.0}, {-15.0, -5.0, -4.0, 1.0}, {-15.0, 5.0, -4.0, 1.0}, {-15.0, 15.0, -4.0, 1.0}, {-15.0, 25.0, -8.0, 1.0}}, {{-5.0, -25.0, -5.0, 1.0}, {-5.0, -15.0, -3.0, 1.0}, {-5.0, -5.0, -8.0, 1.0}, {-5.0, 5.0, -8.0, 1.0}, {-5.0, 15.0, -3.0, 1.0}, {-5.0, 25.0, -5.0, 1.0}}, {{5.0, -25.0, -3.0, 1.0}, {5.0, -15.0, -2.0, 1.0}, {5.0, -5.0, -8.0, 1.0}, {5.0, 5.0, -8.0, 1.0}, {5.0, 15.0, -2.0, 1.0}, {5.0, 25.0, -3.0, 1.0}}, {{15.0, -25.0, -8.0, 1.0}, {15.0, -15.0, -4.0, 1.0}, {15.0, -5.0, -4.0, 1.0}, {15.0, 5.0, -4.0, 1.0}, {15.0, 15.0, -4.0, 1.0}, {15.0, 25.0, -8.0, 1.0}}, {{25.0, -25.0, -10.0, 1.0}, {25.0, -15.0, -5.0, 1.0}, {25.0, -5.0, 2.0, 1.0}, {25.0, 5.0, 2.0, 1.0}, {25.0, 15.0, -5.0, 1.0}, {25.0, 25.0, -10.0, 1.0}}};
    std::array<int,2> deg1={3,3};
    cmmcore::NURBSSurface ns1(pts1,deg1);
    auto pt=cmmcore::vec3(0.,1.,5);
    auto function = [&pt,&ns1](const cmmcore::Vector<2>& v)->double {
        cmmcore::vec3 vv;
        ns1.evaluate(v[0],v[1],vv);

        return (pt-vv).sqLength();
    };
    //std::vector<double> initialPoint = {1.0, 2.0};
    cmmcore::Vector<2> initialPoint = {1.0, 2.0};
    auto timer= cmmcore::Timer();
    timer.start();
    cmmcore::Vector<2> result = cmmcore::newtonsMethod<2>(function, initialPoint);
    timer.stop();
    timer.print("newthon at: ");//232000 215875
    std::cout << "Root found at: (" << result[0] << ", " << result[1] << ")" << std::endl;
    assert(result[0]-1.51302<1e-5&& result[1]-1.79546<1e-5);
    return 0;
}