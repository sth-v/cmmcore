//
// Created by Andrew Astakhov on 09.10.24.
//

#include "cmmcore/ssx.h"
using namespace cmmcore;
inline void printInter(cmmcore::Intersection& ixs) {
    printf("Intersection:\n");
    printf("[");
    for (int i=0; i<ixs.patches.size(); i++) {
        //printf("[[");
        //for (auto& pt :ixs.patches[i].s1.control_points_flat3d()) {
        //    printf("[%f,%f,%f],", pt.x,pt.y,pt.z);
        //}
        //printf("],[");
        //for (auto& pt :ixs.patches[i].s2.control_points_flat3d()) {
        //    printf("[%f,%f,%f],", pt.x,pt.y,pt.z);
        //}
        //printf("]],");
    }
    printf("]\n");
    printf("[");
    for (int i=0; i<ixs.points.size(); i++) {

        printf("[%f,%f,%f],", ixs.points[i].xyz.x, ixs.points[i].xyz.y, ixs.points[i].xyz.z);
    };
    printf("]\n");
    printf("%d\n", cmmcore::counter);
}
inline void printInter2(cmmcore::Intersection& ixs) {
    printf("Intersection Pts:\n");
    printf("[");
    for (int i=0; i<ixs.patches.size(); i++) {
        printf("[");
                   printf("[%f,%f,%f],",  ixs.patches[i].start.point.x,ixs.patches[i].start.point.y,ixs.patches[i].start.point.z);
        printf("[%f,%f,%f]",  ixs.patches[i].end.point.x,ixs.patches[i].end.point.y,ixs.patches[i].end.point.z);
        printf("],");
    }
    printf("]\n");

}
inline void printInterPoints(std::vector<cmmcore::IntersectionPoint>& ixs) {
    printf("Intersection:\n");
    printf("[");
    for (auto& pt :ixs) {
            printf("[%f,%f,%f],", pt.point.x,pt.point.y,pt.point.z);
        }
    printf("]\n");



}
#define flat_index_2d(i,j, ncols) ()

int main() {
        std::vector<std::vector<cmmcore::vec4>> pts1={{{-25.0, -25.0, -10.0, 1.0}, {-25.0, -15.0, -5.0, 1.0}, {-25.0, -5.0, 0.0, 1.0}, {-25.0, 5.0, 0.0, 1.0}, {-25.0, 15.0, -5.0, 1.0}, {-25.0, 25.0, -10.0, 1.0}}, {{-15.0, -25.0, -8.0, 1.0}, {-15.0, -15.0, -4.0, 1.0}, {-15.0, -5.0, -4.0, 1.0}, {-15.0, 5.0, -4.0, 1.0}, {-15.0, 15.0, -4.0, 1.0}, {-15.0, 25.0, -8.0, 1.0}}, {{-5.0, -25.0, -5.0, 1.0}, {-5.0, -15.0, -3.0, 1.0}, {-5.0, -5.0, -8.0, 1.0}, {-5.0, 5.0, -8.0, 1.0}, {-5.0, 15.0, -3.0, 1.0}, {-5.0, 25.0, -5.0, 1.0}}, {{5.0, -25.0, -3.0, 1.0}, {5.0, -15.0, -2.0, 1.0}, {5.0, -5.0, -8.0, 1.0}, {5.0, 5.0, -8.0, 1.0}, {5.0, 15.0, -2.0, 1.0}, {5.0, 25.0, -3.0, 1.0}}, {{15.0, -25.0, -8.0, 1.0}, {15.0, -15.0, -4.0, 1.0}, {15.0, -5.0, -4.0, 1.0}, {15.0, 5.0, -4.0, 1.0}, {15.0, 15.0, -4.0, 1.0}, {15.0, 25.0, -8.0, 1.0}}, {{25.0, -25.0, -10.0, 1.0}, {25.0, -15.0, -5.0, 1.0}, {25.0, -5.0, 2.0, 1.0}, {25.0, 5.0, 2.0, 1.0}, {25.0, 15.0, -5.0, 1.0}, {25.0, 25.0, -10.0, 1.0}}};
    std::vector<std::vector<cmmcore::vec4>>pts2={{{25.0, 14.774795467423544, 5.547618997879466, 1.0}, {25.0, 10.618169208735296, -15.1325103127356, 1.0}, {25.0, 1.8288992061686002, -13.545426491756078, 1.0}, {25.0, 9.871574766108672, 14.261864686419623, 1.0}, {25.0, -15.0, 5.0, 1.0}, {25.0, -25.0, 5.0, 1.0}}, {{15.0, 25.0, 1.8481369394623908, 1.0}, {15.0, 15.0, 5.0, 1.0}, {15.0, 5.0, -1.4589623860307768, 1.0}, {15.0, -5.0, -1.9177595746260625, 1.0}, {15.0, -15.0, -30.948650572598954, 1.0}, {15.0, -25.0, 5.0, 1.0}}, {{5.0, 25.0, 5.0, 1.0}, {5.0, 15.0, -29.589097491066767, 1.0}, {3.802890818198094, 5.0, 5.0, 1.0}, {5.0, -5.0, 5.0, 1.0}, {5.0, -15.0, 5.0, 1.0}, {5.0, -25.0, 5.0, 1.0}}, {{-5.0, 25.0, 5.0, 1.0}, {-5.0, 15.0, 5.0, 1.0}, {-5.0, 5.0, 5.0, 1.0}, {-5.0, -5.0, -27.39452352115122, 1.0}, {-5.0, -15.0, 5.0, 1.0}, {-5.0, -25.0, 5.0, 1.0}}, {{-15.0, 25.0, 5.0, 1.0}, {-15.0, 15.0, -23.968082282285287, 1.0}, {-15.0, 5.0, 5.0, 1.0}, {-15.0, -5.0, 5.0, 1.0}, {-15.0, -15.0, -18.33446589106032, 1.0}, {-15.0, -25.0, 5.0, 1.0}}, {{-25.0, 25.0, 5.0, 1.0}, {-25.0, 15.0, 14.302789083068138, 1.0}, {-25.0, 5.0, 5.0, 1.0}, {-25.0, -5.0, 5.0, 1.0}, {-25.0, -15.0, 5.0, 1.0}, {-25.0, -25.0, 5.0, 1.0}}};
    std::array<int,2> deg1={3,3};
    std::array<int,2> deg2{3,3};
    cmmcore::NURBSSurface ns1(pts1,deg1);
    cmmcore::NURBSSurface ns2(pts2,deg2);
    cmmcore::Intersection intersection;
    findStartIntersectionPoints(ns1,ns2,1e-3,intersection);
    printInter(intersection);
    printInter2(intersection);

}