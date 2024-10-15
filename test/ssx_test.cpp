//
// Created by Andrew Astakhov on 09.10.24.
//
#include "cmmcore/ssx.h"
#include <unordered_set>
void printInter(cmmcore::Intersection& ixs) {
    printf("Intersection:\n");
    printf("[");
    for (int i=0; i<ixs.patches.size(); i++) {
        printf("[[");
        for (auto& pt :ixs.patches[i].s1.control_points_flat3d()) {
            printf("[%f,%f,%f],", pt.x,pt.y,pt.z);
        }
        printf("],[");
        for (auto& pt :ixs.patches[i].s2.control_points_flat3d()) {
            printf("[%f,%f,%f],", pt.x,pt.y,pt.z);
        }
        printf("]],");
    }
    printf("]\n");
    printf("[");
    for (int i=0; i<ixs.points.size(); i++) {

        printf("[%f,%f,%f],", ixs.points[i].xyz.x, ixs.points[i].xyz.y, ixs.points[i].xyz.z);
    };
    printf("]\n");
    printf("%d\n", cmmcore::counter);
}
#define flat_index_2d(i,j, ncols) ()
void test_deffff() {
    std::vector<std::vector<cmmcore::vec4>> pts1={{{-25.0, -25.0, -10.0, 1.0}, {-25.0, -15.0, -5.0, 1.0}, {-25.0, -5.0, 0.0, 1.0}, {-25.0, 5.0, 0.0, 1.0}, {-25.0, 15.0, -5.0, 1.0}, {-25.0, 25.0, -10.0, 1.0}}, {{-15.0, -25.0, -8.0, 1.0}, {-15.0, -15.0, -4.0, 1.0}, {-15.0, -5.0, -4.0, 1.0}, {-15.0, 5.0, -4.0, 1.0}, {-15.0, 15.0, -4.0, 1.0}, {-15.0, 25.0, -8.0, 1.0}}, {{-5.0, -25.0, -5.0, 1.0}, {-5.0, -15.0, -3.0, 1.0}, {-5.0, -5.0, -8.0, 1.0}, {-5.0, 5.0, -8.0, 1.0}, {-5.0, 15.0, -3.0, 1.0}, {-5.0, 25.0, -5.0, 1.0}}, {{5.0, -25.0, -3.0, 1.0}, {5.0, -15.0, -2.0, 1.0}, {5.0, -5.0, -8.0, 1.0}, {5.0, 5.0, -8.0, 1.0}, {5.0, 15.0, -2.0, 1.0}, {5.0, 25.0, -3.0, 1.0}}, {{15.0, -25.0, -8.0, 1.0}, {15.0, -15.0, -4.0, 1.0}, {15.0, -5.0, -4.0, 1.0}, {15.0, 5.0, -4.0, 1.0}, {15.0, 15.0, -4.0, 1.0}, {15.0, 25.0, -8.0, 1.0}}, {{25.0, -25.0, -10.0, 1.0}, {25.0, -15.0, -5.0, 1.0}, {25.0, -5.0, 2.0, 1.0}, {25.0, 5.0, 2.0, 1.0}, {25.0, 15.0, -5.0, 1.0}, {25.0, 25.0, -10.0, 1.0}}};
    std::vector<std::vector<cmmcore::vec4>>pts2={{{25.0, 14.774795467423544, 5.547618997879466, 1.0}, {25.0, 10.618169208735296, -15.1325103127356, 1.0}, {25.0, 1.8288992061686002, -13.545426491756078, 1.0}, {25.0, 9.871574766108672, 14.261864686419623, 1.0}, {25.0, -15.0, 5.0, 1.0}, {25.0, -25.0, 5.0, 1.0}}, {{15.0, 25.0, 1.8481369394623908, 1.0}, {15.0, 15.0, 5.0, 1.0}, {15.0, 5.0, -1.4589623860307768, 1.0}, {15.0, -5.0, -1.9177595746260625, 1.0}, {15.0, -15.0, -30.948650572598954, 1.0}, {15.0, -25.0, 5.0, 1.0}}, {{5.0, 25.0, 5.0, 1.0}, {5.0, 15.0, -29.589097491066767, 1.0}, {3.802890818198094, 5.0, 5.0, 1.0}, {5.0, -5.0, 5.0, 1.0}, {5.0, -15.0, 5.0, 1.0}, {5.0, -25.0, 5.0, 1.0}}, {{-5.0, 25.0, 5.0, 1.0}, {-5.0, 15.0, 5.0, 1.0}, {-5.0, 5.0, 5.0, 1.0}, {-5.0, -5.0, -27.39452352115122, 1.0}, {-5.0, -15.0, 5.0, 1.0}, {-5.0, -25.0, 5.0, 1.0}}, {{-15.0, 25.0, 5.0, 1.0}, {-15.0, 15.0, -23.968082282285287, 1.0}, {-15.0, 5.0, 5.0, 1.0}, {-15.0, -5.0, 5.0, 1.0}, {-15.0, -15.0, -18.33446589106032, 1.0}, {-15.0, -25.0, 5.0, 1.0}}, {{-25.0, 25.0, 5.0, 1.0}, {-25.0, 15.0, 14.302789083068138, 1.0}, {-25.0, 5.0, 5.0, 1.0}, {-25.0, -5.0, 5.0, 1.0}, {-25.0, -15.0, 5.0, 1.0}, {-25.0, -25.0, 5.0, 1.0}}};
    std::array<int,2> deg1={3,3};
    std::array<int,2> deg2{3,3};
    cmmcore::NURBSSurface ns1(pts1,deg1);
    cmmcore::NURBSSurface ns2(pts2,deg2);
    std::vector<cmmcore::NURBSSurface> surfs;
    std::vector<cmmcore::NURBSSurface> surfs2;


    cmmcore::decompose(ns1,surfs);
    cmmcore::decompose(ns2,surfs2);

    std::vector<cmmcore::NURBSSurface> g1;
    std::vector<cmmcore::NURBSSurface> g2;


    cmmcore::Intersection intersection{};
    for (size_t i=0;i<surfs.size();i++)
    {
        g1.push_back(cmmcore::gaussMap(surfs[i]));
    }
    for (size_t i=0;i<surfs2.size();i++)
    {
        g2.push_back(cmmcore::gaussMap(surfs2[i]));
    }

    for (size_t i=0;i<surfs.size();i++)
    {

        for (size_t j=0;j<surfs2.size();j++)
        {

                    cmmcore::detectIntersections(surfs[i],
                        surfs2[j],
                        g1[i],
                        g2[j],
                        1e-3,
                        intersection);


        }
    }



     printInter(intersection);
}
int main() {
    test_deffff();

}