//
// Created by Andrew Astakhov on 08.10.24.
//



#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "cmmcore/monomial.h"
#include <cmmcore/nurbs.h>
using namespace cmmcore;

int main() {
    // Create a simple 3x3 Bézier surface patch
     NURBSSurface surf = NURBSSurface({{{0.62332799, 0.21476556, 0.75989061, 1.        },
       {0.13268676, 0.85092707, 0.19604599, 1.        },
       {0.89074325, 0.64847999, 0.07050099, 1.        },
       {0.91934447, 0.41607874, 0.52330899, 1.        }},
      {{0.50844946, 0.88889217, 0.78646248, 1.        },
       {0.25824996, 0.33773866, 0.46611303, 1.        },
       {0.7197276 , 0.5685781 , 0.05564343, 1.        },
       {0.13142747, 0.81370473, 0.33629716, 1.        }},
      {{0.75692357, 0.64838686, 0.14546997, 1.        },
       {0.35443106, 0.60798288, 0.11322816, 1.        },
       {0.17406329, 0.12349951, 0.91478453, 1.        },
       {0.06059289, 0.53419607, 0.10507294, 1.        }},
      {{0.05756793, 0.59932455, 0.02830079, 1.        },
       {0.17660584, 0.8646054 , 0.62055444, 1.        },
       {0.42505789, 0.73934683, 0.41718578, 1.        },
       {0.33858193, 0.57782595, 0.53923527, 1.        }}}, {3,3});




    auto mono= Monomial2D(surf);

    std::cout << "Monomial coefficients:" << std::endl;
    printf("[" );
    for (auto& mono_coeff : mono.coefficients) {

        printf((format_vec3vec(mono_coeff)+",").c_str());
    }
    printf("],\n" );

    NURBSSurface surf2;
    mono.to_bezier(surf2);

    // Convert monomial back to Bézier form
    Monomial2D monoNormal;
    mono.computeNormal(monoNormal);
    std::cout << "Monomial normal coefficients:" << std::endl;
    printf("[" );
    for (auto& mono_coeff : monoNormal.coefficients) {

        printf((format_vec3vec(mono_coeff)+",").c_str());
    }
    printf("],\n" );

    // Calculate and print the maximum difference between original and reconstructed control points
    double max_diff = 0.0;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {

                double diff = std::abs((surf._control_points[i][j].to_vec3() -surf._control_points[i][j].to_vec3()).length());
                max_diff = std::max(max_diff, diff);

        }
    }

    std::cout << "Maximum difference between original and reconstructed control points: "
              << std::scientific << std::setprecision(6) << max_diff << std::endl;
    assert(max_diff<std::numeric_limits<double>::epsilon());
    return 0;
}