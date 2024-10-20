#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <vector>
#include <array>
#include "cmmcore/vec.h"  // Assuming you have the definitions in a header
#include "cmmcore/nurbs.h"
using namespace emscripten;
std::vector<std::vector<cmmcore::vec4>> create2DVector(int rows, int cols, cmmcore::vec4 defaultValue) {
    return std::vector<std::vector<cmmcore::vec4>>(rows, std::vector<cmmcore::vec4>(cols, defaultValue));
};
void convertJSArrayToCPPVectorDouble(val jsArray,std::vector<double>& cppVector ) {

    unsigned length = jsArray["length"].as<unsigned>();
    cppVector.resize(length);
    for (unsigned i = 0; i < length; ++i) {
        cppVector[i]=jsArray[i].as<double>();
    }

}
void convertJSArrayToNURBSSurfaceControlPoints(val jsArray,std::vector<std::vector<cmmcore::vec4>>& controlPoints ) {

    unsigned length = jsArray["length"].as<unsigned>();
    unsigned length2 = jsArray[0]["length"].as<unsigned>();
    controlPoints.resize(length, std::vector<cmmcore::vec4>(length2, {0,0,0,1}));
    for (unsigned i = 0; i < length; ++i) {
        for (unsigned j = 0; j < length2; ++j) {
        controlPoints[i][j].set(jsArray[i][j][0].as<double>(),jsArray[i][j][1].as<double>(),jsArray[i][j][2].as<double>(),jsArray[i][j][3].as<double>());
        }
    }


}

EMSCRIPTEN_BINDINGS(cmmcore_common_module) {
    function("create2DVector", &create2DVector);
    function("convertJSArrayToNURBSSurfaceControlPoints", &convertJSArrayToNURBSSurfaceControlPoints);
    function("convertJSArrayToCPPVectorDouble", &convertJSArrayToCPPVectorDouble);
}
// Binding for cmmcore::vec3
EMSCRIPTEN_BINDINGS(vec3_module) {
    class_<cmmcore::vec3>("vec3")
        .constructor<>()
        .constructor<double>()
        .constructor<double, double, double>()
        .property("x", &cmmcore::vec3::x)
        .property("y", &cmmcore::vec3::y)
        .property("z", &cmmcore::vec3::z);
}

// Binding for cmmcore::vec4
EMSCRIPTEN_BINDINGS(vec4_module) {
    class_<cmmcore::vec4>("vec4")
        .constructor<>()
        .constructor<double>()
        .constructor<double, double, double, double>()
        .property("x", &cmmcore::vec4::x)
        .property("y", &cmmcore::vec4::y)
        .property("z", &cmmcore::vec4::z)
        .property("w", &cmmcore::vec4::w);
}

// Binding for NURBSSurface
EMSCRIPTEN_BINDINGS(nurbs_surface_module) {
    // Register std::vector<double>
    register_vector<double>("vector_double");

    // Register std::vector<cmmcore::vec4>
    register_vector<cmmcore::vec4>("vector_vec4");

    // Register std::vector<std::vector<cmmcore::vec4>>
    register_vector<std::vector<cmmcore::vec4>>("vector_vector_vec4>");

    // Register std::array<int, 2>
    value_array<std::array<int, 2>>("array_int_2")
        .element(emscripten::index<0>())  // Access first element
        .element(emscripten::index<1>()); // Access second element


    // Bind the NURBSSurface class
    class_<cmmcore::NURBSSurface>("NURBSSurface")
        .constructor<>()
        .constructor<const std::vector<std::vector<cmmcore::vec4>>&, const std::array<int, 2>&, const std::vector<double>&, const std::vector<double>&>()
        .function("evaluate", &cmmcore::NURBSSurface::evaluate);
}
