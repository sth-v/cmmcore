from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from cmmcore.vec cimport vec4,vec3
from cmmcore.nurbs cimport NURBSCurve
cdef extern from "serialization.h"  namespace "cmmcore" nogil:

    # Function declarations
    void serializeInt(vector[uint8_t]& buffer, int value)    
    int deserializeInt(const uint8_t*& data)

    void serializeDouble(vector[uint8_t]& buffer, double value)
    double deserializeDouble(const uint8_t*& data)

    void serializeBool(vector[uint8_t]& buffer, bool value)
    bool deserializeBool(const uint8_t*& data)

    void serializeVectorDouble(vector[uint8_t]& buffer, const vector[double]& vec)
    vector[double] deserializeVectorDouble(const uint8_t*& data)

    void serializeVec4(vector[uint8_t]& buffer, const vec4& v)
    vec4 deserializeVec4(const uint8_t*& data)

    void serializeVectorVec4(vector[uint8_t]& buffer, const vector[vec4]& vec)
    vector[vec4] deserializeVectorVec4(const uint8_t*& data)

    vector[uint8_t] serializeNURBSCurve(NURBSCurve& curve)
    NURBSCurve deserializeNURBSCurve(const uint8_t*& data)

    vector[uint8_t] serializeNURBSCurves(vector[NURBSCurve]& curves)
    vector[NURBSCurve] deserializeNURBSCurves(const vector[uint8_t]& data)

    bool writeNURBSCurvesToFile(const string& filename, vector[NURBSCurve]& curves)
    vector[NURBSCurve] readNURBSCurvesFromFile(const string& filename)