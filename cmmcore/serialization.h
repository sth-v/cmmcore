//
// Created by Andrew Astakhov on 29.09.24.
//

#ifndef SERIALIZATION_H
#define SERIALIZATION_H

#include <vector>
#include <string>
#include <fstream>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include "cmmcore/vec.h"
#include "cmmcore/nurbs.h"

namespace cmmcore {
    // Serialize an integer
    void serializeInt(std::vector<uint8_t> &buffer, int value) {
        uint8_t *data = reinterpret_cast<uint8_t *>(&value);
        buffer.insert(buffer.end(), data, data + sizeof(int));
    }

    // Deserialize an integer
    int deserializeInt(const uint8_t *&data) {
        int value;
        std::memcpy(&value, data, sizeof(int));
        data += sizeof(int);
        return value;
    }

    // Serialize a double
    void serializeDouble(std::vector<uint8_t> &buffer, double value) {
        uint8_t *data = reinterpret_cast<uint8_t *>(&value);
        buffer.insert(buffer.end(), data, data + sizeof(double));
    }

    // Deserialize a double
    double deserializeDouble(const uint8_t *&data) {
        double value;
        std::memcpy(&value, data, sizeof(double));
        data += sizeof(double);
        return value;
    }

    // Serialize a bool
    void serializeBool(std::vector<uint8_t> &buffer, bool value) {
        buffer.push_back(static_cast<uint8_t>(value));
    }

    // Deserialize a bool
    bool deserializeBool(const uint8_t *&data) {
        bool value = static_cast<bool>(*data);
        data += sizeof(uint8_t);
        return value;
    }

    // Serialize a vector of doubles
    void serializeVectorDouble(std::vector<uint8_t> &buffer, const std::vector<double> &vec) {
        serializeInt(buffer, vec.size());
        for (double d: vec) {
            serializeDouble(buffer, d);
        }
    }

    // Deserialize a vector of doubles
    std::vector<double> deserializeVectorDouble(const uint8_t *&data) {
        int size = deserializeInt(data);
        std::vector<double> vec(size);
        for (int i = 0; i < size; ++i) {
            vec[i] = deserializeDouble(data);
        }
        return vec;
    }

    // Serialize a vec4
    void serializeVec4(std::vector<uint8_t> &buffer, const vec4 &v) {
        for (int i = 0; i < 4; ++i) {
            serializeDouble(buffer, v[i]);
        }
    }

    // Deserialize a vec4
    vec4 deserializeVec4(const uint8_t *&data) {
        vec4 v;
        for (int i = 0; i < 4; ++i) {
            v[i] = deserializeDouble(data);
        }
        return v;
    }

    // Serialize a vector of vec4
    void serializeVectorVec4(std::vector<uint8_t> &buffer, const std::vector<vec4> &vec) {
        serializeInt(buffer, vec.size());
        for (const vec4 &v: vec) {
            serializeVec4(buffer, v);
        }
    }

    // Deserialize a vector of vec4
    std::vector<vec4> deserializeVectorVec4(const uint8_t *&data) {
        int size = deserializeInt(data);
        std::vector<vec4> vec(size);
        for (int i = 0; i < size; ++i) {
            vec[i] = deserializeVec4(data);
        }
        return vec;
    }

    // Serialize a NURBSCurve object
    std::vector<uint8_t> serializeNURBSCurve(NURBSCurve &curve) {
        std::vector<uint8_t> buffer;
        serializeInt(buffer, curve.get_degree());
        serializeBool(buffer, curve.is_periodic());
        serializeVectorVec4(buffer, curve.get_control_points());
        serializeVectorDouble(buffer, curve.get_knots());
        return buffer;
    }

    // Deserialize a NURBSCurve object
    NURBSCurve deserializeNURBSCurve(const uint8_t *&data) {
        int degree = deserializeInt(data);
        bool periodic = deserializeBool(data);
        std::vector<vec4> control_points = deserializeVectorVec4(data);
        std::vector<double> knots = deserializeVectorDouble(data);
        return NURBSCurve(control_points, degree, knots, periodic);
    }

    // Serialize a vector of NURBSCurve objects
    std::vector<uint8_t> serializeNURBSCurves(std::vector<NURBSCurve> &curves) {
        std::vector<uint8_t> buffer;
        serializeInt(buffer, curves.size());
        for (NURBSCurve &curve: curves) {
            std::vector<uint8_t> curve_data = serializeNURBSCurve(curve);
            buffer.insert(buffer.end(), curve_data.begin(), curve_data.end());
        }
        return buffer;
    }

    // Deserialize a vector of NURBSCurve objects
    std::vector<NURBSCurve> deserializeNURBSCurves(const std::vector<uint8_t> &data) {
        const uint8_t *data_ptr = data.data();
        int size = deserializeInt(data_ptr);
        std::vector<NURBSCurve> curves(size);
        for (int i = 0; i < size; ++i) {
            curves[i] = deserializeNURBSCurve(data_ptr);
        }
        return curves;
    }

    // Write a vector of NURBSCurve objects to a file
    bool writeNURBSCurvesToFile(const std::string &filename, std::vector<NURBSCurve> &curves) {
        std::vector<uint8_t> data = serializeNURBSCurves(curves);
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            return false;
        }
        file.write(reinterpret_cast<const char *>(data.data()), data.size());
        return file.good();
    }

    // Read a vector of NURBSCurve objects from a file
    inline std::vector<NURBSCurve> readNURBSCurvesFromFile(const std::string &filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Could not open file for reading");
        }
        std::vector<uint8_t> data((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        return deserializeNURBSCurves(data);
    }
}

#endif //SERIALIZATION_H
