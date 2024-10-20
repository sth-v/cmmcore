# cython: language_level=3
# distutils: language = c++
cimport cython
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "vec.h"  namespace "cmmcore" nogil:
    cdef cppclass vec2:
        double x, y
        vec2()
        vec2(double value)
        vec2(double v1, double v2)
        explicit vec2(double[2]& arr)
        vec2(const vec2& from, const vec2& to)
        size_t size() const
        void set(double v1)
        void set(double v1, double v2)
        void set(const vec2 v1)
        void set(const double[2]& arr)
        bool operator==(const vec2& other) const
        vec2 operator-() const
        vec2 operator-(const vec2& other) const
        vec2 operator*(const vec2& other) const
        vec2 operator*(double val) const
        void operator*=(double val)
        void operator+=(const vec2& other)
        void operator/=(double val)
        void operator-=(const vec2& other)
        double& operator[](size_t index)
        double operator[](size_t index) const
        double dot(const vec2& other) const
        double sqLength() const
        double length() const
        vec2 project(const vec2& other) const
        double projection(const vec2& other) const
        void unitize()
        vec2 unit() const
        bool collinear(const vec2& other) const

    double cross(const vec2& a, const vec2& b)
    double dot(const vec2& a, const vec2& b)

    cdef cppclass Vec2Hash:
        size_t operator()(const vec2& v) const

    cdef cppclass Vec2Equal:
        bool operator()(const vec2& a, const vec2& b) const

    cdef cppclass vec3:
        double x
        double y
        double z
        vec3()
        vec3(const double value)
        vec3(const double[ 3]& arr)
        vec3(const double x, const double y, const double z)
        size_t size() const
        void set(double x1, double y1, double z1)
        void set(double[3]& values)
        void set(const vec3& values)
        bool operator==(const vec3& other) const
        double operator[](const size_t i) const
        double& operator[](const size_t i)
        vec3 operator-() const
        double dot(const vec3& b) const
        double dot(const double _x, const double _y, const double _z) const
        vec3 cross(const vec3& b) const
        void cross(const vec3& b, vec3& result) const
        vec3 operator-(const vec3& b) const
        vec3 operator+(const vec3& b) const
        vec3 operator*(const double b) const
        vec3 operator*(const vec3& b) const
        vec3 operator/(const double b) const
        vec3 operator/(const vec3& b) const
        void sub(const vec3& b, vec3& result) const
        void add(const vec3& b, vec3& result) const
        double sqLength() const
        double length() const
        void unitize()
        void normalize()
        vec3 unit() const

    # Vec3Hash
    cdef cppclass Vec3Hash:
        size_t operator()(const vec3& v) const

    # Vec3Equal
    cdef cppclass Vec3Equal:
        bool operator()(const vec3& a, const vec3& b) const

    # vec4
    cdef cppclass vec4:
        double x
        double y
        double z
        double w
        vec4()
        vec4(const double value)
        vec4(const double[4]& arr)
        vec4(const double x, const double y, const double z, const double w)
        vec4(const vec4& other)
        vec3 to_vec3() const
        vec3& to_vec3(vec3& other) const
        size_t size() const
        void set(double x1, double y1, double z1, double w1)
        void set(double[4]& values)
        void set(const vec4& values)
        void set(const vec3& values)
        bool operator==(const vec4& other) const
        double operator[](const size_t i) const
        double& operator[](const size_t i)
        vec4 operator-() const
        double dot(const vec4& b) const
        vec4 cross(const vec4& b) const
        void cross(const vec4& b, vec4& result) const
        void cross(const vec3& b, vec3& result) const
        vec3 cross(const vec3& b) const
        vec4 operator-(const vec4& b) const
        vec4 operator+(const vec4& b) const
        vec4 operator*(const double b) const
        vec4 operator*(const vec4& b) const
        vec4 operator/(const double b) const
        vec4 operator/(const vec4& b) const
        void sub(const vec4& b, vec4& result) const
        void add(const vec4& b, vec4& result) const
        double sqLength() const
        double length() const
        void unitize()
        void normalize()
        vec4 unit() const

    # Function prototypes
    int quantize(double value, int decimals)
    string format_vec3(const vec3& v)
    string format_vec3vec(vector[vec3]& vecs)
    double dot_product(const vec3& a, const vec3& b)
    double mag2(const vec3& a)
    double dot(const vec3& a, const vec3& b)
    double dist(const vec3& a, const vec3& b)
    double dist2(const vec3& a, const vec3& b)
    vec3 operator*(const double v, const vec3& x0)
    void assign(const vec3& a, double& a0, double& a1, double& a2)
