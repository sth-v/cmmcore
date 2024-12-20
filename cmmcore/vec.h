/*
Copyright (c) 2024 Andrew Astakhov <aa@contextmachine.ru>. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#ifndef VEC_H
#define VEC_H
#include <stdexcept>
#include <limits>
#include <iostream>
#include <array>
#include <math.h>
#include <cstddef>
#include <iomanip>
#include <random>
#include <unordered_map>
#include <sstream>
#include <memory>
#include <cstddef>
#ifdef CYTHON_ABI
#include "cm_limits.h"
#include "cm_hash.h"
#else
#include "cmmcore/cm_limits.h"
#include "cmmcore/cm_hash.h"
#endif
#ifdef __ARM_NEON
#include <arm_neon.h>
#elif defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif
#include <cmath>
namespace cmmcore {

  // SIMD-friendly 3D vector
  /**
   * @brief Represents a 3D vector with support for SIMD alignment.
   *
   * Vec3 is a basic structure representing a 3D vector with SIMD-friendly 16-byte
   * alignment. It supports default construction, initialization with a scalar
   * value, and explicit initialization with individual x, y, z components.
   */
  struct Vectors3 {
    size_t size;
    double* x;
    double* y;
    double* z;

  };
#define CMMCORE_VEC_PROJECTION(self,other) ((self).dot(other) / (other).dot(other))
#define CMMCORE_VEC_PROJECT(self,other) other*((self).dot(other) / (other).dot(other))
#define CMMCORE_VEC_COLLINEAR(self,other) std::abs((self).unit().dot((other).unit()))==1.
#define CMMCORE_VEC_PARALLEL(self,other) (self).unit().dot((other).unit())==1.
#define CMMCORE_VEC_PERP(self,other) (self).unit().dot((other).unit())==0
  
  struct vec2 {
    double x,y;
    vec2() =default;
    vec2(const double value) : x(value), y(value) {}
    vec2(const double v1,const double v2 ) : x(v1), y(v2) {}
    explicit vec2(const std::array<double,2>& arr) : x(arr[0]), y(arr[1]) {}
    vec2(const vec2& from, const vec2& to) : x(to.x - from.x), y(to.y - from.y) {}
    size_t size() const {
      return 2;
    }
    void set(const double v1) {
      x = v1;
      y=v1;
    }
    void set(const double v1,const double v2) {
      x = v1;
      y=v2;
    }
    void set(const vec2 v1) {
      x = v1.x;
      y=v1.y;
    }
    void set(const std::array<double, 2>& arr) {
      x = arr[0];
      y=arr[1];
    }
    double distance(const vec2& other) const
    {
      return (*this-other).length();
    }
    bool operator==(const vec2& other) const {

      return std::abs(x - other.x) < std::numeric_limits<double>::epsilon() &&
             std::abs(y - other.y) < std::numeric_limits<double>::epsilon();

    }
    vec2 operator-() const {
      return {-x, -y};
    }
    vec2 operator-(const vec2& other) const {
      return {x-other.x, y-other.y};
    }
    vec2 operator*(const vec2& other) const {
      return {x*other.x, y*other.y};
    }
    vec2 operator*(const double& val) const {
      return {x*val, y*val};
    }
    void operator*=(const double& val) {
      x*=val;
      y*=val;
    }
    void operator+=(const vec2& other) {
      x+=other.x;
      y+=other.y;
    }
    void operator/=(const double& val) {
      if (val == 0) {
        throw std::invalid_argument("Division by zero");
      }
      x/=val;
      y/=val;
    }

    void operator-=(const vec2& other) {
      x-=other.x;
      y-=other.y;
    }
    double& operator[](const size_t index) {
      switch (index) {
        case 0:
          return x;
        case 1:
          return y;
        default:
          throw std::out_of_range("index out of range");
      }
    }
    double operator[](const size_t index) const {
      switch (index) {
        case 0:
          return x;
        case 1:
          return y;
        default:
          throw std::out_of_range("index out of range");
      }
    }
    double dot(const vec2& other) const {
      return x * other.x + y * other.y;
    }
    double distanceSq(const vec2& other) const
    {
      return (*this-other).sqLength();
    }
    double sqLength() const {
      return dot(*this);
    }
    double length() const {
      return std::sqrt(dot(*this));
    }
    vec2 project(const vec2& other) const {
      return CMMCORE_VEC_PROJECT(*this,other);
    }
    double projection(const vec2& other) const {
      return CMMCORE_VEC_PROJECTION(*this,other);
    }
    void unitize() {
      double l = length();
      if (l==0) {
        throw std::invalid_argument("Division by zero");
      }
      else {
        x/=l;
        y/=l;

      }
    }
    vec2 unit() const{
      double l = length();
      if (l==0) {
        throw std::invalid_argument("Division by zero");
      }


      return {x/l, y/l};


    }
    bool collinear(const vec2& other) const {
      return CMMCORE_VEC_COLLINEAR(*this,other);
    }
    void operator*=(const double val) {
      x*=val;
      y*=val;

    }
    friend std::ostream&  operator<<(std::ostream& os, const vec2& obj);

  };




inline double cross(const vec2& a, const vec2& b) {
    return a.x * b.y - a.y * b.x;
  }

  inline double dot(const vec2& a, const vec2& b) {
    return a.dot(b);
  }
struct vec3 {
  double x, y, z;

  /// @brief Default constructor, initializes all components to 0.
#ifdef __ARM_NEON

#endif

  vec3() : x(0), y(0), z(0) {}

  /**
   * @brief Initializes all components to the same value.
   *
   * @param value The value to initialize x, y, and z.
   */
  vec3(double value) : x(value), y(value), z(value) {}
  vec3(const std::array<double,3>& arr) : x(arr[0]), y(arr[1]), z(arr[2]) {}
  /**
   * @brief Initializes the vector with individual x, y, z components.
   *
   * @param x The x component.
   * @param y The y component.
   * @param z The z component.
   */
  vec3(double x, double y, double z) : x(x), y(y), z(z) {}
  size_t size() const {
    return 3;
  }
  void set(double x1, double y1, double z1) {
    x=x1;
    y=y1;
    z=z1;
  }
  void set(std::array<double, 3>& values) {
    set(values[0], values[1], values[2]);
  }
  void set(const vec3& values) {
    x=values.x;
    y=values.y;
    z=values.z;

  }
  double distance(const vec3& other) const
  {
    return (*this-other).length();
  }
  void operator*=(const double val) {
    x*=val;
    y*=val;
    z*=val;
  }
  void operator/=(const double val) {
    x/=val;
    y/=val;
    z/=val;
  }
  bool operator==(const vec3& other) const {

    return std::abs(x - other.x) < std::numeric_limits<double>::epsilon() &&
           std::abs(y - other.y) < std::numeric_limits<double>::epsilon()&&
           std::abs(z - other.z) < std::numeric_limits<double>::epsilon();
  }
  /**
 * @brief Access vector components by index.
 *
 * @param i The index (0 for x, 1 for y, 2 for z).
 * @return The value of the specified component.
 * @throws std::runtime_error if the index is out of bounds (i < 0 or i > 2).
 */
  double operator[](const size_t i) const{
    switch (i) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      default:
        throw std::invalid_argument("Index must be 0-2");
    }
  }
   double& operator[](const size_t i){
    switch (i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      throw std::invalid_argument("Index must be 0-2");
    }
  }
  vec3 operator-() const {
    return {-x, -y, -z};
  }
  void operator-=( const vec3& b) {
   x-=b.x;
    y-=b.y;
    z-=b.z;
  }
  void operator+=( const vec3& b) {
    x+=b.x;
    y+=b.y;
    z+=b.z;
  }

  double dot( const vec3& b) const {
#ifdef __ARM_NEON
    float64x2_t va =  vld1q_f64(&x);float64x2_t vb =  vld1q_f64(&b.x);
    float64x2_t _vml=va*vb;
    double xzz=(z*b.z);
    float64x1_t zdot=vld1_f64(&xzz);
    float64x1_t res=vget_low_f64(_vml)+ vget_high_f64(_vml)+zdot;




    return  vget_lane_f64(res,0);
#else

    return x * b.x + y * b.y + z * b.z;
#endif
  }
  bool operator==(const vec2& other) const{
    return other.x==x && other.y==y;
  }
  double dot( const double _x,const double _y,const double _z) const {

    return x * _x + y * _y + z * _z;
  }
  vec3 cross( const vec3& b) const {
    return {
      y * b.z - z * b.y,
      z * b.x - x * b.z,
      x * b.y - y * b.x
  };
  }
  void cross( const vec3& b, vec3& result ) const {
    result.x = y * b.z - z * b.y;
    result.y=  z * b.x - x * b.z;
    result.z=  x * b.y - y * b.x;
  }

  vec3 operator-(const vec3& b) const {
    return {x - b.x, y - b.y, z - b.z};
  }
  vec3 operator+(const vec3& b) const {
    return {x + b.x, y + b.y, z +b.z};
  }
  vec3 operator*(const double b) const {
    return {x * b, y * b, z *b};
  }
  vec3 operator*(const vec3& b) const {
    return {x * b.x, y * b.y, z *b.z};
  }
  vec3 operator/(const double b) const {
    return {x / b, y / b, z /b};
  }
  vec3 operator/(const  vec3& b) const {
    return {x / b.x, y / b.y, z /b.z};
  }
  void sub( const vec3& b, vec3& result ) const {
    result.x=x - b.x;
    result.y=y - b.y;
    result.z=z - b.z;
  }
  void add( const vec3& b, vec3& result ) const {
    result.x=x + b.x;
    result.y=y + b.y;
    result.z=z + b.z;
  }

  double distanceSq(const vec3& other) const
  {
    return (*this-other).sqLength();
  }
  double sqLength() const {
    #ifdef ___ARM_NEON //__ARM_NEON
    float32x4_t v =  vld1q_f32(&x);
    double sq =vaddvq_f32(vmulq_f32(v,v));
    #else
    double sq=x*x+y*y+z*z;
    #endif
    return sq;
  }
  double length() const {
    #ifdef ___ARM_NEON //__ARM_NEON
    float32x4_t v =  vld1q_f32(&x);
    double sq =vaddvq_f32(vmulq_f32(v,v));
    #else
    double sq=x*x+y*y+z*z;
    #endif
    return sqrtf(sq);
  }
  void unitize() {
    if (const double l = length(); l <= std::numeric_limits<double>::epsilon() ||
        l - 1 <= std::numeric_limits<double>::epsilon()) {
    } else {
      x/=l;
      y/=l;
      z/=l;
    }
  }
  void normalize() {
    unitize();
  }
  vec3 unit() const {
    vec3 result(x,y,z);
    result.unitize();
    return result;

    }
  vec3 project(const vec3& other) const {
    return CMMCORE_VEC_PROJECT(*this,other);
  }
  bool collinear(const vec3& other) const {
    return CMMCORE_VEC_COLLINEAR(*this,other);
  }
  friend std::ostream& operator<<(std::ostream& os, const vec3& obj);

};

  inline  void vectorProjection(const vec3&a, const vec3&b,vec3& result)
  {
    const double x0 = pow(b.x, 2);
    const double x1 = pow(b.y, 2);
    const double x2 = pow(b.z, 2);
    const double x3 = 1.0/(x0 + x1 + x2);
    const double x4 = a.y*x3;
    const double x5 = b.x*b.y;
    const double x6 = a.z*x3;
    const double x7 = b.z*x6;
    const double x8 = a.x*x3;
    result.x = b.x*x7 + x0*x8 + x4*x5;
    result.y = b.y*x7 + x1*x4 + x5*x8;
    result.z = b.x*b.z*x8 + b.y*b.z*x4 + x2*x6;
  }

  inline vec3 vectorProjection(const vec3&a, const vec3&b)
  {
    vec3 result;
    vectorProjection(a,b,result);
    return result;
  }

  struct Vec2Hash {
    size_t operator()(const vec2& v) const {
      int qx = quantize(v.x, CMMCORE_DECIMALS);
      int qy = quantize(v.y,CMMCORE_DECIMALS);
      size_t h1 = std::hash<int>()(qx);
      size_t h2 = std::hash<int>()(qy);
      return h1 ^ (h2 << 1) ;
    }


  };

  struct Vec2Equal {
    bool operator()(const vec2& a, const vec2& b) const {
      return(
          quantize(a.x,CMMCORE_DECIMALS) == quantize(b.x,CMMCORE_DECIMALS) &&
          quantize(a.y,CMMCORE_DECIMALS) == quantize(b.y,CMMCORE_DECIMALS) );

    }

  };

struct Vec3Hash {
  size_t operator()(const vec3& v) const {
    int qx = quantize(v.x, CMMCORE_DECIMALS);
    int qy = quantize(v.y,CMMCORE_DECIMALS);
    int qz = quantize(v.z,CMMCORE_DECIMALS);
    size_t h1 = std::hash<int>()(qx);
    size_t h2 = std::hash<int>()(qy);
    size_t h3 = std::hash<int>()(qz);
    return h1 ^ (h2 << 1) ^ (h3 << 2);
  }


};
struct Vec3Equal {
  bool operator()(const vec3& a, const vec3& b) const {
    return
        quantize(a.x,CMMCORE_DECIMALS) == quantize(b.x,CMMCORE_DECIMALS) &&
        quantize(a.y,CMMCORE_DECIMALS) == quantize(b.y,CMMCORE_DECIMALS) &&
        quantize(a.z,CMMCORE_DECIMALS) == quantize(b.z,CMMCORE_DECIMALS);
  }

};


using Vec3UnorderedMap=std::unordered_map<vec3, int, Vec3Hash, Vec3Equal>;

inline std::string format_vec3(const vec3 &v) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2);
  oss << "[" << v.x << "," << v.y << "," << v.z << "]";
  return oss.str();
}
inline std::string format_vec3vec(const std::vector<vec3> &vecs) {
  std::string repr("[");
  for (const auto& v : vecs) {
    repr+=format_vec3(v);
    repr+=",";
  }
  repr+="]";
  return repr;
}
#ifdef ___ARM_NEON
inline double dot_product( const vec3& a, const vec3& b ) {float32x2_t va =  vld1_f32(&a.x);float32x2_t vb =  vld1_f32(&b.x);float32_t zdot=a.z*b.z;return vpadds_f32(vmul_f32(va,vb))+zdot;}
inline void dot_product(const vec3& v, const Vectors3& vecs, double* result) {
  // Load the Vec3 elements into NEON registers
  float32x4_t v_x = vdupq_n_f32(v.x);
  float32x4_t v_y = vdupq_n_f32(v.y);
  float32x4_t v_z = vdupq_n_f32(v.z);

  // Process in chunks of 4 using NEON
  size_t i = 0;
  for (; i + 4 <= vecs.size; i += 4) {
    // Load 4 elements from vecs.x, vecs.y, vecs.z into NEON registers
    float32x4_t vec_x = vld1q_f32(vecs.x + i);
    float32x4_t vec_y = vld1q_f32(vecs.y + i);
    float32x4_t vec_z = vld1q_f32(vecs.z + i);

    // Compute the dot product for 4 vectors at a time
    float32x4_t dot_x = vmulq_f32(v_x, vec_x);
    float32x4_t dot_y = vmulq_f32(v_y, vec_y);
    float32x4_t dot_z = vmulq_f32(v_z, vec_z);

    // Sum the partial products to get the dot products
    float32x4_t dot_product = vaddq_f32(dot_x, vaddq_f32(dot_y, dot_z));

    // Store the result in the result array
    vst1q_f32(result + i, dot_product);
  }

  // Handle the remaining elements that don't fit in a group of 4
  for (; i < vecs.size; ++i) {
    result[i] = v.x * vecs.x[i] + v.y * vecs.y[i] + v.z * vecs.z[i];
  }
}
#elif defined(__x86_64__) || defined(_M_X64)
inline void dot_product(const Vec3& v, const Vectors3& vecs, double* result) {
  // Load the Vec3 elements into SSE registers
  __m128 v_x = _mm_set1_ps(v.x);  // Duplicate v.x across the 4 lanes of a 128-bit register
  __m128 v_y = _mm_set1_ps(v.y);  // Duplicate v.y across the 4 lanes of a 128-bit register
  __m128 v_z = _mm_set1_ps(v.z);  // Duplicate v.z across the 4 lanes of a 128-bit register

  // Process in chunks of 4 using SSE
  size_t i = 0;
  for (; i + 4 <= vecs.size; i += 4) {
    // Load 4 elements from vecs.x, vecs.y, vecs.z into SSE registers
    __m128 vec_x = _mm_loadu_ps(vecs.x + i);  // Unaligned load
    __m128 vec_y = _mm_loadu_ps(vecs.y + i);  // Unaligned load
    __m128 vec_z = _mm_loadu_ps(vecs.z + i);  // Unaligned load

    // Compute the dot product for 4 vectors at a time
    __m128 dot_x = _mm_mul_ps(v_x, vec_x);  // v.x * vecs.x
    __m128 dot_y = _mm_mul_ps(v_y, vec_y);  // v.y * vecs.y
    __m128 dot_z = _mm_mul_ps(v_z, vec_z);  // v.z * vecs.z

    // Sum the partial products to get the dot products
    __m128 dot_product = _mm_add_ps(dot_x, _mm_add_ps(dot_y, dot_z));  // dot_x + dot_y + dot_z

    // Store the result in the result array
    _mm_storeu_ps(result + i, dot_product);  // Unaligned store
  }

  // Handle the remaining elements that don't fit in a group of 4
  for (; i < vecs.size; ++i) {
    result[i] = v.x * vecs.x[i] + v.y * vecs.y[i] + v.z * vecs.z[i];
  }
}

#else
inline void dot_product(const vec3& v, const Vectors3& vecs, double* result) {
  for (int i=0; i  <= vecs.size; i ++) {
    result[i]=v.dot(vecs.x[i],vecs.y[i],vecs.z[i]);
  }
}
#endif
#define sqr(x) ((x)*(x))
double mag2(const vec3 &a)
{
  double l=std::sqrt(a[0]);
  for(unsigned int i=1; i<3; ++i) l+=sqr(a[i]);
  return l;
}

inline double dot(const vec3& a, const vec3& b) {
  return a.dot(b);
}
inline double dist(const vec3& a, const vec3& b) {

  return  std::sqrt((a.x-b.x)+(a.y-b.y)+(a.z-b.z));;
}
inline double dist2(const vec3& a, const vec3& b) {

  return (a-b).sqLength();
}

inline vec3 operator*(const double v,const vec3 &x0) {
  return {x0.x*v,x0.y*v,x0.z*v};
}

inline void assign(const vec3 &a, double &a0, double &a1, double &a2)
{

  a0=a.x; a1=a.y; a2=a.z;
}


inline void cartesian_to_spherical(const vec3& xyz,vec3& rtp ) {



    double XsqPlusYsq = xyz.x *  xyz.x +  xyz.y *  xyz.y;

    rtp.x= std::sqrt(XsqPlusYsq +  xyz.z *  xyz.z);        // Radius r
    rtp.y = std::atan2(std::sqrt(XsqPlusYsq), xyz.z); // Polar angle θ (theta)
    rtp.z = std::atan2(xyz.y, xyz.x);                     // Azimuthal angle φ (phi)
  }

  inline void spherical_to_cartesian(const vec3& rtp ,vec3& xyz) {
  double r = rtp.x;
  double theta = rtp.y; // Polar angle
  double phi = rtp.z;   // Azimuthal angle φ
  xyz.x = r * std::sin(theta) * std::cos(phi);
    xyz.y = r * std::sin(theta) * std::sin(phi);
   xyz.z = r * std::cos(theta);
}

  inline void cartesian_to_spherical(const vec3& xyz,vec2& tp ) {
  tp.x = std::atan2(std::sqrt(xyz.x *  xyz.x +  xyz.y *  xyz.y), xyz.z); // Polar angle θ (theta)
  tp.y = std::atan2(xyz.y, xyz.x);
}

  inline void spherical_to_cartesian(const vec2& tp ,vec3& xyz) {
  double r = 1;
  double theta = tp.x; // Polar angle
  double phi = tp.y;   // Azimuthal angle φ
  xyz.x = r * std::sin(theta) * std::cos(phi);
  xyz.y = r * std::sin(theta) * std::sin(phi);
  xyz.z = r * std::cos(theta);
}
  inline void cartesian_to_spherical(const std::vector<vec3>& xyz,std::vector<vec2>& tp ) {
  tp.resize(xyz.size());
  for (size_t i = 0; i < xyz.size(); ++i) {
    tp[i].x = std::atan2(std::sqrt(xyz[i].x *  xyz[i].x +  xyz[i].y *  xyz[i].y), xyz[i].z); // Polar angle θ (theta)
    tp[i].y = std::atan2(xyz[i].y, xyz[i].x);
  }

}


struct vec4 {
  double x, y, z, w;

  /// @brief Default constructor, initializes all components to 0.
#ifdef __ARM_NEON

#endif

  vec4() : x(0), y(0), z(0), w(0) {}

  /**
   * @brief Initializes all components to the same value.
   *
   * @param value The value to initialize x, y, and z.
   */

  vec4(double value) : x(value), y(value), z(value), w(value) {}

  /**
   * @brief Initializes the vector with individual x, y, z components.
   *
   * @param x The x component.
   * @param y The y component.
   * @param z The z component.
   */
  vec4(const std::array<double,4>& arr) : x(arr[0]), y(arr[1]), z(arr[2]),w(arr[3]) {}

  vec4(double x, double y, double z,double w) : x(x), y(y), z(z),w(w) {}
  vec4(const vec4& other) : x(other.x), y(other.y), z(other.z),w(other.w) {}
  vec3 to_vec3() const {


    return vec3(x/w, y/w, z/w);
  }
  vec3& to_vec3(vec3& other) const {

      other.x=x/w;
      other.y=y/w;
      other.z=z/w;
      return other;

  }

  size_t size() const {
    return 4;
  }
  void set(double x1, double y1, double z1,double w1) {
    x=x1;
    y=y1;
    z=z1;
    w=w1;
  }
  void set(const std::array<double, 4>& values) {
    set(values[0], values[1], values[2], values[3]);
  }
  void set(const vec4& values) {
    x=values.x;
    y=values.y;
    z=values.z;
    w=values.w;

  }
  void set(const vec3& values) {
    x=values.x;
    y=values.y;
    z=values.z;
    w=1.;

  }
  vec4 project(const vec4& other) const {
    return CMMCORE_VEC_PROJECT(*this,other);
  }
  double distance(const vec4& other) const
  {
    return (*this-other).length();
  }
  double distanceSq(const vec4& other) const
  {
    return (*this-other).sqLength();
  }
  bool operator==(const vec4& other) const {

    return std::abs(x - other.x) < std::numeric_limits<double>::epsilon() &&
           std::abs(y - other.y) < std::numeric_limits<double>::epsilon()&&
           std::abs(z - other.z) < std::numeric_limits<double>::epsilon();
  }
  /**
 * @brief Access vector components by index.
 *
 * @param i The index (0 for x, 1 for y, 2 for z).
 * @return The value of the specified component.
 * @throws std::runtime_error if the index is out of bounds (i < 0 or i > 2).
 */
  double operator[](const size_t i) const{
    switch (i) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      case 3:
        return w;
      default:
        throw std::invalid_argument("Index must be 0-3");
    }
  }
  double& operator[](const size_t i){
    switch (i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    case 3:
      return w;
    default:
      throw std::invalid_argument("Index must be 0-3");
    }
  }
  vec4 operator-() const {
    return {-x, -y, -z, -w};
  }
  void operator-=( const vec4& b) {
   x-=b.x;
    y-=b.y;
    z-=b.z;
    w-=b.w;
  }
  void operator+=( const vec4& b) {
    x+=b.x;
    y+=b.y;
    z+=b.z;
    w+=b.w;
  }
  double dot( const vec4& b) const {
#ifdef __ARM_NEON
    float64x2_t va =  vld1q_f64(&x);float64x2_t vb =  vld1q_f64(&b.x);
    float64x2_t _vml=va*vb;

    float64x2_t vaa =  vld1q_f64(&z);float64x2_t vb1 =  vld1q_f64(&b.z);
    float64x2_t _vml2=vaa*vb1;


    float64x1_t res=vget_low_f64(_vml)+ vget_high_f64(_vml)+vget_low_f64(_vml2)+ vget_high_f64(_vml2);







    //float32x4_t _vmul=vmul_f64(va,vb);
    // First, add the two pairs of elements


    // Then add the remaining two elements

    // Extract the final scalar result
    return    vget_lane_f64(res,0);

#else

    return x * b.x + y * b.y + z * b.z;
#endif
  }

  vec4 cross( const vec4& b) const {
    return {
      y * b.z - z * b.y,
      z * b.x - x * b.z,
      x * b.y - y * b.x,
      1.
  };
  }
  void cross( const vec4& b, vec4& result ) const {
    auto v1=to_vec3();
    auto v2=b.to_vec3();
    result.x = v1.y * v2.z -  v2.z *  v2.y;
    result.y=  v1.z * v2.x  - v2.x * v2.z;
    result.z=  v1.x * v2.y  - v2.y * v2.x;
    result.w=1.;
  }
  void cross( const vec3& b, vec3& result ) const {
    auto _x=x/w;
    auto _y=y/w;
    auto _z=z/w;
       
     
     
    result.x = _y * b.z - _z * b.y;
    result.y=  _z * b.x - _x * b.z;
    result.z=  _x * b.y - _y * b.x;
  }
  vec3 cross( const vec3& b) const {

    auto _x=x/w;
    auto _y=y/w;
    auto _z=z/w;
    return {_y * b.z - _z * b.y,_z * b.x - _x * b.z,_x * b.y - _y * b.x};
  }

  vec4 operator-(const vec4& b) const {
    return {x - b.x, y - b.y, z - b.z, w - b.w};
  }
  vec4 operator+(const vec4& b) const {
    vec4 result;
#ifdef ___ARM_NEON
    float32x4_t va =  vld1q_f32(&x);
    float32x4_t vb =  vld1q_f32(&b.x);
    float32x4_t vn = vpaddq_f32(va,vb);
    vst1q_f32(&result.x, vn);
    return result;
#else

    return {x + b.x, y + b.y, z +b.z,w+b.w};
#endif

  }
  vec4 operator*(const double b) const {
    return {x * b, y * b, z *b, w*b};
  }
  vec4 operator*(const vec4& b) const {
    return {x * b.x, y * b.y, z *b.z,w*b.w};
  }
  vec4 operator/(const double b) const {
    return {x / b, y / b, z /b, w /b};
  }
  vec4 operator/(const  vec4& b) const {
    return {x / b.x, y / b.y, z /b.z, w /b.w};
  }
  void sub( const vec4& b, vec4& result ) const {
    result.x=x - b.x;
    result.y=y - b.y;
    result.z=z - b.z;
    result.w=w - b.w;
  }
  void add( const vec4& b, vec4& result ) const {
#ifdef __ARM_NEON
    float64x2x2_t va = vld2q_f64(&x);
    float64x2x2_t vb =  vld2q_f64(&b.x);

    float64x2x2_t sm={vaddq_f64(va.val[0],vb.val[0]),vaddq_f64(va.val[1],vb.val[1])};
    vst2q_f64(&result.x, sm);

#else

    result.x=x + b.x;
    result.y=y + b.y;
    result.z=z + b.z;
#endif

  }


  double sqLength() const {

    double sq=x*x+y*y+z*z+w*w;
    return sq;
  }
  double length() const {
    #ifdef ___ARM_NEON
    auto sq =sqLength();

    auto res=std::sqrt(vld1_f64(&sq));

    return res;
    #else
    double sq=x*x+y*y+z*z+w*w;
    return sqrtf(sq);
    #endif

  }
  void unitize() {
    if (const double l = length(); l <= std::numeric_limits<double>::epsilon() ||
        l - 1 <= std::numeric_limits<double>::epsilon()) {
    } else
    {
      x/=(w*l);
      y/=(w*l);
      z/=(w*l);

      w=1.;

    }
  }
  void normalize() {
    unitize();
  }
  vec4 unit() const {
    vec4 result(x,y,z,w);
    result.unitize();
    return result;

    }
  bool collinear(const vec4& other) const {
    return CMMCORE_VEC_COLLINEAR(*this,other);
  }
  friend std::ostream& operator<<(std::ostream& os, const vec4& obj);
};
  inline void cross(const vec3& a, const vec3& b, vec4& result) {
    result.x=a.y*b.z - a.z*b.y;
    result.y=a.z*b.x - a.x*b.z;
    result.z=a.x*b.y - a.y*b.x;
    result.w=1.0;
  }



  std::ostream& operator<<(std::ostream& os, const vec2& obj) {
    os << "[" << obj.x << "," << obj.y << "]";
    return os;
  }
  std::ostream& operator<<(std::ostream& os, const vec3& obj) {
    os << "[" << obj.x << "," << obj.y  << obj.z <<"]";
    return os;
  }


  std::ostream& operator<<(std::ostream& os, const vec4& obj) {
    os << "[" << obj.x << "," << obj.y  << obj.z << obj.w << "]";
    return os;
  }

}
#endif //VEC_H
