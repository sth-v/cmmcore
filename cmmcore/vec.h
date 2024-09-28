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

#include <cstddef>
#include <iomanip>
#include <random>
#include <unordered_map>
#include <sstream>
#include <memory>
#ifdef __APPLE__
  #include <sys/_types/_ssize_t.h>
#else
#include <cstddef>
#endif
#ifdef __ARM_NEON
#include <arm_neon.h>
#elif defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif
#include <cmath>
namespace cmmcore {
#define CMMCORE_DECIMALS 8
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
  [[nodiscard]] double dot( const vec3& b) const {
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
  [[nodiscard]] double dot( const double _x,const double _y,const double _z) const {

    return x * _x + y * _y + z * _z;
  }
  [[nodiscard]] vec3 cross( const vec3& b) const {
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


  [[nodiscard]] double sqLength() const {
    #ifdef ___ARM_NEON //__ARM_NEON
    float32x4_t v =  vld1q_f32(&x);
    double sq =vaddvq_f32(vmulq_f32(v,v));
    #else
    double sq=x*x+y*y+z*z;
    #endif
    return sq;
  }
  [[nodiscard]] double length() const {
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
  [[nodiscard]] vec3 unit() const {
    vec3 result(x,y,z);
    result.unitize();
    return result;

    }

};

  // Quantize to CMMCORE_DECIMALS decimal places
 inline int quantize(double value, int decimals) {
    return static_cast<int>(value *  decimals);


  }
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
inline std::string format_vec3vec(std::vector<vec3> &vecs) {
  std::string repr("[");
  for (auto& v : vecs) {
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
template<class T>
inline void assign(const std::array<T, 3> &a, T &a0, T &a1, T &a2)
{


  a0=a[0]; a1=a[1]; a2=a[2];
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

  vec4(double value) : x(value), y(value), z(value), w(1) {}

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
  double operator[](const size_t i)const{
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
  [[nodiscard]] double dot( const vec4& b) const {
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


  [[nodiscard]] double sqLength() const {

    double sq=x*x+y*y+z*z+w*w;
    return sq;
  }
  [[nodiscard]] double length() const {
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
    } else {
      x/=l;
      y/=l;
      z/=l;
      w/=l;
    }
  }
  void normalize() {
    unitize();
  }
  [[nodiscard]] vec4 unit() const {
    vec4 result(x,y,z,w);
    result.unitize();
    return result;

    }

};

}
#endif //VEC_H
