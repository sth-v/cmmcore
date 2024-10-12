//
// Created by Andrew Astakhov on 11.10.24.
//

#ifndef AABB_H
#define AABB_H
#ifdef CYTHON_ABI
#include "vec.h"

#else
#include "cmmcore/vec.h"

#endif
namespace cmmcore{
  struct alignas(16) AABB {
    vec3 min, max;

    AABB()
      : min(std::numeric_limits<double>::max()),
        max(std::numeric_limits<double>::lowest()) {
    }

    AABB(const vec3 &min, const vec3 &max) : min(min), max(max) {
    }

    AABB(const std::vector<vec3>& pts) {
      min.set(pts[0]);
      max.set(pts[0]);
      for (size_t i=1; i<pts.size(); i++) {
        expand(pts[i]);
      }
    }
    bool infinity() const {
      return (max.x == max.y == max.z == std::numeric_limits<double>::min() && min.x == min.y == min.z ==
              std::numeric_limits<double>::max());
    }

    /**
     * @brief Checks if this AABB intersects with another.
     *
     * @param other The other AABB to check for intersection.
     * @return True if the AABBs intersect, false otherwise.
     */
    bool intersects(const AABB &other) const {
#ifdef __ARM_NEON
      // ARM NEON double-precision support is available in ARMv8-A and later.
      // We use float64x2_t for vectors of two doubles.

      // Load x and y components
      float64x2_t min1_xy = vld1q_f64(&min.x);
      float64x2_t max1_xy = vld1q_f64(&max.x);
      float64x2_t min2_xy = vld1q_f64(&other.min.x);
      float64x2_t max2_xy = vld1q_f64(&other.max.x);

      // Compare min1 <= max2 and max1 >= min2 for x and y
      uint64x2_t cmp1_xy = vcleq_f64(min1_xy, max2_xy);
      uint64x2_t cmp2_xy = vcgeq_f64(max1_xy, min2_xy);
      uint64x2_t result_xy = vandq_u64(cmp1_xy, cmp2_xy);

      // Extract comparison results
      uint64_t cmp_results[2];
      vst1q_u64(cmp_results, result_xy);
      bool intersects_xy = (cmp_results[0] != 0) && (cmp_results[1] != 0);

      // Process z component using scalar operations
      bool intersects_z = min.z <= other.max.z && max.z >= other.min.z;

      return intersects_xy && intersects_z;

#elif defined(__x86_64__) || defined(_M_X64)
    // Use SSE2 intrinsics for double-precision operations

    // Load x and y components
    __m128d min1_xy = _mm_load_pd(&min.x);
    __m128d max1_xy = _mm_load_pd(&max.x);
    __m128d min2_xy = _mm_load_pd(&other.min.x);
    __m128d max2_xy = _mm_load_pd(&other.max.x);

    // Compare min1 <= max2 and max1 >= min2 for x and y
    __m128d cmp1_xy = _mm_cmple_pd(min1_xy, max2_xy);
    __m128d cmp2_xy = _mm_cmpge_pd(max1_xy, min2_xy);
    __m128d result_xy = _mm_and_pd(cmp1_xy, cmp2_xy);

    // Extract comparison results
    int mask_xy = _mm_movemask_pd(result_xy); // Returns a 2-bit mask
    bool intersects_xy = (mask_xy == 0x3);    // Both bits set

    // Process z component using scalar operations
    bool intersects_z = min.z <= other.max.z && max.z >= other.min.z;

    return intersects_xy && intersects_z;
#else
    // Scalar fallback
    return min.x <= other.max.x && max.x >= other.min.x &&
           min.y <= other.max.y && max.y >= other.min.y &&
           min.z <= other.max.z && max.z >= other.min.z;
#endif
    }

    bool intersection(const AABB &other, AABB &result) const {
#ifdef __ARM_NEON
      // Load x and y components
      float64x2_t min1_xy = vld1q_f64(&min.x);
      float64x2_t max1_xy = vld1q_f64(&max.x);
      float64x2_t min2_xy = vld1q_f64(&other.min.x);
      float64x2_t max2_xy = vld1q_f64(&other.max.x);

      // Compare min1 <= max2 and max1 >= min2 for x and y
      uint64x2_t cmp1_xy = vcleq_f64(min1_xy, max2_xy);
      uint64x2_t cmp2_xy = vcgeq_f64(max1_xy, min2_xy);
      uint64x2_t result_xy = vandq_u64(cmp1_xy, cmp2_xy);

      // Extract comparison results
      uint64_t cmp_results[2];
      vst1q_u64(cmp_results, result_xy);
      bool intersects_xy = (cmp_results[0] != 0) && (cmp_results[1] != 0);

      // Process z component using scalar operations
      bool intersects_z = min.z <= other.max.z && max.z >= other.min.z;

      if (intersects_xy && intersects_z) {
        // Compute new min and max for x and y
        float64x2_t new_min_xy = vmaxq_f64(min1_xy, min2_xy);
        float64x2_t new_max_xy = vminq_f64(max1_xy, max2_xy);

        // Store results
        vst1q_f64(&result.min.x, new_min_xy);
        vst1q_f64(&result.max.x, new_max_xy);

        // Compute new min and max for z
        result.min.z = std::max(min.z, other.min.z);
        result.max.z = std::min(max.z, other.max.z);

        return true;
      }
      return false;

#elif defined(__x86_64__) || defined(_M_X64)
    // Load x and y components
    __m128d min1_xy = _mm_load_pd(&min.x);
    __m128d max1_xy = _mm_load_pd(&max.x);
    __m128d min2_xy = _mm_load_pd(&other.min.x);
    __m128d max2_xy = _mm_load_pd(&other.max.x);

    // Compare min1 <= max2 and max1 >= min2 for x and y
    __m128d cmp1_xy = _mm_cmple_pd(min1_xy, max2_xy);
    __m128d cmp2_xy = _mm_cmpge_pd(max1_xy, min2_xy);
    __m128d result_xy = _mm_and_pd(cmp1_xy, cmp2_xy);

    // Extract comparison results
    int mask_xy = _mm_movemask_pd(result_xy); // Returns a 2-bit mask
    bool intersects_xy = (mask_xy == 0x3);    // Both bits set

    // Process z component using scalar operations
    bool intersects_z = min.z <= other.max.z && max.z >= other.min.z;

    if (intersects_xy && intersects_z) {
      // Compute new min and max for x and y
      __m128d new_min_xy = _mm_max_pd(min1_xy, min2_xy);
      __m128d new_max_xy = _mm_min_pd(max1_xy, max2_xy);

      // Store results
      _mm_store_pd(&result.min.x, new_min_xy);
      _mm_store_pd(&result.max.x, new_max_xy);

      // Compute new min and max for z
      result.min.z = std::max(min.z, other.min.z);
      result.max.z = std::min(max.z, other.max.z);

      return true;
    }
    return false;
#else
    // Scalar fallback
    if (min.x <= other.max.x && max.x >= other.min.x &&
        min.y <= other.max.y && max.y >= other.min.y &&
        min.z <= other.max.z && max.z >= other.min.z) {
      result.min = vec3(std::max(min.x, other.min.x),
                        std::max(min.y, other.min.y),
                        std::max(min.z, other.min.z));
      result.max = vec3(std::min(max.x, other.max.x),
                        std::min(max.y, other.max.y),
                        std::min(max.z, other.max.z));
      return true;
    }
    return false;
#endif
    }

    /**
     * @brief Merges this AABB with another AABB and returns the result.
     *
     * @param other The other AABB to merge with.
     * @return A new AABB that encompasses both this and the other AABB.
     */
    AABB merge(const AABB &other) const {
#ifdef __ARM_NEON
      // Load x and y components
      float64x2_t min1_xy = vld1q_f64(&min.x);
      float64x2_t max1_xy = vld1q_f64(&max.x);
      float64x2_t min2_xy = vld1q_f64(&other.min.x);
      float64x2_t max2_xy = vld1q_f64(&other.max.x);

      // Compute new min and max for x and y
      float64x2_t new_min_xy = vminq_f64(min1_xy, min2_xy);
      float64x2_t new_max_xy = vmaxq_f64(max1_xy, max2_xy);

      AABB result;
      vst1q_f64(&result.min.x, new_min_xy);
      vst1q_f64(&result.max.x, new_max_xy);

      // Compute new min and max for z
      result.min.z = std::min(min.z, other.min.z);
      result.max.z = std::max(max.z, other.max.z);

      return result;

#elif defined(__x86_64__) || defined(_M_X64)
    // Load x and y components
    __m128d min1_xy = _mm_load_pd(&min.x);
    __m128d max1_xy = _mm_load_pd(&max.x);
    __m128d min2_xy = _mm_load_pd(&other.min.x);
    __m128d max2_xy = _mm_load_pd(&other.max.x);

    // Compute new min and max for x and y
    __m128d new_min_xy = _mm_min_pd(min1_xy, min2_xy);
    __m128d new_max_xy = _mm_max_pd(max1_xy, max2_xy);

    AABB result;
    _mm_store_pd(&result.min.x, new_min_xy);
    _mm_store_pd(&result.max.x, new_max_xy);

    // Compute new min and max for z
    result.min.z = std::min(min.z, other.min.z);
    result.max.z = std::max(max.z, other.max.z);

    return result;
#else
    // Scalar fallback
    return AABB(vec3(std::min(min.x, other.min.x), std::min(min.y, other.min.y),
                     std::min(min.z, other.min.z)),
                vec3(std::max(max.x, other.max.x), std::max(max.y, other.max.y),
                     std::max(max.z, other.max.z)));
#endif
    }

    /**
     * @brief Computes the volume of the AABB.
     *
     * @return The volume of the AABB.
     */
    double volume() const {
#ifdef __ARM_NEON
      // Load x and y components
      float64x2_t min_vec = vld1q_f64(&min.x);
      float64x2_t max_vec = vld1q_f64(&max.x);
      float64x2_t diff = vsubq_f64(max_vec, min_vec);

      // Multiply x and y differences
      double diff_array[2];
      vst1q_f64(diff_array, diff);
      double area_xy = diff_array[0] * diff_array[1];

      // Compute z difference
      double diff_z = max.z - min.z;

      return area_xy * diff_z;

#elif defined(__x86_64__) || defined(_M_X64)
    // Load x and y components
    __m128d min_vec = _mm_load_pd(&min.x);
    __m128d max_vec = _mm_load_pd(&max.x);
    __m128d diff = _mm_sub_pd(max_vec, min_vec);

    // Multiply x and y differences
    double diff_array[2];
    _mm_storeu_pd(diff_array, diff);
    double area_xy = diff_array[0] * diff_array[1];

    // Compute z difference
    double diff_z = max.z - min.z;

    return area_xy * diff_z;
#else
    // Scalar fallback
    vec3 diff = {max.x - min.x, max.y - min.y, max.z - min.z};
    return diff.x * diff.y * diff.z;
#endif
    }

    void expand(const vec3 &point) {
      min.x = std::min(min.x, point.x);
      min.y = std::min(min.y, point.y);
      min.z = std::min(min.z, point.z);
      max.x = std::max(max.x, point.x);
      max.y = std::max(max.y, point.y);
      max.z = std::max(max.z, point.z);
    }
  };

    class SphericalAABB {
  public:
    vec3 start;
    vec3 end;

    SphericalAABB() = default;

    SphericalAABB(const std::vector<vec3> &pts, const vec3 &center = 0): start(0), end(0) {
      std::vector<vec3> sphere_points1;
      project_to_sphere(pts, center, sphere_points1);
      for (int axis = 0; axis < 3; ++axis) {
        std::vector<vec2> plane_coords1;
        plane_coords1.reserve(sphere_points1.size());

        for (const auto &pt: sphere_points1) {
          if (axis == 0) {
            plane_coords1.emplace_back(pt.y, pt.z);
          } else if (axis == 1) {
            plane_coords1.emplace_back(pt.z, pt.x);
          } else {
            plane_coords1.emplace_back(pt.x, pt.y);
          }
        }
        find_smallest_wedge(plane_coords1, axis);
      }
    }

    bool separable(const SphericalAABB &other) const {
      for (int axis = 0; axis < 3; ++axis) {
        if ((end[axis] < other.start[axis] && other.end[axis] > start[axis]) || (
              other.end[axis] < start[axis] && end[axis] > other.start[axis])) {
          return true;
        }
      }
      return false;
    };

    bool intersects(const SphericalAABB &other) const {
      return !separable(other);
    }

  private:
    static void project_to_sphere(const std::vector<vec3> &points, const vec3 &center, std::vector<vec3> &result) {
      result.resize(points.size());
      for (int i = 0; i < points.size(); ++i) {
        const auto &point = points[i];
        point.sub(center, result[i]);

        result[i].unitize();
      }
    }

    void find_smallest_wedge(const std::vector<vec2> &points_2d, size_t dim) {
      std::vector<double> angles(points_2d.size());
      for (int i = 0; i < points_2d.size(); ++i) {
        angles[i] = std::atan2(points_2d[i].y, points_2d[i].x);
      }
      std::sort(angles.begin(), angles.end());

      std::vector<double> angle_diffs(angles.size());
      std::adjacent_difference(angles.begin(), angles.end(), angle_diffs.begin());
      angle_diffs[0] = angles[0] + 2 * M_PI - angles.back(); // Handle circular difference

      auto max_gap_it = std::max_element(angle_diffs.begin(), angle_diffs.end());
      size_t max_gap_idx = std::distance(angle_diffs.begin(), max_gap_it);

      start[dim] = angles[max_gap_idx];
      end[dim] = angles[(max_gap_idx + 1) % angles.size()];
    }
  };
}
#endif //AABB_H
