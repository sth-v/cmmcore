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

| Platform           | Compiler        | Command |
|--------------------|------------x-----|-------------------------------------------------------------------------------|
| **macOS (M1)**     | `clang++`       | `clang++ -std=c++17 -o bvh_program
main.cpp -march=armv8-a+simd -O3`          | | **Linux (ARM64)**  | `g++` | `g++
-std=c++17 -o bvh_program main.cpp -march=armv8-a+simd -O3`              | |
**Linux (AMD)**    | `g++`           | `g++ -std=c++17 -o bvh_program main.cpp
-march=x86-64-v3 -O3`                 | | **Windows (MSVC)** | `cl.exe`        |
`cl /EHsc /std:c++17 /O2 /arch:AVX2 main.cpp /Fe:bvh_program.exe`             |

### Explanation for **Linux (ARM64)**:
- This uses the same architecture flag as macOS for ARM, targeting `armv8-a`
with SIMD instructions using `NEON`.

### Explanation for **Windows (MSVC)**
- cl: Microsoft C++ compiler.
- /EHsc: Enables C++ exception handling.
- /std:c++17: Specifies C++17 as the standard.
- /O2: Enables optimization for speed (equivalent to -O3 in GCC/Clang).
- /arch:AVX2: Targets AVX2 instructions, common for modern AMD and Intel
processors.
- main.cpp: The name of the C++ source file.
- /Fe:bvh_program.exe: Specifies the output executable name as bvh_program.exe.

*/

#ifndef BVH_HPP
#define BVH_HPP
#include "cmmcore/vec.h"
#include <algorithm>
#include <stack>
#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>
#include <cstddef>
#include <iomanip>
#include <random>
#include <cstring>
#include <sstream>
#ifdef __ARM_NEON
#include <arm_neon.h>
#elif defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>

#endif
namespace cmmcore {
#ifdef __ARM_NEON
#include <arm_neon.h>
#elif defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif
// SIMD-optimized AABB (Axis-Aligned Bounding Box)
/**
 * @brief Represents an axis-aligned bounding box with SIMD optimization.
 *
 * This structure is designed to work efficiently with SIMD intrinsics on ARM
 * and x86 architectures.
 */
struct alignas(16) AABB {
  vec3 min, max;

  AABB()
      : min(std::numeric_limits<double>::max()),
        max(std::numeric_limits<double>::lowest()) {}
  AABB(const vec3 &min, const vec3 &max) : min(min), max(max) {}

  bool infinity() const{
    return (max.x==max.y==max.z==std::numeric_limits<double>::min() &&  min.x==min.y==min.z==std::numeric_limits<double>::max());

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
  void expand(const vec3& point) {
    min.x = std::min(min.x, point.x);
    min.y = std::min(min.y, point.y);
    min.z = std::min(min.z, point.z);
    max.x = std::max(max.x, point.x);
    max.y = std::max(max.y, point.y);
    max.z = std::max(max.z, point.z);
  }
};
// 3D Object Representation
/**
 * @brief Represents a 3D object with an associated bounding box.
 *
 * Each object is assigned a unique ID and has an associated axis-aligned
 * bounding box.
 */
class Object3D {
public:
  AABB bounding_box;
  int id = 0;

  /// @brief Default constructor, initializes the bounding box to an invalid
  /// state and ID to 0.
  Object3D() {
    bounding_box = AABB();
    id = 0;
  }
  /**
   * @brief Constructs an object with a specific bounding box.
   *
   * @param bb The bounding box for the object.
   */
  Object3D(const AABB &bb) : bounding_box(bb) { id = 0; }
  /**
   * @brief Constructs an object with a specific bounding box and ID.
   *
   * @param bb The bounding box for the object.
   * @param id The unique identifier for the object.
   */
  Object3D(const AABB &bb, int id) : bounding_box(bb), id(id) {}
};

// BVH Node Representation
/**
 * @brief Represents a node in a Bounding Volume Hierarchy (BVH).
 *
 * Each node contains a bounding box, pointers to left and right child nodes,
 * and may contain an object.
 */
class BVHNode {
public:
  AABB bounding_box;
  std::unique_ptr<BVHNode> left;
  std::unique_ptr<BVHNode> right;
  Object3D *object;
  BVHNode() : left(nullptr), right(nullptr), object(nullptr) {

  }
  BVHNode(const BVHNode& other) : bounding_box(other.bounding_box), object(nullptr) {
    if (other.left) {
      left = std::make_unique<BVHNode>(*other.left);  // Recursively copy left subtree
    }
    if (other.right) {
      right = std::make_unique<BVHNode>(*other.right);  // Recursively copy right subtree
    }
  }

  explicit BVHNode(const AABB &bb, Object3D *obj = nullptr)
      : bounding_box(bb), object(obj) {}
  void setLeft(std::unique_ptr<BVHNode> &node) { left.swap(node); }
  void setRight(std::unique_ptr<BVHNode> &node) { left.swap(node); }
  // Custom copy assignment operator for deep copy
  BVHNode& operator=(const BVHNode& other) {
    if (this == &other) return *this; // Handle self-assignment

    bounding_box = other.bounding_box;

    object=other.object? object : nullptr;

    if (other.left) {
      left = std::make_unique<BVHNode>(*other.left);  // Recursively copy left subtree
    } else {
      left.reset();  // Ensure left is null if other.left is null
    }
    if (other.right) {
      right = std::make_unique<BVHNode>(*other.right);  // Recursively copy right subtree
    } else {
      right.reset();  // Ensure right is null if other.right is null
    }

    return *this;
  }
};
inline vec3 compute_extents(const AABB &bbox) {
  return vec3{bbox.max.x - bbox.min.x, bbox.max.y - bbox.min.y,
              bbox.max.z - bbox.min.z};
}
// The rest of the functions (split_objects, build_bvh, intersect_bvh,
// traverse_bvh, contains_point) remain the same as in the previous optimized
// version

// Function to split objects along the longest axis
std::pair<std::vector<Object3D *>,
          std::vector<Object3D *>> inline split_objects(std::vector<Object3D *>
                                                            &objects) {
  if (objects.size() <= 1) {
    return {objects, {}};
  }

  AABB total_bounds;
  for (const auto &obj : objects) {
    total_bounds = total_bounds.merge(obj->bounding_box);
  }

  vec3 extents = compute_extents(total_bounds);

  // Check if all extents are zero
  if (extents.x == 0 && extents.y == 0 && extents.z == 0) {
    // Split the objects arbitrarily (e.g., in half)
    size_t mid = objects.size() / 2;
    return {std::vector<Object3D *>(objects.begin(), objects.begin() + mid),
            std::vector<Object3D *>(objects.begin() + mid, objects.end())};
  }

  // Sort axes by extent, descending order
  std::array<int, 3> axes = {0, 1, 2};
  std::sort(axes.begin(), axes.end(),
            [&extents](int a, int b) { return extents[a] > extents[b]; });

  for (int axis : axes) {
    if (extents[axis] > 0) {
      double split_pos =
          (total_bounds.min[axis] + total_bounds.max[axis]) * 0.5f;

      auto it = std::partition(objects.begin(), objects.end(),
                               [axis, split_pos](const Object3D *obj) {
                                 return obj->bounding_box.min[axis] < split_pos;
                               });

      std::vector<Object3D *> left_objects(objects.begin(), it);
      std::vector<Object3D *> right_objects(it, objects.end());

      // Check if we've actually split the objects
      if (!left_objects.empty() && !right_objects.empty()) {
        return {left_objects, right_objects};
      }
    }
  }

  // If we couldn't split along any axis, split arbitrarily in half
  size_t mid = objects.size() / 2;
  return {std::vector<Object3D *>(objects.begin(), objects.begin() + mid),
          std::vector<Object3D *>(objects.begin() + mid, objects.end())};
}

enum BVHMethod {
  TOP_DOWN,
  BOTTOM_UP
};
/**
 * @brief Recursively builds a BVH tree from a list of 3D objects.
 *
 * This method implements a **top-down recursive splitting** approach to
 * construct the BVH. It uses the **spatial median split** strategy, where
 * objects are divided along the longest axis of their combined bounding box.
 * The median split position is chosen as the midpoint of the bounding box on
 * the selected axis.
 *
 * @param objects A vector of pointers to Object3D instances.
 * @return A unique pointer to the root BVHNode.
 */
// Top-down BVH construction (existing method)
inline std::unique_ptr<BVHNode>
build_bvh_top_down(std::vector<Object3D *> &objects) {
  if (objects.empty()) {
    return nullptr;
  }

  if (objects.size() == 1) {
    return std::make_unique<BVHNode>(objects[0]->bounding_box, objects[0]);
  }

  auto [left_objects, right_objects] = split_objects(objects);

  // If splitting didn't reduce the problem size, create a leaf node with
  // multiple objects
  if (left_objects.size() == objects.size() ||
      right_objects.size() == objects.size()) {
    auto node = std::make_unique<BVHNode>(AABB());
    for (auto *obj : objects) {
      node->bounding_box = node->bounding_box.merge(obj->bounding_box);
    }
    // Store all objects in this node
    node->object = objects[0]; // Store the first object (you might want to
                               // modify this based on your needs)
    return node;
  }

  auto node = std::make_unique<BVHNode>(AABB());
  node->left = build_bvh_top_down(left_objects);
  node->right = build_bvh_top_down(right_objects);

  if (node->left) {
    node->bounding_box = node->bounding_box.merge(node->left->bounding_box);
  }
  if (node->right) {
    node->bounding_box = node->bounding_box.merge(node->right->bounding_box);
  }

  return node;
}

/**
 * @brief Builds a BVH (Bounding Volume Hierarchy) tree using a bottom-up
 * approach.
 *
 * This function constructs a BVH tree by iteratively merging leaf nodes over
 * the axis with the largest bounding box spread. It starts with individual
 * objects as leaf nodes and progressively merges them until a single root node
 * is formed.
 *
 * @param objects A reference to a vector of pointers to Object3D, representing
 *        the objects to be included in the BVH.
 * @return A unique pointer to the root BVHNode representing the constructed
 * BVH. Returns nullptr if the input vector is empty.
 */
inline std::unique_ptr<BVHNode>
build_bvh_bottom_up(std::vector<Object3D *> &objects) {
  if (objects.empty()) {
    return nullptr;
  }

  // Local struct for leaf nodes during construction
  struct BVHLeaf {
    AABB bbox;
    Object3D *object;
    std::unique_ptr<BVHNode> node;
    BVHLeaf() : object(nullptr), node(nullptr) { bbox = AABB(); }
    explicit BVHLeaf(std::unique_ptr<BVHNode> &node1)
        : bbox(node1->bounding_box), object(node1->object),node(std::move(node1))  {}
    BVHLeaf(Object3D *obj)
        : bbox(obj->bounding_box), object(obj),
          node(std::make_unique<BVHNode>(obj->bounding_box, obj)) {}
    bool empty() const {
      return ((object == nullptr) or (node == nullptr));
    }
  };

  // Helper function to find the axis with the largest spread
  auto find_best_axis = [](const std::vector<BVHLeaf> &leaves) {
    vec3 min_bounds(std::numeric_limits<double>::max());
    vec3 max_bounds(std::numeric_limits<double>::lowest());

    for (const auto &leaf : leaves) {

      min_bounds.x = std::min(min_bounds.x, leaf.bbox.min.x);
      min_bounds.y = std::min(min_bounds.y, leaf.bbox.min.y);
      min_bounds.z = std::min(min_bounds.z, leaf.bbox.min.z);
      max_bounds.x = std::max(max_bounds.x, leaf.bbox.max.x);
      max_bounds.y = std::max(max_bounds.y, leaf.bbox.max.y);
      max_bounds.z = std::max(max_bounds.z, leaf.bbox.max.z);
    }

    vec3 spread = {max_bounds.x - min_bounds.x, max_bounds.y - min_bounds.y,
                   max_bounds.z - min_bounds.z};

    return (spread.x > spread.y && spread.x > spread.z) ? 0
           : (spread.y > spread.z)                      ? 1
                                                        : 2;
  };

  std::vector<BVHLeaf> leaves;
  leaves.reserve(objects.size());
  for (auto *obj : objects) {
    leaves.emplace_back(obj);
  }
  while (leaves.size() > 1) {
    int best_axis = find_best_axis(leaves);
    std::sort(leaves.begin(), leaves.end(),
              [best_axis](const BVHLeaf &a, const BVHLeaf &b) {
                return a.bbox.min[best_axis] < b.bbox.min[best_axis];
              });

    std::vector<BVHLeaf> new_leaves;
    new_leaves.reserve(leaves.size() / 2 + 1);

    // new_leaves.resize(leaves.size() / 2 + 1);

    for (size_t i = 0; i < leaves.size(); i += 2) {
      if (i + 1 < leaves.size()) {
        auto new_node =
            std::make_unique<BVHNode>(leaves[i].bbox.merge(leaves[i + 1].bbox));
        new_node->left = std::move(leaves[i].node);
        new_node->right = std::move(leaves[i + 1].node);
        new_leaves.emplace_back(new_node);
      } else {
        new_leaves.push_back(std::move(leaves[i]));
      }
    }

    leaves = std::move(new_leaves);
  }

  return std::move(leaves[0].node);
}

// Wrapper function to allow choosing between top-down and bottom-up methods
std::unique_ptr<BVHNode> inline build_bvh(std::vector<Object3D *> &objects,
                                          bool use_bottom_up) {
  if (use_bottom_up) {
    return build_bvh_bottom_up(objects);
  } else {
    return build_bvh_top_down(objects);
  }
}


inline void intersect_bvh_iterative(
    const BVHNode* node1, const BVHNode* node2,
    std::vector<std::pair<const BVHNode*, const BVHNode*>>& intersections) {
    std::stack<std::pair<const BVHNode*, const BVHNode*>> stack;
    stack.emplace(node1, node2);

    while (!stack.empty()) {
        auto [current_node1, current_node2] = stack.top();
        stack.pop();
        
        if (!current_node1 || !current_node2 ||
            !current_node1->bounding_box.intersects(current_node2->bounding_box)) {
            continue;
        }

        if (current_node1->object && current_node2->object) {
            intersections.emplace_back(current_node1, current_node2);
            continue;
        }

        if (!current_node1->left && !current_node1->right) {
            if (current_node2->left) {
                stack.emplace(current_node1, current_node2->left.get());
            }
            if (current_node2->right) {
                stack.emplace(current_node1, current_node2->right.get());
            }
        } else if (!current_node2->left && !current_node2->right) {
            if (current_node1->left) {
                stack.emplace(current_node1->left.get(), current_node2);
            }
            if (current_node1->right) {
                stack.emplace(current_node1->right.get(), current_node2);
            }
        } else {
            if (current_node1->left && current_node2->left) {
                stack.emplace(current_node1->left.get(), current_node2->left.get());
            }
            if (current_node1->left && current_node2->right) {
                stack.emplace(current_node1->left.get(), current_node2->right.get());
            }
            if (current_node1->right && current_node2->left) {
                stack.emplace(current_node1->right.get(), current_node2->left.get());
            }
            if (current_node1->right && current_node2->right) {
                stack.emplace(current_node1->right.get(), current_node2->right.get());
            }
        }
    }
}
/**
 * @brief Finds intersections between two BVH trees.
 *
 * @param node1 Pointer to the root of the first BVH tree.
 * @param node2 Pointer to the root of the second BVH tree.
 * @param intersections A vector of intersecting BVH node pairs.
 */
void inline intersect_bvh(
    const BVHNode *node1, const BVHNode *node2,
    std::vector<std::pair<const BVHNode *, const BVHNode *>> &intersections) {
  if (!node1 || !node2 ) {
    return;
  }
  AABB aabb{};
  bool res=node1->bounding_box.intersection(node2->bounding_box,aabb) ;
  if (!res){return;}
  if (aabb.volume()<1e-6){
    return;
  }

  if (node1->object && node2->object) {
    intersections.emplace_back(node1, node2);
    return;
  }

  if (!node1->left && !node1->right) {
    intersect_bvh(node1, node2->left.get(), intersections);
    intersect_bvh(node1, node2->right.get(), intersections);
  } else if (!node2->left && !node2->right) {
    intersect_bvh(node1->left.get(), node2, intersections);
    intersect_bvh(node1->right.get(), node2, intersections);
  } else {
    intersect_bvh(node1->left.get(), node2->left.get(), intersections);
    intersect_bvh(node1->left.get(), node2->right.get(), intersections);
    intersect_bvh(node1->right.get(), node2->left.get(), intersections);
    intersect_bvh(node1->right.get(), node2->right.get(), intersections);
  }

}

/**
 * @brief Traverses a BVH tree to find objects intersecting with a target
 * bounding box.
 *
 * @param node Pointer to the current BVHNode.
 * @param target_bbox The target AABB to check for intersections.
 * @param results A vector to store the resulting objects that intersect with
 * the target bounding box.
 */
void traverse_bvh(const BVHNode *node, const AABB &target_bbox,
                  std::vector<const Object3D *> &results) {
  if (!node || !node->bounding_box.intersects(target_bbox)) {
    return;
  }

  if (node->object) {
    results.push_back(node->object);
    return;
  }

  traverse_bvh(node->left.get(), target_bbox, results);
  traverse_bvh(node->right.get(), target_bbox, results);
}

// Function to check if a point is contained in any object in the BVH
inline void contains_point(const BVHNode *bvh_root, const vec3 &pt,
                           std::vector<const Object3D *> &results) {

  AABB pt_box(pt, pt);
  traverse_bvh(bvh_root, pt_box, results);
}

// Helper function to create a formatted string for a vec3


// Recursive function to output BVH nodes information
inline std::string output_bvh_nodes(const BVHNode *node) {
  if (!node)
    return "";

  std::ostringstream oss;
  oss << "[" << format_vec3(node->bounding_box.min) << ","
      << format_vec3(node->bounding_box.max) << "]";

  if (node->left || node->right) {
    oss << ",[";
    if (node->left)
      oss << output_bvh_nodes(node->left.get());
    if (node->left && node->right)
      oss << ",";
    if (node->right)
      oss << output_bvh_nodes(node->right.get());
    oss << "]";
  }

  return oss.str();
}
inline std::string output_aabb(const AABB &obj) {

  return "[" + format_vec3(obj.min) + "," + format_vec3(obj.max) + "]";
}

inline std::string output_object3d(const Object3D *obj) {
  if (!obj)
    return "";
  return output_aabb(obj->bounding_box);
}

inline std::string output_objects3d(const std::vector<Object3D *> &obj) {
  std::string result = "";

  for (auto o : obj) {
    result += (output_object3d(o) + ",");
  }

  return "[" + result + "]";
}



}
#endif // BVH_HPP