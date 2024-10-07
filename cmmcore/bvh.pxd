# cython: language_level=3
# distutils: language = c++
cimport cython

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp cimport bool

cdef extern from "bvh.hpp":


    cdef cppclass AABB:
        AABB()
        AABB(const vec3 &min, const vec3 &max)
        bool infinity() const
        bool intersects(const AABB &other) const
        bool intersection(const AABB &other, AABB &result) const
        AABB merge(const AABB &other) const
        double volume() const
        void expand(const vec3& point)

    cdef cppclass Object3D:
        Object3D()
        Object3D(const AABB &bb)
        Object3D(const AABB &bb, int id)

    cdef cppclass BVHNode:
        BVHNode()
        BVHNode(const BVHNode& other)
        BVHNode(const AABB &bb, Object3D *obj = nullptr)
        void setLeft(unique_ptr[BVHNode] node)
        void setRight(unique_ptr[BVHNode] node)
        BVHNode& operator=(const BVHNode& other)

    vec3 compute_extents(const AABB &bbox)
    ctypedef  Object3D* Object3D_Ptr
    ctypedef  const BVHNode* const_BVHNode_Ptr
    pair[vector[Object3D_Ptr], vector[Object3D_Ptr]] split_objects(vector[Object3D *] &objects)

    unique_ptr[BVHNode] build_bvh_top_down(vector[Object3D *] &objects)
    unique_ptr[BVHNode] build_bvh_bottom_up(vector[Object3D *] &objects)

    unique_ptr[BVHNode] build_bvh(vector[Object3D *] &objects, bool use_bottom_up)

    void intersect_bvh_iterative(const BVHNode* node1, const BVHNode* node2, vector[pair[const_BVHNode_Ptr, const_BVHNode_Ptr]]& intersections)
    void intersect_bvh(const BVHNode *node1, const BVHNode *node2, vector[pair[const_BVHNode_Ptr, const_BVHNode_Ptr]] &intersections)

    void traverse_bvh(const BVHNode *node, const AABB &target_bbox, vector[const Object3D *] &results)
    void contains_point(const BVHNode *bvh_root, const vec3 &pt, vector[const Object3D *] &results)

    string output_bvh_nodes(const BVHNode *node)
    string output_aabb(const AABB &obj)
    string output_object3d(const Object3D *obj)
    string output_objects3d(const vector[Object3D *] &obj)

