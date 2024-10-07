# cython: language_level=3
# distutils: language = c++
cimport cython
from libcpp.vector cimport vector
from libcpp.array cimport array
from libcpp.pair cimport pair
from libcpp cimport bool
from cmmcore.vec cimport vec3,vec4
from cmmcore.aabb cimport AABB



cdef extern from "nurbs.h" nogil:
    cdef cppclass NURBSCurve:
        NURBSCurve()
        NURBSCurve(const vector[vec4] &control_points)
        NURBSCurve(const vector[vec4] &control_points, bool periodic)
        NURBSCurve(const vector[vec3] &control_points, int degree)
        NURBSCurve(const vector[vec4] &control_points, int degree)
        NURBSCurve(const vector[vec3] &control_points, int degree, bool periodic)
        NURBSCurve(const vector[vec4] &control_points, int degree, bool periodic)
        NURBSCurve(const vector[vec4] &control_points, int degree, const vector[double] &knots, bool periodic)
        NURBSCurve(const NURBSCurve &other)

        NURBSCurve &operator=(const NURBSCurve &other)

        bool is_periodic() const
        void _update_interval()
        void knots_update_hook()
        void normalize_knots()
        void generate_knots()
        void make_periodic()
        void generate_knots_periodic()
        void insert_knot(double t, int count)
        pair[NURBSCurve, NURBSCurve] split(double param, bool normalize_knots=*) const
        void evaluate(double t, vec3 &result) const
        const vector[vec4] get_control_points()
        void set_control_points(vector[vec4] &cpts)
        int get_degree()
        void set_degree(int val)
        const vector[double] get_knots()
        void set_knots(vector[double] &knots)
        array[double, 2] interval()
        AABB &aabb()
        void rebuildAABB()

    cdef cppclass NURBSSurface:
        NURBSSurface(const vector[vector[vec4]] &control_points, const array[int, 2] &degree, const vector[double] &knots_u = *, const vector[double] &knots_v = *)

        void generate_knots_u()
        void generate_knots_v()
        void evaluate(double u, double v, vec4 &result) const noexcept
        void insert_knot_u(double t, int count)
        void insert_knot_v(double t, int count)
        pair[NURBSSurface, NURBSSurface] split_surface_u(double param, double tol=*) const
        pair[NURBSSurface, NURBSSurface] split_surface_v(double param, double tol=*) const

        void _update_interval() noexcept