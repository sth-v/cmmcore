cdef extern from "quickhull.hpp" namespace "cmmcore" nogil:

    cdef cppclass quick_hull[point_iterator, value_type]:
        quick_hull(size_type, value_type const&&) except * # Deleted constructor
        quick_hull(size_type const _dimension, value_type const & _eps)

        value_type cos_of_dihedral_angle(facet const & _first, facet const & _second) const

        void make_facet(facet & _facet, point_array const & _vertices, size_type const _against, point_iterator const _apex, size_type const _neighbour)
        void reuse_facet(facet & _facet, point_array const & _vertices, size_type const _against, point_iterator const _apex, size_type const _neighbour)

        void copy_point(point_iterator const _from, vrow _to) const
        void subtract(vrow _minuend, crow _subtrahend) const

        void gshift(vrow _augend, value_type const & _addend) const
        void divide(vrow _dividend, value_type const & _divisor) const
        void subtract_and_assign(vrow _assignee, vrow _minuend, crow _subtrahend) const
        void multiply_and_add(vrow _assignee, crow _multiplicand, value_type const & _factor) const
        void scale_and_shift(vrow _multiplicand, crow _direction, value_type const & _factor) const

        void matrix_transpose_copy(point_array const & _vertices)
        void matrix_restore(size_type const _identity)
        void matrix_restore()
        void matrix_sqr(size_type const _size)
        value_type det(matrix const & _matrix, size_type const _dimension)
        value_type det()

        void set_hyperplane_equation(facet & _facet)
        bool orthonormalize(point_list const & _affine_space, size_type const _rank, crow const _origin)
        void forward_transformation(size_type const _rank)
        bool steal_best(point_list & _basis)

    cdef cppclass facet:
        template typename iterator value_type distance(iterator const _point) const

    # Add additional methods following the established pattern

