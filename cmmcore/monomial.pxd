# cython: language_level=3
# distutils: language = c++
cimport cython
from libcpp.pair cimport pair
from cmmcore.matrix cimport Matrix, Tensor3D

cdef extern from "monomial.h" namespace "cmmcore" nogil:
    Matrix bpmat(int n)
    Tensor3D bezier_to_monomial(const Tensor3D& control_points)
    Tensor3D monomial_to_bezier(const Tensor3D& monomial_coeffs)
    Tensor3D cross_product_monomial(const Tensor3D& a_coeffs, const Tensor3D& b_coeffs)
    pair[Tensor3D, Tensor3D] monomial_partial_derivatives(const Tensor3D& coeffs)
    Tensor3D normal_vector_monomial(const Tensor3D& coeffs)

