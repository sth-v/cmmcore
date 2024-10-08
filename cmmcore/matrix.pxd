# cython: language_level=3
# distutils: language = c++
cimport cython
from libcpp.vector cimport vector
cdef extern from "matrix.h" namespace "cmmcore" nogil:

    cdef cppclass Matrix:
        Matrix(size_t r, size_t c)

        double& operator()(size_t i, size_t j)
        const double& operator()(size_t i, size_t j) const

        size_t getRows() const
        size_t getCols() const

        Matrix operator*(const Matrix& other) const

        Matrix transpose() const

        Matrix inverse() const

    ctypedef vector[Matrix] Tensor3D

