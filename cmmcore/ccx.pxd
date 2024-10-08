# cython: language_level=3
# distutils: language = c++
cimport cython
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from cmmcore.nurbs cimport NURBSCurve

cdef extern from "ccx.h" namespace "cmmcore" nogil:
    cdef void ccx(NURBSCurve& c1, NURBSCurve& c2, double tol, vector[pair[double, double]]& result)
