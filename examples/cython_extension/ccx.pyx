cimport cython
from libcpp.vector cimport vector
from libcpp.pair cimport pair

from cmmcore.vec cimport vec3,vec4
from cmmcore.nurbs cimport NURBSCurve
from cmmcore.ccx cimport ccx


def pyccx(double[:,:] pts1,double[:,:] pts2,int degree1=3,int degree2=3, double tol=1e-5):
    cdef vector[vec4] ptsc1=vector[vec4](pts1.shape[0])
    cdef vector[vec4] ptsc2=vector[vec4](pts2.shape[0])
    cdef int i
    for i in range(pts1.shape[0]):

        ptsc1[i].x=pts1[i,0]
        ptsc1[i].y=pts1[i,1]
        ptsc1[i].z=pts1[i,2]
        if pts1.shape[1]==4:
            ptsc1[i].w=pts1[i,3]
        else:
            ptsc1[i].w=1
    for i in range(pts2.shape[0]):

        ptsc2[i].x=pts2[i, 0]
        ptsc2[i].y=pts2[i, 1]
        ptsc2[i].z=pts2[i, 2]
        if pts2.shape[1] == 4:
            ptsc2[i].w=pts2[i, 3]
        else:
            ptsc2[i].w=1
    cdef NURBSCurve nc1=NURBSCurve(ptsc1, degree1)
    cdef NURBSCurve nc2=NURBSCurve(ptsc2, degree2)
    cdef vector[pair[double,double]] intersections
    ccx(nc1,nc2,tol,intersections)
    cdef list py_ints=[]
    cdef tuple temp
    for i in range(intersections.size()):
        temp=(intersections[i].first, intersections[i].second)
        py_ints.append(temp)
    return py_ints

