# Cython Extension example
This example illustrates how to build a simple python cython extension using cmmcore. 
## Project structure
A simple python extension structure might look like this:
```
.
├─── ccx.pyx # Cython source file. 
├─── main.py # Python program.
├─── setup.py # Build script for the extension
└─── requirenments.txt # python requirenments

```
In this example, there is additionally a `build.sh` file containing all the commands for building. But since it re-creates a virtual environment, you probably don't want to take it into the actual project.

## ccx.pyx
The extension file, contains just one function `pyccx`, which implements a wrapper over the NURBS two-curve intersection function `cmmcore::ccx`. 

```python

def pyccx( pts1, pts2, degree1=3, degree2=3,tol=1e-5):
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


```
In this example, we will not create a wrapper of the `cmmcore::NURBSCurve` class, so our function will only take control points and degree for each of the curves to be intersected. 
If you are familiar with cython you can easily extend these definitions.

## Build
1. Navigate to the given folder 
    ```bash
    cd examples/cython_extension ## If you are in the root of the project
    ```
2. Create and activate the virtual environment
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    ```
3. Install the necessary dependencies
    ```bash
    python3 -m pip install -r requirenments.txt
    ```
4. Build the extension
    ```bash
    python3 setup.py build_ext --inplace
    ```
## Run
Make sure everything works:
```bash
python3 main.py
```
Output:
```
[(0.3237280994653702, 0.22155483439564705), (1.2614987790584564, 0.786673478782177), (2.3095895648002625, 3.0814577788114548), (4.2597944140434265, 6.622736647725105), (2.8862611055374146, 7.53760227560997), (8.019260108470917, 10.729096576571465), (10.308691799640656, 12.926283463835716), (0.06767705082893372, 18.62513067573309), (6.138576149940491, 5.043478101491928)]
ccx calculated at: 0.5918329999999999 ms. 
```
**0.59** milliseconds - pretty fast for my machine.

