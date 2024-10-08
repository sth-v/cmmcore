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
5. Make sure everything works
    ```bash
    python3 main.py
    ```
    Output:
    ```
    [(0.3237280994653702, 0.22155483439564705), (1.2614987790584564, 0.786673478782177), (2.3095895648002625, 3.0814577788114548), (4.2597944140434265, 6.622736647725105), (2.8862611055374146, 7.53760227560997), (8.019260108470917, 10.729096576571465), (10.308691799640656, 12.926283463835716), (0.06767705082893372, 18.62513067573309), (6.138576149940491, 5.043478101491928)]
    ccx calculated at: 0.5918329999999999 ms. 
    ```
    0.59 milliseconds - pretty fast for my machine.

