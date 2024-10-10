# cmmcore
The concise, cutting-edge CAD engine.

## Overview

## Build

### Native
cmmcore is a header-only library and therefore does not need to be compiled into a static or dynamic library. You can still do this with meson, but it is not required. 
The files in the test directories demonstrate how to use cmmcore in your own native executables. 
To build them run :
```bash
make 
```
### Cython/Python
**cmmcore** reveals declarations in `.pxd` files. Which makes it also a cython library that you can use in your cython projects. An example of such use can be found in [`./examples/cython_extension`](./examples/cython_extension/README.md) directory.

### WASM
**cmmcore** can be compiled into wasm using emscripten.
> WASM buindings are not complete at the moment, a lot of work is needed to make them usable. The bindings currently contains a few working examples and will be extended in the future.

To build and run the wasm module, follow these steps:
1. Make sure emscripten is installed and available in your environment
2. Go to the ./wasm directory
    ```bash
    cd wasm
    ```
3. Compile cmmcore.wasm and cmmcore.js:
    ```bash
    make
    ```
4. Run file server
    ```bash
    python3 -m http.server
    ```
   You'll see something like this:
   ```
   Serving HTTP on :: port 8000 (http://[::]:8000/) ...
   ```

5. In your browser, go to http://localhost:8000/index.html and open the browser console. In the browser console, you should see:
    ```
    Evaluated surface point: 0.386633295312500070.53927819546875010.3723660509375
    ```
