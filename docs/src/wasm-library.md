#  Building WASM library 
**cmmcore** can be compiled into wasm using emscripten.
> WASM buindings are not complete at the moment, a lot of work is needed to make them usable. The bindings currently contains a few working examples and will be extended in the future.

```js,editable
   cmmcore().then(cmm => {
        const vec3 = cmm.vec3;
        const vec4 = cmm.vec4;
        const NURBSSurface = cmm.NURBSSurface;
        const controlPointsArray = [
            [[0.62, 0.21, 0.75, 1.], [0.13, 0.85, 0.19, 1.], [0.89, 0.64, 0.07, 1.], [0.91, 0.41, 0.52, 1.]],
            [[0.50, 0.88, 0.78, 1.], [0.25, 0.33, 0.46, 1.], [0.71, 0.56, 0.05, 1.], [0.13, 0.81, 0.33, 1.]],
            [[0.75, 0.64, 0.14, 1.], [0.35, 0.60, 0.11, 1.], [0.17, 0.12, 0.91, 1.], [0.06, 0.53, 0.10, 1.]],
            [[0.05, 0.59, 0.02, 1.], [0.17, 0.86, 0.62, 1.], [0.42, 0.73, 0.41, 1.], [0.33, 0.57, 0.53, 1.]]
        ];


        // Create some control points (2x2 grid of vec4)
        const controlPoints = cmm.create2DVector(4, 4, new vec4(0.0, 0.0, 0.0, 1.0))

        // Create control points as vector<vector<vec4>>

        // Row 1 of control points

        cmm.convertJSArrayToNURBSSurfaceControlPoints(controlPointsArray, controlPoints)

        // Define the degree of the surface
        const degree = [3, 3];

        // Optionally, define the knot vectors (empty in this case)
        const knotsU = new cmm.vector_double();
        const knotsV = new cmm.vector_double();

        // Create a NURBSSurface instance
        const surface = new NURBSSurface(controlPoints, degree, knotsU, knotsV);

        // Evaluate the surface at u=0.5, v=0.5
        const result = new vec3();
        surface.evaluate(0.5, 0.5, result);

        console.log("Evaluated surface point:", result.x, result.y, result.z);
    });

```

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
    Evaluated surface point: 0.38093750,532656250,36734375
    ```
