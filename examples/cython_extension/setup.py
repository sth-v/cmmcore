import sys
from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy as np
from pathlib import Path
root=Path(__file__).resolve().parent.parent.parent.resolve()
print(root)
include_dirs=[str(root),*sys.path]
compile_args = ["-O3", "-flto", "-std=c++20"]
sys.path.append(str(root))

if sys.platform=="darwin" :
    sys.path.extend(['/opt/homebrew/opt/llvm/lib','/opt/homebrew/opt/llvm/bin','/opt/homebrew/opt/llvm/include'])
    compile_args+=["-mcpu=native","-march=armv8-a+simd"]


else:
    compile_args+=[
    '-msse', '-msse2', '-msse3', '-mssse3',
    '-msse4.1', '-msse4.2', '-mavx', '-mavx2'
  ]


compiler_directives = dict(
    boundscheck=False,
    wraparound=False,
    cdivision=True,
    nonecheck=False,
    overflowcheck=False,
    initializedcheck=False,
    embedsignature=False,
    language_level="3str",
)

setup(ext_modules = cythonize(Extension(
           "ccx",                       # the extension name
           sources=["ccx.pyx"],
           include_dirs=include_dirs,         # the Cython source and
            extra_compile_args=compile_args,  # additional C++ source files
           language="c++",                    # generate and compile C++ code
      ),
      compiler_directives=compiler_directives
    )
)