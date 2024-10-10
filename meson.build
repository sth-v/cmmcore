project('cmmcore', 'cpp', default_options : ['cpp_std=c++17', 'optimization=3', 'buildtype=release'])
# Create an environment object
env = environment()
my_inc = include_directories('.')

cmmcore_sources=[   'cmmcore/vec.h',
                                'cmmcore/binom.h',
                                'cmmcore/nurbs_utils.h',
                                'cmmcore/nurbs.h',
                                'cmmcore/ccx.h',
                                'cmmcore/convexhull.h',
                                'cmmcore/polygon.h',
                                'cmmcore/polyline.h',
                                'cmmcore/monomial.h',
                                'cmmcore/gauss_map.h']
cmmcore_wasm_binding_sources=['wasm/cmmcore.cpp']

# Extend the PATH variable
if host_machine.system() != 'emscripten'
    base_compiler_args=['-flto','-lm']
    base_compiler_args=[]
    # Detect the host operating system
    if host_machine.system() == 'darwin'
        simd_args = ['-mcpu=apple-m1','-march=armv8-a+simd']



    else
        simd_args = [
            '-msse', '-msse2', '-msse3', '-mssse3',
            '-msse4.1', '-msse4.2', '-mavx', '-mavx2'
        ]
    endif

    # Add the SIMD arguments to the project
    add_project_arguments(base_compiler_args+simd_args, language : 'cpp')





    libcmmcore = shared_library('cmmcore',
                               cmmcore_sources,

                                version : '0.1.4',
                                 install : false,
                                include_directories : my_inc)

    libcmmcore_dep = declare_dependency(link_with : libcmmcore, include_directories : my_inc)


    subdir('test')
endif


# Build the Wasm (Emscripten) executable
if host_machine.system() == 'emscripten'

    cmmcore_wasm = executable(
        'cmmcore_wasm',
        cmmcore_sources+cmmcore_wasm_binding_sources,

        override_options: ['cpp_args=-O3','cpp_args=-s', 'cpp_args=MODULARIZE=1', 'cpp_args=-s', 'cpp_args=EXPORT_NAME="cmmcore"', 'cpp_args=--bind'],
        include_directories : my_inc,
        install : false

    )


endif


