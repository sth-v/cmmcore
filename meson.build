project('cmmcore', 'cpp', default_options : ['cpp_std=c++17', 'optimization=3', 'buildtype=release'])
# Create an environment object
env = environment()
# Extend the PATH variable

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


my_inc = include_directories('.')


libcmmcore = shared_library('cmmcore',
                            'cmmcore/nurbs_utils.h',
                            'cmmcore/nurbs.h',
   'cmmcore/ccx.h',
                            'cmmcore/binom.h',
                            'cmmcore/vec.h',
                            'cmmcore/convexhull.h',
                            'cmmcore/polygon.h',
                            'cmmcore/polyline.h',
          'cmmcore/monomial.h',
          'cmmcore/gauss_map.h',
                            version : '0.1.4',
                             install : false,
                            include_directories : my_inc)

libcmmcore_dep = declare_dependency(link_with : libcmmcore, include_directories : my_inc)


subdir('test')






