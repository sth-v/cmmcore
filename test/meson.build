
#cmm_mesh_test=executable('cmm_mesh_test',
#                         'cmm_mesh_test.cpp',      override_options : ['cpp_std=c++17'] ,link_with : libcmmcore,  include_directories : my_inc, install:false)
#cmm_numeric_test=executable('cmm_numeric_test',
#                         'cmm_numeric_test.cpp',      override_options : ['cpp_std=c++17'] ,link_with : libcmmcore,  include_directories : my_inc,install:false)
#cmm_cpt_test=executable('cmm_cpt_test',
#                         'cmm_cpt_test.cpp',      override_options : ['cpp_std=c++17'] ,link_with : libcmmcore,  include_directories : my_inc, install:false)
#
#
#deboor_tests=executable('deboor_tests',
#                         'deboor_tests.c',      override_options : ['cpp_std=c17'] ,link_with : libcmmcore,  include_directories : my_inc, install:false)
#
#cmm_algo_test=executable('cmm_algo_test',
#                         'cmm_algo_test.cpp',      override_options : ['cpp_std=c++17'] ,link_with : libcmmcore,  include_directories : my_inc, install:false)
#vec_perf=executable('vec_perf',
#                         'vecPerf.cpp',      override_options : ['cpp_std=c++17'] ,link_with : libcmmcore,  include_directories : my_inc, install:false)
#test('cmm_algo_test', cmm_algo_test)
#test('vec_perf', vec_perf)
#test('cmm_cpt_test', cmm_cpt_test)
#test('cmm_numeric_test', cmm_numeric_test)
INSTALL=false
INC_DIRS=my_inc

# find_span_test
find_span_test=executable('find_span_test',
                          'find_span_test.cpp',
                          include_directories : INC_DIRS,
                          install:INSTALL
)
test('find_span_test',
     find_span_test)

# evaluate_nurbs_test
evaluate_nurbs_test = executable('evaluate_nurbs_test',
                                 'evaluate_nurbs_test.cpp',
                                 include_directories : my_inc,
                                 install : false
)
test('evaluate_nurbs_test',
     evaluate_nurbs_test)

# curve_derivative_test
curve_derivative_test = executable('curve_derivative_test',
                                   'curve_derivative_test.cpp',
                                   include_directories : my_inc,
                                   install : false
)
test('curve_derivative_test',
     curve_derivative_test)

# surface_derivative_test
surface_derivative_test = executable('surface_derivative_test',
                                     'surface_derivative_test.cpp',
                                     include_directories : my_inc,
                                     install : false
)
test('surface_derivative_test',
     surface_derivative_test)

# knots_test
knots_test = executable('knots_test',
                        'knots_test.cpp',
                        include_directories : my_inc,
                        install : false
)
test('knots_test',
     knots_test)

# cxx_test
cxx_test = executable('cxx_test',
               'cxx_test.cpp',
               include_directories : my_inc,
               install : false,link_with: libcmmcore)

test('cxx_test',
     cxx_test)

# monomial_test
monomial_test = executable('monomial_test',
               'monomial_test.cpp',
               include_directories : my_inc,
               install : false
)
test('monomial_test',
     monomial_test)
# chull2d_test
chull2d_test = executable('chull2d_test',
                          'chull2d_test.cpp',
                          include_directories : my_inc,
                          install : false
)
test('chull2d_test',
     chull2d_test)

# polygon_test
polygon_test = executable('polygon_test',
                          'polygon_test.cpp',
                          include_directories : my_inc,
                          install : false
)
test('polygon_test',
     polygon_test)
# gauss_map_test
gauss_map_test = executable('gauss_map_test',
                            'gauss_map_test.cpp',
                            include_directories : my_inc,
                            install : false
)
test('gauss_map_test',
     gauss_map_test)

# ssx_test
ssx_test = executable('ssx_test',
                      'ssx_test.cpp',
                      include_directories : my_inc,
                      install : false
)
test('ssx_test',
     ssx_test)


# spatial_separability_test
spatial_separability_test = executable('spatial_separability_test',
                                       'spatial_separability_test.cpp',
                                       include_directories : my_inc,
                                       install : false
)
test('spatial_separability_test',
     spatial_separability_test)
# spherical_separability_test
spherical_separability_test = executable('spherical_separability_test',
                                         'spherical_separability_test.cpp',
                                         include_directories : my_inc,
                                         install : false
)
test('spherical_separability_test',
     spherical_separability_test)
