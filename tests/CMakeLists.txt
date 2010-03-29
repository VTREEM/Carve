include_directories("${carve_SOURCE_DIR}/include")
include_directories("${carve_SOURCE_DIR}/external/GLOOP/include")
include_directories("${carve_SOURCE_DIR}/external/GLEW/include")
include_directories("${carve_SOURCE_DIR}/external/GLUI/include")
include_directories("${carve_SOURCE_DIR}/common")
include_directories(${OPENGL_INCLUDE_DIR})

link_directories("${gloop_BINARY_DIR}")
link_directories("${carve_BINARY_DIR}/common")

add_executable       (test_geom          test_geom.cpp)
target_link_libraries(test_geom          carve)

add_executable       (test_eigen         test_eigen.cpp)
target_link_libraries(test_eigen         carve)

add_executable       (test_spacetree     test_spacetree.cpp)
target_link_libraries(test_spacetree     carve)

add_executable       (test_aabb          test_aabb.cpp)
target_link_libraries(test_aabb          carve)

add_executable       (test_rescale       test_rescale.cpp)
target_link_libraries(test_rescale       carve)

add_executable       (tetrahedron        tetrahedron.cpp)
target_link_libraries(tetrahedron        carve)

add_executable       (test_intersect     test_intersect.cpp)
target_link_libraries(test_intersect     carve carve_fileformats carve_ui carve_misc glui gloop_model ${OPENGL_LIBRARY})

add_executable       (test_interpolate   test_interpolate.cpp)
target_link_libraries(test_interpolate   carve carve_ui glui ${OPENGL_LIBRARY} ${GLUT_LIBRARY})

add_executable       (test_csg_interpolate   test_csg_interpolate.cpp)
target_link_libraries(test_csg_interpolate   carve carve_ui carve_misc glui ${OPENGL_LIBRARY} ${GLUT_LIBRARY})

add_executable       (test_slice             test_slice.cpp)
target_link_libraries(test_slice             carve carve_ui carve_misc glui ${OPENGL_LIBRARY})

add_executable       (test_slice_classify    test_slice_classify.cpp)
target_link_libraries(test_slice_classify    carve carve_ui carve_misc glui ${OPENGL_LIBRARY})

add_executable       (test_triangulate       test_triangulate.cpp)
target_link_libraries(test_triangulate       carve carve_ui carve_misc glui ${OPENGL_LIBRARY})

add_executable       (test_hole_incorporate  test_hole_incorporate.cpp)
target_link_libraries(test_hole_incorporate  carve carve_ui carve_misc glui ${OPENGL_LIBRARY})
