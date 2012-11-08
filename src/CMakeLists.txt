include_directories("${carve_SOURCE_DIR}/include")
include_directories("${carve_SOURCE_DIR}/external/GLOOP/include")
include_directories("${carve_SOURCE_DIR}/external/GLEW/include")
include_directories("${carve_SOURCE_DIR}/external/GLUI/include")
include_directories("${carve_SOURCE_DIR}/common")
include_directories(${OPENGL_INCLUDE_DIR})

link_directories("${gloop_BINARY_DIR}")
link_directories("${carve_BINARY_DIR}/common")

if(CARVE_INTERSECT_GLU_TRIANGULATOR AND CARVE_WITH_GUI)
  add_executable       (intersect   glu_triangulator.cpp intersect.cpp)
  target_link_libraries(intersect   carve_fileformats carve_misc carve gloop_model ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
  add_executable       (slice       glu_triangulator.cpp slice.cpp)
  target_link_libraries(slice       carve_fileformats carve_misc carve gloop_model ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
else(CARVE_INTERSECT_GLU_TRIANGULATOR AND CARVE_WITH_GUI)
  add_definitions(-DDISABLE_GLU_TRIANGULATOR)
  add_executable       (intersect   intersect.cpp)
  target_link_libraries(intersect   carve_fileformats carve_misc carve gloop_model)
  add_executable       (slice       slice.cpp)
  target_link_libraries(slice       carve_fileformats carve_misc carve gloop_model)
endif(CARVE_INTERSECT_GLU_TRIANGULATOR AND CARVE_WITH_GUI)

if(CARVE_WITH_GUI)
  add_executable       (view   view.cpp)
  target_link_libraries(view   carve_fileformats carve_misc carve_ui carve gloop_model glui ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
  install(TARGETS view RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
endif(CARVE_WITH_GUI)

add_executable       (face_merge face_merge.cpp)
target_link_libraries(face_merge carve_fileformats carve gloop_model)

add_executable       (mesh_simplify mesh_simplify.cpp)
target_link_libraries(mesh_simplify carve_fileformats carve gloop_model)

add_executable       (triangulate triangulate.cpp)
target_link_libraries(triangulate carve_fileformats carve gloop_model)

add_executable       (cutgraph    cut.cpp)
target_link_libraries(cutgraph    carve_fileformats carve gloop_model)

add_executable       (extrude     extrude.cpp)
target_link_libraries(extrude     carve_fileformats carve gloop_model)

add_executable       (convert     convert.cpp)
target_link_libraries(convert     carve_fileformats carve gloop_model)

add_executable       (close_manifold     close_manifold.cpp)
target_link_libraries(close_manifold     carve_fileformats carve gloop_model)

add_executable       (selfintersect     selfintersect.cpp)
target_link_libraries(selfintersect     carve_fileformats carve gloop_model)

foreach(tgt slice intersect triangulate convert)
  install(TARGETS ${tgt}
          RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
endforeach(tgt)
