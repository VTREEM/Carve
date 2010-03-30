include_directories("${carve_SOURCE_DIR}/include")
include_directories("${carve_SOURCE_DIR}/external/GLOOP/include")
include_directories("${carve_SOURCE_DIR}/external/GLEW/include")
include_directories("${carve_SOURCE_DIR}/external/GLUI/include")
include_directories("${carve_SOURCE_DIR}/common")
include_directories(${OPENGL_INCLUDE_DIR})

link_directories("${gloop_BINARY_DIR}")
link_directories("${carve_BINARY_DIR}/common")

if(ENABLE_GLU_TRI)
  add_executable       (intersect   glu_triangulator.cpp intersect.cpp)
  target_link_libraries(intersect   carve_fileformats carve_misc carve gloop_model ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
else(ENABLE_GLU_TRI)
  add_definitions(-DDISABLE_GLU_TRIANGULATOR)
  add_executable       (intersect   intersect.cpp)
  target_link_libraries(intersect   carve_fileformats carve_misc carve gloop_model)
endif(ENABLE_GLU_TRI)

add_executable       (view   view.cpp)
target_link_libraries(view   carve_fileformats carve_misc carve_ui carve gloop_model glui ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

add_executable       (triangulate triangulate.cpp)
target_link_libraries(triangulate carve_fileformats carve gloop_model)

add_executable       (cutgraph    cut.cpp)
target_link_libraries(cutgraph    carve_fileformats carve gloop_model)

add_executable       (convert     convert.cpp)
target_link_libraries(convert     carve_fileformats carve gloop_model)

foreach(tgt intersect view triangulate convert)
  install(TARGETS ${tgt}
          RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
endforeach(tgt)
