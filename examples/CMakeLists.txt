include_directories("${carve_SOURCE_DIR}/include")
include_directories("${carve_SOURCE_DIR}/external/GLOOP/include")
include_directories("${carve_SOURCE_DIR}/external/GLEW/include")
include_directories("${carve_SOURCE_DIR}/external/GLUI/include")
include_directories("${carve_SOURCE_DIR}/common")
include_directories(${OPENGL_INCLUDE_DIR})

link_directories("${gloop_BINARY_DIR}")
link_directories("${carve_BINARY_DIR}/common")

if(CARVE_WITH_GUI)
  add_executable       (custom_collector   custom_collector.cpp)
  target_link_libraries(custom_collector   carve carve_ui carve_misc glui carve_fileformats gloop_model ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

  add_executable       (texture_example    texture_example.cpp)
  target_link_libraries(texture_example    carve carve_ui carve_misc glui carve_fileformats gloop_model ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

  add_executable       (detail_clip        detail_clip.cpp)
  target_link_libraries(detail_clip        carve carve_misc carve_fileformats gloop_model)
endif(CARVE_WITH_GUI)
