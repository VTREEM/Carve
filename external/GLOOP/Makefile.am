noinst_LTLIBRARIES=libgloop-math.la libgloop-image.la libgloop-model.la
noinst_HEADERS=

AM_CPPFLAGS=@CPPFLAGS@ -Iinclude -I@top_srcdir@/external/GLEW/include

libgloop_math_la_SOURCES=src/matrix.cpp src/quaternion.cpp

libgloop_image_la_SOURCES=src/radiance.cpp

libgloop_model_la_SOURCES=src/model/ply_format.cpp src/model/stream.cpp src/model/obj_format.cpp src/model/vtk_format.cpp

if with_GUI
  noinst_LTLIBRARIES+=libgloop.la

  libgloop_la_CXXFLAGS= -Iinclude -I@top_srcdir@/external/GLEW/include @GL_CFLAGS@
  libgloop_la_SOURCES=src/shader.cpp src/surface.cpp src/texparam.cpp src/vbo.cpp src/fbo.cpp
endif
