bin_PROGRAMS = intersect triangulate convert cut
noinst_PROGRAMS = test_geom test_eigen test_spacetree test_aabb test_aabb_tri test_rescale tetrahedron
noinst_LTLIBRARIES = libcarve_fileformats.la libcarve_misc.la

noinst_HEADERS= geom_draw.hpp geometry.hpp glu_triangulator.hpp		\
	opts.hpp read_ply.hpp rgb.hpp scene.hpp stringfuncs.hpp		\
	write_ply.hpp mersenne_twister.h coords.h brick_texture.h	\
	carve_texture.h leaf_texture.h

AM_CPPFLAGS=@CPPFLAGS@ -I$(top_srcdir)/include @GL_CFLAGS@ @GLUT_CFLAGS@
AM_CPPFLAGS += -I$(top_srcdir)/external/GLOOP/include
AM_LDFLAGS=@LDFLAGS@

libcarve_fileformats_la_SOURCES = read_ply.cpp write_ply.cpp
libcarve_fileformats_la_LIBADD = $(top_srcdir)/external/GLOOP/libgloop-model.la
libcarve_fileformats_la_CPPFLAGS = -I$(top_srcdir)/external/GLOOP/include

libcarve_misc_la_SOURCES=geometry.cpp

if with_GUI
  bin_PROGRAMS += view
  AM_CPPFLAGS += -I$(top_srcdir)/external/GLEW/include
  AM_CPPFLAGS += -I$(top_srcdir)/external/GLUI/include
  noinst_LTLIBRARIES += libcarve_ui.la
  noinst_PROGRAMS += custom_collector test_intersect test_interpolate	\
	  test_csg_interpolate test_slice test_slice_classify		\
	  test_triangulate texture_example problem
  libcarve_ui_la_SOURCES = geom_draw.cpp scene.cpp
  libcarve_ui_la_LIBADD = $(top_srcdir)/external/GLUI/libglui.la
endif


test_eigen_SOURCES=test_eigen.cpp
test_eigen_LDADD=../lib/libintersect.la

test_spacetree_SOURCES=test_spacetree.cpp
test_spacetree_LDADD=../lib/libintersect.la

test_aabb_SOURCES=test_aabb.cpp
test_aabb_LDADD=../lib/libintersect.la

test_aabb_tri_SOURCES=test_aabb_tri.cpp
test_aabb_tri_LDADD=../lib/libintersect.la

test_triangulate_SOURCES=test_triangulate.cpp
test_triangulate_LDADD=libcarve_ui.la libcarve_misc.la ../lib/libintersect.la @GL_LIBS@ @GLUT_LIBS@

custom_collector_SOURCES=custom_collector.cpp
custom_collector_LDADD=libcarve_fileformats.la libcarve_ui.la libcarve_misc.la ../lib/libintersect.la @GL_LIBS@ @GLUT_LIBS@

test_intersect_SOURCES=test_intersect.cpp
test_intersect_LDADD=libcarve_fileformats.la libcarve_ui.la libcarve_misc.la ../lib/libintersect.la @GL_LIBS@ @GLUT_LIBS@

triangulate_SOURCES=triangulate.cpp
triangulate_LDADD=libcarve_fileformats.la ../lib/libintersect.la

cut_SOURCES=cut.cpp
cut_LDADD=libcarve_misc.la ../lib/libintersect.la libcarve_fileformats.la

intersect_SOURCES=intersect.cpp
intersect_LDADD=libcarve_fileformats.la libcarve_misc.la ../lib/libintersect.la

test_geom_SOURCES=test_geom.cpp
test_geom_LDADD=../lib/libintersect.la

convert_SOURCES=convert.cpp
convert_LDADD=libcarve_fileformats.la libcarve_misc.la ../lib/libintersect.la

view_SOURCES=view.cpp
view_LDADD=libcarve_fileformats.la libcarve_ui.la libcarve_misc.la ../lib/libintersect.la @GL_LIBS@ @GLUT_LIBS@

test_interpolate_SOURCES=test_interpolate.cpp
test_interpolate_LDADD=libcarve_ui.la ../lib/libintersect.la @GL_LIBS@ @GLUT_LIBS@

test_csg_interpolate_SOURCES=test_csg_interpolate.cpp
test_csg_interpolate_LDADD=libcarve_ui.la libcarve_misc.la ../lib/libintersect.la @GL_LIBS@ @GLUT_LIBS@

test_slice_SOURCES=test_slice.cpp
test_slice_LDADD=libcarve_ui.la libcarve_misc.la ../lib/libintersect.la @GL_LIBS@ @GLUT_LIBS@

test_slice_classify_SOURCES=test_slice_classify.cpp
test_slice_classify_LDADD=libcarve_ui.la libcarve_misc.la ../lib/libintersect.la @GL_LIBS@ @GLUT_LIBS@

texture_example_SOURCES=texture_example.cpp
texture_example_LDADD=libcarve_fileformats.la libcarve_ui.la libcarve_misc.la ../lib/libintersect.la @GL_LIBS@ @GLUT_LIBS@

problem_SOURCES=problem.cpp
problem_LDADD=libcarve_ui.la ../lib/libintersect.la @GL_LIBS@ @GLUT_LIBS@

test_rescale_SOURCES=test_rescale.cpp
test_rescale_LDADD=../lib/libintersect.la

tetrahedron_SOURCES=tetrahedron.cpp
tetrahedron_LDADD=../lib/libintersect.la
