// Begin License:
// Copyright (C) 2006-2011 Tobias Sargeant (tobias.sargeant@gmail.com).
// All rights reserved.
//
// This file is part of the Carve CSG Library (http://carve-csg.com/)
//
// This file may be used under the terms of the GNU General Public
// License version 2.0 as published by the Free Software Foundation
// and appearing in the file LICENSE.GPL2 included in the packaging of
// this file.
//
// This file is provided "AS IS" with NO WARRANTY OF ANY KIND,
// INCLUDING THE WARRANTIES OF DESIGN, MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE.
// End:


#if defined(HAVE_CONFIG_H)
#  include <carve_config.h>
#endif

#include <carve/carve.hpp>
#include <carve/interpolator.hpp>
#include <carve/csg.hpp>

#include <fstream>
#include <string>
#include <utility>
#include <set>

#include <time.h>

#include "geom_draw.hpp"
#include "geometry.hpp"
#include "scene.hpp"
#include "rgb.hpp"

#if defined(__APPLE__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

static inline void glVertex(const carve::geom3d::Vector &v) {
  glVertex3f(g_scale * (v.x + g_translation.x),
             g_scale * (v.y + g_translation.y),
             g_scale * (v.z + g_translation.z));
}

static inline void glColor(const cRGB &c, float alpha) {
  glColor4f(c.r, c.g, c.b, alpha);
}

carve::interpolate::FaceEdgeAttr<cRGB> fe_edgeflag;

void drawMeshSetEdgeFlags(const carve::mesh::MeshSet<3> *meshset, double alpha) {
  glBegin(GL_LINES);

  for (carve::mesh::MeshSet<3>::const_face_iter i = meshset->faceBegin(); i != meshset->faceEnd(); ++i) {
    const carve::mesh::MeshSet<3>::face_t *f = *i;
    for (carve::mesh::MeshSet<3>::face_t::const_edge_iter_t e = f->begin(); e != f->end(); ++e) {
      if (fe_edgeflag.hasAttribute(f, e.idx())) {
        glColor(fe_edgeflag.getAttribute(f, e.idx()), alpha);
        glVertex(e->v1()->v);
        glVertex(e->v2()->v);
      }
    }
  }

  glEnd();
}

void drawMeshSetEdgeFlags(const carve::mesh::MeshSet<3> *meshset) {
  glDisable(GL_LIGHTING);
  glDepthMask(GL_FALSE);
  glDisable(GL_DEPTH_TEST);

  drawMeshSetEdgeFlags(meshset, 0.2);

  glEnable(GL_DEPTH_TEST);

  drawMeshSetEdgeFlags(meshset, 0.8);

  glDepthMask(GL_TRUE);
  glEnable(GL_LIGHTING);
}



carve::mesh::MeshSet<3> *edgeFlagCube(const carve::math::Matrix &transform = carve::math::Matrix::IDENT()) {
  std::vector<carve::mesh::MeshSet<3>::vertex_t> v;
  v.push_back(carve::mesh::MeshSet<3>::vertex_t(transform * carve::geom::VECTOR(+1.0, +1.0, +1.0)));
  v.push_back(carve::mesh::MeshSet<3>::vertex_t(transform * carve::geom::VECTOR(-1.0, +1.0, +1.0)));
  v.push_back(carve::mesh::MeshSet<3>::vertex_t(transform * carve::geom::VECTOR(-1.0, -1.0, +1.0)));
  v.push_back(carve::mesh::MeshSet<3>::vertex_t(transform * carve::geom::VECTOR(+1.0, -1.0, +1.0)));
  v.push_back(carve::mesh::MeshSet<3>::vertex_t(transform * carve::geom::VECTOR(+1.0, +1.0, -1.0)));
  v.push_back(carve::mesh::MeshSet<3>::vertex_t(transform * carve::geom::VECTOR(-1.0, +1.0, -1.0)));
  v.push_back(carve::mesh::MeshSet<3>::vertex_t(transform * carve::geom::VECTOR(-1.0, -1.0, -1.0)));
  v.push_back(carve::mesh::MeshSet<3>::vertex_t(transform * carve::geom::VECTOR(+1.0, -1.0, -1.0)));

  std::vector<carve::mesh::MeshSet<3>::face_t *> faces;
  faces.reserve(6);

  faces.push_back(new carve::mesh::MeshSet<3>::face_t(&v[0], &v[1], &v[2], &v[3]));
  faces.push_back(new carve::mesh::MeshSet<3>::face_t(&v[7], &v[6], &v[5], &v[4]));
  faces.push_back(new carve::mesh::MeshSet<3>::face_t(&v[0], &v[4], &v[5], &v[1]));
  faces.push_back(new carve::mesh::MeshSet<3>::face_t(&v[1], &v[5], &v[6], &v[2]));
  faces.push_back(new carve::mesh::MeshSet<3>::face_t(&v[2], &v[6], &v[7], &v[3]));
  faces.push_back(new carve::mesh::MeshSet<3>::face_t(&v[3], &v[7], &v[4], &v[0]));

  fe_edgeflag.setAttribute(faces[0], 0, cRGB(1,1,1));
  fe_edgeflag.setAttribute(faces[0], 1, cRGB(1,1,1));
  fe_edgeflag.setAttribute(faces[0], 2, cRGB(1,1,1));
  fe_edgeflag.setAttribute(faces[0], 3, cRGB(1,1,1));

  fe_edgeflag.setAttribute(faces[1], 0, cRGB(0,0,0));
  fe_edgeflag.setAttribute(faces[1], 1, cRGB(0,0,0));
  fe_edgeflag.setAttribute(faces[1], 2, cRGB(0,0,0));
  fe_edgeflag.setAttribute(faces[1], 3, cRGB(0,0,0));



  return new carve::mesh::MeshSet<3>(faces);
}

struct TestScene : public Scene {
  GLuint draw_list_base;
  std::vector<bool> draw_flags;

  virtual bool key(unsigned char k, int x, int y) {
    const char *t;
    static const char *l = "1234567890!@#$%^&*()";
    t = strchr(l, k);
    if (t != NULL) {
      int layer = t - l;
      if (layer < draw_flags.size()) {
        draw_flags[layer] = !draw_flags[layer];
      }
    }
    return true;
  }

  virtual GLvoid draw() {
    for (int i = 0; i < draw_flags.size(); ++i) {
      if (draw_flags[i]) glCallList(draw_list_base + i);
    }
  }

  TestScene(int argc, char **argv, int n_dlist) : Scene(argc, argv) {
    draw_list_base = glGenLists(n_dlist);

    draw_flags.resize(n_dlist, false);
  }

  virtual ~TestScene() {
    glDeleteLists(draw_list_base, draw_flags.size());
  }
};

int main(int argc, char **argv) {
  installDebugHooks();

  TestScene *scene = new TestScene(argc, argv, 4);

  g_scale = 10.0;

  glNewList(scene->draw_list_base, GL_COMPILE);
  carve::mesh::MeshSet<3> *a = edgeFlagCube(carve::math::Matrix::ROT(1.6, -.2, .3, .4));
  carve::mesh::MeshSet<3> *b = makeTorus(20, 20, 1.4, .5);
  carve::csg::CSG csg;
  fe_edgeflag.installHooks(csg);
  carve::mesh::MeshSet<3> *c = csg.compute(a, b, carve::csg::CSG::A_MINUS_B);
  glEndList();

  glNewList(scene->draw_list_base + 1, GL_COMPILE);
  drawMeshSet(a, .6, .6, .6, 1.0);
  drawMeshSetEdgeFlags(a);
  glEndList();

  glNewList(scene->draw_list_base + 2, GL_COMPILE);
  drawMeshSet(b, .6, .6, .6, 1.0);
  drawMeshSetEdgeFlags(b);
  glEndList();

  glNewList(scene->draw_list_base + 3, GL_COMPILE);
  drawMeshSet(c, .6, .6, .6, 1.0);
  drawMeshSetEdgeFlags(c);
  glEndList();

  scene->draw_flags[0] = false;
  scene->draw_flags[1] = false;
  scene->draw_flags[2] = false;
  scene->draw_flags[3] = true;

  scene->run();

  delete scene;

  return 0;
}
