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

#include <carve/interpolator.hpp>
#include <carve/csg_triangulator.hpp>

#include <carve/csg.hpp>

#include "carve_texture.h"
#include "brick_texture.h"
#include "leaf_texture.h"

#include "scene.hpp"

#if defined(__APPLE__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <fstream>
#include <string>
#include <utility>
#include <set>

#include <time.h>

#if defined(__GNUC__)
#define __stdcall
#endif

#if defined(GLU_TESS_CALLBACK_VARARGS)
  typedef GLvoid (__stdcall *GLUTessCallback)(...);
#else
  typedef void (__stdcall *GLUTessCallback)();
#endif

carve::geom3d::Vector g_translation;
double g_scale;

static inline void glVertex(const carve::geom3d::Vector &v) {
  glVertex3f(g_scale * (v.x + g_translation.x),
             g_scale * (v.y + g_translation.y),
             g_scale * (v.z + g_translation.z));
}

struct tex_t {
  float u;
  float v;

  tex_t() : u(0.0f), v(0.0f) { }
  tex_t(float _u, float _v) : u(_u), v(_v) { }
};

tex_t operator*(double s, const tex_t &t) {
  return tex_t(t.u * s, t.v * s);
}

tex_t &operator+=(tex_t &t1, const tex_t &t2) {
  t1.u += t2.u;
  t1.v += t2.v;
  return t1;
}

struct vt_t {
  double x, y, z;
  float u, v;
};

void __stdcall tess_vertex(vt_t *v, bool *is_textured) {
  if (*is_textured) {
    glTexCoord2f(v->u, v->v);
  }
  glVertex3d(v->x, v->y, v->z);
}

void drawTexturedPolyhedron(carve::mesh::MeshSet<3> *poly,
                            carve::interpolate::FaceVertexAttr<tex_t> &fv_tex,
                            carve::interpolate::FaceAttr<GLuint> &f_tex_num) {
  glEnable(GL_TEXTURE_2D);

  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

  GLUtesselator *tess = gluNewTess();

  gluTessCallback(tess, GLU_TESS_BEGIN, (GLUTessCallback)glBegin);
  gluTessCallback(tess, GLU_TESS_VERTEX_DATA, (GLUTessCallback)tess_vertex);
  gluTessCallback(tess,  GLU_TESS_END, (GLUTessCallback)glEnd);

  for (carve::mesh::MeshSet<3>::face_iter i = poly->faceBegin(); i != poly->faceEnd(); ++i) {
    carve::mesh::MeshSet<3>::face_t *f = *i;
    std::vector<vt_t> vc(f->nVertices());

    bool textured = true;
    for (carve::mesh::MeshSet<3>::face_t::edge_iter_t e = f->begin(); e != f->end(); ++e) {
      vc[e.idx()].x = g_scale * (e->vert->v.x + g_translation.x);
      vc[e.idx()].y = g_scale * (e->vert->v.y + g_translation.y);
      vc[e.idx()].z = g_scale * (e->vert->v.z + g_translation.z);

      if (fv_tex.hasAttribute(f, e.idx())) {
        tex_t t = fv_tex.getAttribute(f, e.idx());
        vc[e.idx()].u = t.u;
        vc[e.idx()].v = t.v;
      } else {
        textured = false;
      }
    }

    if (textured) {
      GLuint tex_num = f_tex_num.getAttribute(f, 0);
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, tex_num);
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    } else {
      glColor4f(0.5f, 0.6f, 0.7f, 1.0f);
    }

    glNormal3dv(f->plane.N.v);

    gluTessBeginPolygon(tess, (void *)&textured);
    gluTessBeginContour(tess);

    for (size_t j = 0; j != vc.size(); ++j) {
      gluTessVertex(tess, (GLdouble *)&vc[j], (GLvoid *)&vc[j]);
    }

    gluTessEndContour(tess);
    gluTessEndPolygon(tess);

  }

  gluDeleteTess(tess);

  glDisable(GL_TEXTURE_2D);
}

void drawWireframePolyhedron(carve::mesh::MeshSet<3> *poly) {
  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
  glDisable(GL_LIGHTING);

  for (carve::mesh::MeshSet<3>::face_iter i = poly->faceBegin(); i != poly->faceEnd(); ++i) {
    carve::mesh::MeshSet<3>::face_t *f = *i;

    glBegin(GL_LINE_LOOP);
    for (carve::mesh::MeshSet<3>::face_t::edge_iter_t e = f->begin(); e != f->end(); ++e) {
      double x = g_scale * (e->vert->v.x + g_translation.x);
      double y = g_scale * (e->vert->v.y + g_translation.y);
      double z = g_scale * (e->vert->v.z + g_translation.z);
      glVertex3d(x, y, z);
    }
    glEnd();
  }

  glEnable(GL_LIGHTING);
}

carve::mesh::MeshSet<3> *texturedCube(
    carve::interpolate::FaceVertexAttr<tex_t> &fv_tex,
    carve::interpolate::FaceAttr<GLuint> &f_tex_num,
    GLuint tex,
    const carve::math::Matrix &transform = carve::math::Matrix::IDENT()) {

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

  for (size_t i = 0; i < 6; ++i) {
    fv_tex.setAttribute(faces[i], 0, tex_t(0.0f, 1.0f));
    fv_tex.setAttribute(faces[i], 1, tex_t(1.0f, 1.0f));
    fv_tex.setAttribute(faces[i], 2, tex_t(1.0f, 0.0f));
    fv_tex.setAttribute(faces[i], 3, tex_t(0.0f, 0.0f));
    f_tex_num.setAttribute(faces[i], tex);
  }

  carve::mesh::MeshSet<3> *poly = new carve::mesh::MeshSet<3>(faces);

  return poly;
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
    for (size_t i = 0; i < draw_flags.size(); ++i) {
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

GLuint initTexture(GLuint w, GLuint h, const unsigned char *data) {
   GLuint tex;
   glGenTextures(1, &tex);
   glBindTexture(GL_TEXTURE_2D, tex);
   glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
   glTexImage2D(GL_TEXTURE_2D,
                0,
                GL_RGB,
                w,
                h,
                0,
                GL_RGB,
                GL_UNSIGNED_BYTE,
                data);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

   return tex;
}

void destroyTexture(GLuint tex) {
  glDeleteTextures(1, &tex);
}

int main(int argc, char **argv) {
  TestScene *scene = new TestScene(argc, argv, 2);

  GLuint tex_1 = initTexture(128, 128, carve_texture);
  GLuint tex_2 = initTexture(128, 128, brick_texture);
  GLuint tex_3 = initTexture(128, 128, leaf_texture);

  g_scale = 10.0;

  carve::interpolate::FaceVertexAttr<tex_t> fv_tex;
  carve::interpolate::FaceAttr<GLuint> f_tex_num;
  carve::mesh::MeshSet<3> *base = NULL;

  bool b = true;
  for (int x = -10; x <= +10; x += 5) {
    for (int y = -10; y <= +10; y += 5) {
      for (int z = -10; z <= +10; z += 5) {
        double rot = x * .17 + y * .06 + z * .09;
        carve::mesh::MeshSet<3> *r = texturedCube(fv_tex, f_tex_num, b ? tex_2 : tex_3,
                                 carve::math::Matrix::TRANS(x/2.5, y/2.5, z/2.5) *
                                 carve::math::Matrix::ROT(rot, 1,2,3));
        b = !b;
        if (base) {
          carve::mesh::MeshSet<3> *temp = base;
          carve::csg::CSG csg;
          fv_tex.installHooks(csg);
          f_tex_num.installHooks(csg);

          base = csg.compute(temp, r, carve::csg::CSG::UNION);
          delete temp;
          delete r;
        } else {
          base = r;
        }
      }
    }
  }

  carve::mesh::MeshSet<3> *r1 = texturedCube(fv_tex, f_tex_num, tex_1,
                            carve::math::Matrix::TRANS(0,0,4) *
                            carve::math::Matrix::SCALE(4,4,4));

  carve::mesh::MeshSet<3> *r2 = texturedCube(fv_tex, f_tex_num, tex_1,
                            carve::math::Matrix::TRANS(0,0,5) *
                            carve::math::Matrix::SCALE(2, 2, 2));

  carve::csg::CSG csg;
  fv_tex.installHooks(csg);
  f_tex_num.installHooks(csg);

  carve::mesh::MeshSet<3> *r3 = csg.compute(base, r1, carve::csg::CSG::INTERSECTION, NULL, carve::csg::CSG::CLASSIFY_EDGE);
  carve::mesh::MeshSet<3> *r4 = csg.compute(r3, r2, carve::csg::CSG::UNION, NULL, carve::csg::CSG::CLASSIFY_EDGE);

  glNewList(scene->draw_list_base, GL_COMPILE);
  drawTexturedPolyhedron(r4, fv_tex, f_tex_num);
  glEndList();

  glNewList(scene->draw_list_base+1, GL_COMPILE);
  drawWireframePolyhedron(r3);
  glEndList();

  scene->draw_flags[0] = true;
  scene->draw_flags[1] = true;

  scene->run();

  destroyTexture(tex_1);
  destroyTexture(tex_2);
  destroyTexture(tex_3);

  delete scene;

  return 0;
}
