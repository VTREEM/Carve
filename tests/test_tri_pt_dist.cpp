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

#include <carve/csg.hpp>
#include <carve/csg_triangulator.hpp>

#include "scene.hpp"
#include "geom_draw.hpp"
#include "geometry.hpp"
#include <carve/input.hpp>

#include <gloop/gl.hpp>
#include <gloop/glu.hpp>
#include <gloop/glut.hpp>

#include <fstream>
#include <string>
#include <utility>
#include <set>

#include <time.h>

#include BOOST_INCLUDE(random.hpp)

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

uint32_t getseed() {
#if defined(__APPLE__)
    srandomdev();
#endif
      return random();
}

boost::mt19937 rng(getseed());
boost::uniform_on_sphere<double> distrib(3);
boost::variate_generator<boost::mt19937 &, boost::uniform_on_sphere<double> > gen(rng, distrib);

carve::geom::vector<3> randomUnitVector() {
  carve::geom::vector<3> vec;
  vec = gen();
  return vec;
}

int main(int argc, char **argv) {
  TestScene *scene = new TestScene(argc, argv, 1);

  glNewList(scene->draw_list_base + 0, GL_COMPILE);

  glEnable(GL_DEPTH_TEST);
  glDepthMask(GL_TRUE);
  glEnable(GL_CULL_FACE);
  glColor4f(.7, .7, .7, 1.0);
  glEnable(GL_LIGHTING);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  carve::geom::tri<3> tri(
    randomUnitVector() * 10.0,
    randomUnitVector() * 10.0,
    randomUnitVector() * 10.0);

  carve::geom::vector<3> p = randomUnitVector() * 20.0;
  carve::geom::vector<3> tp = carve::geom::closestPoint(tri, p);
  double r = carve::geom::distance(p, tp);
  carve::geom::sphere<3> sphere(p, r);

  drawTri(tri);
  drawSphere(sphere);

  glEndList();

  scene->run();

  delete scene;

  return 0;
}
