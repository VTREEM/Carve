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

#include <gtest/gtest.h>

#if defined(HAVE_CONFIG_H)
#  include <carve_config.h>
#endif

#include <carve/carve.hpp>
#include <carve/triangle_intersection.hpp>

#include <fstream>

typedef carve::geom::vector<2> vec2;
typedef carve::geom::vector<2> vec3;

std::istream &operator>>(std::istream &in, vec2 &vec) {
  char c = 0;
  in >> vec.x;
  if (in >> c && c != ',') in.setstate(std::ios_base::failbit);
  in >> vec.y;
  return in;
}

std::ostream &operator<<(std::ostream &out, const vec2 &vec) {
  out << vec.x << "," << vec.y;
  return out;
}

TEST(TriangleIntersectionTest, Test3D) {
  vec3 tri_a[3], tri_b[3];
  tri_a[0] = carve::geom::VECTOR(0,   1,    0);
  tri_a[1] = carve::geom::VECTOR(1,   0.5,  0);
  tri_a[2] = carve::geom::VECTOR(0.5, 0.75, 1);

  tri_b[0] = carve::geom::VECTOR(2,   1,    1);
  tri_b[1] = carve::geom::VECTOR(1,   0.5,  1);
  tri_b[2] = carve::geom::VECTOR(1.5, 0.75, 0);
}

TEST(TriangleIntersectionTest, Test2D) {
  std::ifstream in("intersection_2d.txt");

  while(in.good()) {
    char c;
    vec2 t1[3], t2[3];

    in >> c >> t1[0] >> t1[1] >> t1[2] >> t2[0] >> t2[1] >> t2[2];
    if (!in.eof()) {
      ASSERT_EQ(carve::geom::triangle_intersection(t1, t2), (c == 't') ? carve::geom::TR_INT_INT : carve::geom::TR_INT_NONE);
    }
  }
}
