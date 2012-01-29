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
#include <carve/triangulator.hpp>
#include <carve/geom.hpp>

TEST(Triangulate, Test2) {
  std::vector<carve::geom::vector<2> > poly;
  std::vector<carve::triangulate::tri_idx> result;

  poly.push_back(carve::geom::VECTOR(0.0614821249999999985824672,0.0192046249999999994739763));
  poly.push_back(carve::geom::VECTOR(0.0480376250000000004636291,0.0355867500000000003268497));
  poly.push_back(carve::geom::VECTOR(0.0316555000000000030802028,0.0490312499999999984456878));
  poly.push_back(carve::geom::VECTOR(0.0480376250000000004636291,0.0355867500000000003268497));

  carve::triangulate::triangulate(poly, result);
}

TEST(Triangulate, Test1) {
  std::vector<carve::geom::vector<2> > poly;
  std::vector<carve::triangulate::tri_idx> result;

  poly.push_back(carve::geom::VECTOR(-0.0197657499999999985984545,-0.00112325000000000618793905));
  poly.push_back(carve::geom::VECTOR(-0.0197657499999999985984545,-0.0562291249999999978581577));
  poly.push_back(carve::geom::VECTOR(0.000514374999999999966797393,-0.0562291249999999978581577));
  poly.push_back(carve::geom::VECTOR(0.000514374999999999966797393,-0.00112324999999999924904515));

  carve::triangulate::triangulate(poly, result);
}
