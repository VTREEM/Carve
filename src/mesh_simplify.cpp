// Begin License:
// Copyright (C) 2006-2008 Tobias Sargeant (tobias.sargeant@gmail.com).
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
#include <carve/poly.hpp>
#include <carve/polyline.hpp>
#include <carve/pointset.hpp>
#include <carve/rtree.hpp>

#include "read_ply.hpp"

#include "opts.hpp"

#include <fstream>
#include <string>
#include <utility>
#include <set>
#include <algorithm>

#include <time.h>
#include <sys/time.h>

int main(int argc, char **argv) {
  carve::input::Input inputs;
  readPLY(std::string(argv[1]), inputs);
  carve::mesh::MeshSet<3> *p;
  p = carve::input::Input::create<carve::mesh::MeshSet<3> >(*inputs.input.begin());

  typedef carve::geom::RTreeNode<3, carve::mesh::Face<3> *> face_rtree_t;

  face_rtree_t *tree = face_rtree_t::construct_STR(p->faceBegin(), p->faceEnd(), 50, 4);

  return 0;
}
