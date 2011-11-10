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
#include <carve/poly.hpp>
#include <carve/polyline.hpp>
#include <carve/pointset.hpp>
#include <carve/rtree.hpp>
#include <carve/mesh_ops.hpp>
#include <carve/mesh_simplify.hpp>
#include <carve/geom2d.hpp>
#include <carve/heap.hpp>

#include "read_ply.hpp"
#include "write_ply.hpp"

#include "opts.hpp"

#include <fstream>
#include <string>
#include <utility>
#include <set>
#include <algorithm>

#include <time.h>
#include <sys/time.h>

typedef carve::mesh::MeshSet<3> meshset_t;
typedef carve::mesh::Mesh<3> mesh_t;
typedef mesh_t::vertex_t vertex_t;
typedef mesh_t::edge_t edge_t;
typedef mesh_t::face_t face_t;

int main(int argc, char **argv) {
  try {
    carve::input::Input inputs;
    readPLY(std::string(argv[1]), inputs);
    carve::mesh::MeshSet<3> *p;
    p = carve::input::Input::create<carve::mesh::MeshSet<3> >(*inputs.input.begin());

    carve::mesh::MeshSimplifier simplifier;

    simplifier.mergeCoplanarFaces(p, 1e-2);
    // simplifier.snap(p, -5, 0, 0);

    writePLY(std::cout, p, true);
    return 0;
  } catch (carve::exception e) {
    std::cerr << "exception: " << e.str() << std::endl;
  }
}
