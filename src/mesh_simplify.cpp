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
#include "write_ply.hpp"

#include "opts.hpp"

#include <fstream>
#include <string>
#include <utility>
#include <set>
#include <algorithm>

#include <time.h>
#include <sys/time.h>

typedef carve::geom::RTreeNode<3, carve::mesh::Face<3> *> face_rtree_t;

void simplify(carve::mesh::Mesh<3> *mesh,
              const double min_edgelen) {
  face_rtree_t *tree = face_rtree_t::construct_STR(mesh->faces.begin(), mesh->faces.end(), 4, 4);

  for (size_t i = 0; i < mesh->faces.size(); ++i) {
    std::cerr << "i=" << i << " mesh->faces.size()=" << mesh->faces.size() << std::endl;
    carve::mesh::Mesh<3>::face_t *f = mesh->faces[i];
    carve::mesh::Mesh<3>::edge_t *e = f->edge;

    if (e != NULL) {
      std::cerr << "in" << " f->nEdges(): " << f->nEdges() << std::endl;
      do {
        if (e->length2() < min_edgelen * min_edgelen) {
          std::cerr << "remove: " << e << std::endl;
          e = e->collapse();
        } else {
          std::cerr << "skip" << std::endl;
          e = e->next;
        }
      } while (e != f->edge);
    }
  }

  for (size_t i = 0; i < mesh->faces.size(); ) {
    carve::mesh::Mesh<3>::face_t *f = mesh->faces[i];

    if (f->nEdges() < 3) {
      std::cerr << "here." << std::endl;
      if (f->nEdges() == 1) {
        f->edge->collapse();
      } else if (f->nEdges() == 2) {
        carve::mesh::Mesh<3>::edge_t *r1 = f->edge->rev;
        carve::mesh::Mesh<3>::edge_t *r2 = f->edge->next->rev;
        if (r1) r1->rev = r2;
        if (r2) r2->rev = r1;
      }

      delete f;
      mesh->faces[i] = mesh->faces.back();
      mesh->faces.pop_back();
      std::cerr << "done." << std::endl;
    } else {
      ++i;
    }
  }

  delete tree;

  mesh->cacheEdges();
}



int main(int argc, char **argv) {
  carve::input::Input inputs;
  readPLY(std::string(argv[1]), inputs);
  carve::mesh::MeshSet<3> *p;
  p = carve::input::Input::create<carve::mesh::MeshSet<3> >(*inputs.input.begin());

  p->transform(carve::geom::quantize<10,3>());

  for (size_t i = 0; i < p->meshes.size(); ++i) {
    simplify(p->meshes[i], 2e-3);
  }

  carve::poly::Polyhedron *poly = carve::polyhedronFromMesh(p, -1);
  writePLY(std::cout, poly, true);
  return 0;
}
