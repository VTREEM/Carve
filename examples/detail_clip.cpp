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

#include <iostream>
#include <carve/csg.hpp>

#include "read_ply.hpp"
#include "write_ply.hpp"

#include <carve/input.hpp>

#include <fstream>
#include <string>
#include <utility>
#include <set>

#include <time.h>

class DetailClip : public carve::csg::CSG::Collector {
  DetailClip();
  DetailClip(const DetailClip &);
  DetailClip &operator=(const DetailClip &);

public:
  std::list<carve::poly::Face<3> > faces;
  std::set<const carve::poly::Face<3> *> seen;
  const carve::poly::Polyhedron *src_a;
  const carve::poly::Polyhedron *src_b;
  
  DetailClip(const carve::poly::Polyhedron *_src_a,
             const carve::poly::Polyhedron *_src_b) : carve::csg::CSG::Collector(), src_a(_src_a), src_b(_src_b) {
  }

  virtual ~DetailClip() {
  }

  virtual void collect(carve::csg::FaceLoopGroup *grp, carve::csg::CSG::Hooks &hooks) {
    if (grp->face_loops.head->orig_face->owner != src_b) return;
    if (grp->classificationAgainst(src_a, -1) == carve::csg::FACE_IN) return;
    if (grp->classificationAgainst(src_a, 0) == carve::csg::FACE_IN) return;

    for (carve::csg::FaceLoop *f = grp->face_loops.head; f; f = f->next) {
      if (seen.find(f->orig_face) == seen.end()) {
        faces.push_back(carve::poly::Face<3>());
        faces.back().init(f->orig_face, f->orig_face->vbegin(), f->orig_face->vend(), false);
        seen.insert(f->orig_face);
      }
    }
  }

  virtual carve::poly::Polyhedron *done(carve::csg::CSG::Hooks &hooks) {
    return new carve::poly::Polyhedron(faces);
  }
};


int main(int argc, char **argv) {
  carve::poly::Polyhedron *a, *b;
  a = readPLY(argv[1]);
  b = readPLY(argv[2]);
  DetailClip detail_clip_collector(a, b);
  carve::poly::Polyhedron *c = carve::csg::CSG().compute(a, b, detail_clip_collector, NULL, carve::csg::CSG::CLASSIFY_EDGE);
  writePLY(std::cout, c, false);
  return 0;
}
