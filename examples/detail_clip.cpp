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
  std::list<carve::mesh::MeshSet<3>::face_t *> faces;
  std::set<const carve::mesh::MeshSet<3>::face_t *> seen;
  const carve::mesh::MeshSet<3> *src_a;
  const carve::mesh::MeshSet<3> *src_b;
  
  DetailClip(const carve::mesh::MeshSet<3> *_src_a,
             const carve::mesh::MeshSet<3> *_src_b) : carve::csg::CSG::Collector(), src_a(_src_a), src_b(_src_b) {
  }

  virtual ~DetailClip() {
  }

  virtual void collect(carve::csg::FaceLoopGroup *grp, carve::csg::CSG::Hooks &hooks) {
    if (grp->face_loops.head->orig_face->mesh->meshset != src_b) return;
    if (grp->classificationAgainst(NULL) == carve::csg::FACE_IN) return;
    if (grp->classificationAgainst(src_a->meshes[0]) == carve::csg::FACE_IN) return;

    for (carve::csg::FaceLoop *f = grp->face_loops.head; f; f = f->next) {
      if (seen.find(f->orig_face) == seen.end()) {
        std::vector<carve::mesh::MeshSet<3>::vertex_t *> vertices;
        f->orig_face->getVertices(vertices);
        faces.push_back(f->orig_face->create(vertices.begin(), vertices.end(), false));
        seen.insert(f->orig_face);
      }
    }
  }

  virtual carve::mesh::MeshSet<3> *done(carve::csg::CSG::Hooks &hooks) {
    return new carve::mesh::MeshSet<3>(faces);
  }
};


int main(int argc, char **argv) {
  carve::mesh::MeshSet<3> *a, *b;
  a = readPLYasMesh(argv[1]);
  b = readPLYasMesh(argv[2]);
  DetailClip detail_clip_collector(a, b);
  carve::mesh::MeshSet<3> *c = carve::csg::CSG().compute(a, b, detail_clip_collector, NULL, carve::csg::CSG::CLASSIFY_EDGE);
  writePLY(std::cout, c, false);
  return 0;
}
