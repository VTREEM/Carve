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

#include <carve/csg.hpp>
#include <carve/triangulator.hpp>

namespace carve {
  namespace csg {

    class CarveTriangulator : public carve::csg::CSG::Hook {

    public:
      CarveTriangulator() {
      }
      virtual ~CarveTriangulator() {
      }
      virtual void processOutputFace(std::vector<carve::poly::Face<3> *> &faces,
                                     const carve::poly::Face<3> *orig,
                                     bool flipped) {
        size_t f = 0;
        while (f < faces.size()) {
          carve::poly::Face<3> *face = faces[f];
          if (face->vertices.size() == 3) {
            ++f;
            continue;
          }
          std::vector<carve::geom2d::P2> projected;
          projected.reserve(face->vertices.size());
          for (size_t i = 0; i < face->vertices.size(); ++i) {
            projected.push_back(face->project(face->vertices[i]->v));
          }
          std::vector<carve::triangulate::tri_idx> result;
          carve::triangulate::triangulate(projected, result);
          
          faces.erase(faces.begin() + f);
          faces.insert(faces.begin() + f, result.size(), NULL);
          std::vector<const carve::poly::Vertex<3> *> fv;
          fv.resize(3);
          for (size_t i = 0; i < result.size(); ++i) {
            fv[0] = face->vertices[result[i].a];
            fv[1] = face->vertices[result[i].b];
            fv[2] = face->vertices[result[i].c];
            faces[f++] = face->create(fv, false);
          }
        }
      }
    };

  }
}
