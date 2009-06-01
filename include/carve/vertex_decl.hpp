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


#pragma once

#include <carve/carve.hpp>

#include <carve/geom2d.hpp>
#include <carve/vector.hpp>
#include <carve/matrix.hpp>
#include <carve/geom3d.hpp>
#include <carve/aabb.hpp>
#include <carve/tag.hpp>

#include <vector>
#include <list>
#include <map>

namespace carve {
  namespace poly {



    struct Polyhedron;
    class Edge;



    class Vertex : public tagable {
    public:
      carve::geom3d::Vector v;
      Polyhedron *owner;

      Vertex();
      ~Vertex();
      Vertex(const carve::geom3d::Vector &);
    };



    struct hash_vertex_ptr {
      size_t operator()(const Vertex * const &v) const {
        return (size_t)v;
      }

      size_t operator()(const std::pair<const Vertex *, const Vertex *> &v) const {
        size_t r = (size_t)v.first;
        size_t s = (size_t)v.second;
        return r ^ ((s >> 16) | (s << 16));
      }

    };



    struct vec_adapt_vertex_ref {
      const carve::geom3d::Vector &operator()(const Vertex &v) { return v.v; }
      carve::geom3d::Vector &operator()(Vertex &v) { return v.v; }
    };



    struct vec_adapt_vertex_ptr {
      const carve::geom3d::Vector &operator()(const Vertex * const &v) { return v->v; }
      carve::geom3d::Vector &operator()(Vertex *&v) { return v->v; }
    };



  }
}
