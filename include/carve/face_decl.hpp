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



    class Face : public tagable {
    public:
      std::vector<const Vertex *> vertices; // pointer into polyhedron.vertices
      std::vector<const Edge *> edges; // pointer into polyhedron.edges

      Polyhedron *owner;

      carve::geom3d::AABB aabb;
      carve::geom3d::Plane plane_eqn;
      carve::geom3d::Vector centroid;
      int manifold_id;

      carve::geom2d::P2 (*project)(const carve::geom3d::Vector &);
      carve::geom3d::Vector (*unproject)(const carve::geom2d::P2 &, const carve::geom3d::Plane &);

      Face(const std::vector<const Vertex *> &_vertices, bool delay_recalc = false);
      Face(const Vertex *v1, const Vertex *v2, const Vertex *v3, bool delay_recalc = false);
      Face(const Vertex *v1, const Vertex *v2, const Vertex *v3, const Vertex *v4, bool delay_recalc = false);

      Face(const Face *base, const std::vector<const Vertex *> &_vertices, bool flipped) {
        init(base, _vertices, flipped);
      }

      Face() {}
      ~Face() {}

      bool recalc();

      Face *init(const Face *base, const std::vector<const Vertex *> &_vertices, bool flipped);
      Face *create(const std::vector<const Vertex *> &_vertices, bool flipped) const {
        return (new Face)->init(this, _vertices, flipped);
      }
      Face *clone() const {
        return (new Face)->init(this, vertices, false);
      }
      void invert();

      bool containsPoint(const carve::geom3d::Vector &p) const;
      bool simpleLineSegmentIntersection(const carve::geom3d::LineSegment &line,
                                         carve::geom3d::Vector &intersection) const;
      IntersectionClass lineSegmentIntersection(const carve::geom3d::LineSegment &line,
                                                carve::geom3d::Vector &intersection) const;

    };



    struct hash_face_ptr {
      size_t operator()(const Face * const &f) const {
        return (size_t)f;
      }
    };



    typedef std::vector<Face *> FacePtrVector;



    namespace face {



      static inline carve::geom2d::P2 project(const Face *f, const carve::geom3d::Vector &v) {
        return f->project(v);
      }



      static inline carve::geom2d::P2 project(const Face &f, const carve::geom3d::Vector &v) {
        return f.project(v);
      }



      static inline carve::geom3d::Vector unproject(const Face *f, const carve::geom2d::P2 &p) {
        return f->unproject(p, f->plane_eqn);
      }



      static inline carve::geom3d::Vector unproject(const Face &f, const carve::geom2d::P2 &p) {
        return f.unproject(p, f.plane_eqn);
      }



    }



    struct p2_adapt_project {
      typedef carve::geom2d::P2 (*proj_t)(const carve::geom3d::Vector &);
      proj_t proj;
      p2_adapt_project(proj_t _proj) : proj(_proj) { }
      carve::geom2d::P2 operator()(const carve::geom3d::Vector &v) const { return proj(v); }
      carve::geom2d::P2 operator()(const carve::geom3d::Vector *v) const { return proj(*v); }
      carve::geom2d::P2 operator()(const Vertex &v) const { return proj(v.v); }
      carve::geom2d::P2 operator()(const Vertex *v) const { return proj(v->v); }
    };



  }
}
