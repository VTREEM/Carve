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

#include <carve/poly.hpp>

double CALC_X(const carve::geom::plane<3> &p, double y, double z) { return -(p.d + p.N.y * y + p.N.z * z) / p.N.x; }
double CALC_Y(const carve::geom::plane<3> &p, double x, double z) { return -(p.d + p.N.x * x + p.N.z * z) / p.N.y; }
double CALC_Z(const carve::geom::plane<3> &p, double x, double y) { return -(p.d + p.N.x * x + p.N.y * y) / p.N.z; }

namespace carve {
  namespace poly {

    carve::geom2d::P2 _project_1(const carve::geom3d::Vector &v) {
      return carve::geom::VECTOR(v.z, v.y);
    }

    carve::geom2d::P2 _project_2(const carve::geom3d::Vector &v) {
      return carve::geom::VECTOR(v.x, v.z);
    }

    carve::geom2d::P2 _project_3(const carve::geom3d::Vector &v) {
      return carve::geom::VECTOR(v.y, v.x);
    }

    carve::geom2d::P2 _project_4(const carve::geom3d::Vector &v) {
      return carve::geom::VECTOR(v.y, v.z);
    }

    carve::geom2d::P2 _project_5(const carve::geom3d::Vector &v) {
      return carve::geom::VECTOR(v.z, v.x);
    }

    carve::geom2d::P2 _project_6(const carve::geom3d::Vector &v) {
      return carve::geom::VECTOR(v.x, v.y);
    }


    carve::geom3d::Vector _unproject_1(const carve::geom2d::P2 &p, const carve::geom3d::Plane &plane_eqn) {
      return carve::geom::VECTOR(CALC_X(plane_eqn, p.y, p.x), p.y, p.x);
    }

    carve::geom3d::Vector _unproject_2(const carve::geom2d::P2 &p, const carve::geom3d::Plane &plane_eqn) {
      return carve::geom::VECTOR(p.x, CALC_Y(plane_eqn, p.x, p.y), p.y);
    }

    carve::geom3d::Vector _unproject_3(const carve::geom2d::P2 &p, const carve::geom3d::Plane &plane_eqn) {
      return carve::geom::VECTOR(p.y, p.x, CALC_Z(plane_eqn, p.y, p.x));
    }

    carve::geom3d::Vector _unproject_4(const carve::geom2d::P2 &p, const carve::geom3d::Plane &plane_eqn) {
      return carve::geom::VECTOR(CALC_X(plane_eqn, p.x, p.y), p.x, p.y);
    }

    carve::geom3d::Vector _unproject_5(const carve::geom2d::P2 &p, const carve::geom3d::Plane &plane_eqn) {
      return carve::geom::VECTOR(p.y, CALC_Y(plane_eqn, p.y, p.x), p.x);
    }

    carve::geom3d::Vector _unproject_6(const carve::geom2d::P2 &p, const carve::geom3d::Plane &plane_eqn) {
      return carve::geom::VECTOR(p.x, p.y, CALC_Z(plane_eqn, p.x, p.y));
    }

    Face::Face(const std::vector<const carve::poly::Vertex *> &_vertices,
               bool delay_recalc) : tagable() {
      vertices = _vertices;
      if (!delay_recalc && !recalc()) { }
    }

    Face::Face(const carve::poly::Vertex *a,
               const carve::poly::Vertex *b,
               const carve::poly::Vertex *c,
               bool delay_recalc) : tagable() {
      vertices.reserve(3);
      vertices.push_back(a);
      vertices.push_back(b);
      vertices.push_back(c);
      if (!delay_recalc && !recalc()) { }
    }

    Face::Face(const carve::poly::Vertex *a,
               const carve::poly::Vertex *b,
               const carve::poly::Vertex *c,
               const carve::poly::Vertex *d,
               bool delay_recalc) : tagable() {
      vertices.reserve(4);
      vertices.push_back(a);
      vertices.push_back(b);
      vertices.push_back(c);
      vertices.push_back(d);
      if (!delay_recalc && !recalc()) { }
    }

    static carve::geom2d::P2 (*project_tab[2][3])(const carve::geom3d::Vector &) = {
      { &_project_1, &_project_2, &_project_3 },
      { &_project_4, &_project_5, &_project_6 }
    };
    static carve::geom3d::Vector (*unproject_tab[2][3])(const carve::geom2d::P2 &, const carve::geom3d::Plane &) = {
      { &_unproject_1, &_unproject_2, &_unproject_3 },
      { &_unproject_4, &_unproject_5, &_unproject_6 }
    };

    void Face::invert() {
      size_t n_verts = vertices.size();
      for (size_t i = 0; i < n_verts / 2; ++i) std::swap(vertices[i], vertices[n_verts - 1 - i]);


      if (project != NULL) {
        plane_eqn.negate();

        int da = carve::geom::dominantAxis(plane_eqn.N);

        project = project_tab[plane_eqn.N.v[da] > 0][da];
        unproject = unproject_tab[plane_eqn.N.v[da] > 0][da];
      }

      if (edges.size() == n_verts) {
        for (size_t i = 0; i < (n_verts - 1) / 2; ++i) std::swap(edges[i], edges[n_verts - 2 - i]);
        for (size_t i = 0; i < n_verts; i++) {
          const carve::poly::Vertex *v1 = vertices[i];
          const carve::poly::Vertex *v2 = vertices[(i+1) % n_verts];
          ASSERT((edges[i]->v1 == v1 && edges[i]->v2 == v2) || (edges[i]->v1 == v2 && edges[i]->v2 == v1));
        }
      }
    }

    bool Face::recalc() {
      aabb.fit(vertices.begin(), vertices.end(), vec_adapt_vertex_ptr());

      if (!carve::geom3d::fitPlane(vertices.begin(), vertices.end(), vec_adapt_vertex_ptr(), centroid, plane_eqn)) {
        return false;
      }

      int da = carve::geom::dominantAxis(plane_eqn.N);
      project = project_tab[0][da];

      double A = carve::geom2d::signedArea(vertices, p2_adapt_project(project));
      if ((A < 0.0) ^ (plane_eqn.N.v[da] < 0.0)) {
        plane_eqn.negate();
      }

      project = project_tab[plane_eqn.N.v[da] > 0][da];
      unproject = unproject_tab[plane_eqn.N.v[da] > 0][da];

      return true;
    }

    Face *Face::init(const Face *base, const std::vector<const carve::poly::Vertex *> &_vertices, bool flipped) {
      vertices.reserve(_vertices.size());

      if (flipped) {
        std::copy(_vertices.rbegin(), _vertices.rend(), std::back_insert_iterator<std::vector<const carve::poly::Vertex *> >(vertices));
        plane_eqn = -base->plane_eqn;
      } else {
        std::copy(_vertices.begin(), _vertices.end(), std::back_insert_iterator<std::vector<const carve::poly::Vertex *> >(vertices));
        plane_eqn = base->plane_eqn;
      }

      aabb.fit(vertices.begin(), vertices.end(), vec_adapt_vertex_ptr());
      carve::geom::centroid(vertices.begin(), vertices.end(), vec_adapt_vertex_ptr(), centroid);
      untag();

      int da = carve::geom::dominantAxis(plane_eqn.N);

      project = project_tab[plane_eqn.N.v[da] > 0][da];
      unproject = unproject_tab[plane_eqn.N.v[da] > 0][da];

      return this;
    }

    bool Face::containsPoint(const carve::geom3d::Vector &p) const {
      if (!carve::math::ZERO(carve::geom::distance(plane_eqn, p))) return false;
      // return pointInPolySimple(vertices, p2_adapt_project(project), (this->*project)(p));
      return carve::geom2d::pointInPoly(vertices, p2_adapt_project(project), face::project(this, p)).iclass != POINT_OUT;
    }

    bool Face::simpleLineSegmentIntersection(const carve::geom3d::LineSegment &line,
                                             carve::geom3d::Vector &intersection) const {
      if (!line.OK()) return false;

      carve::geom3d::Vector p;
      IntersectionClass intersects = carve::geom3d::lineSegmentPlaneIntersection(plane_eqn,
                                                                  line,
                                                                  p);
      if (intersects == INTERSECT_NONE || intersects == INTERSECT_BAD) {
        return false;
      }

      carve::geom2d::P2 proj_p(face::project(this, p));
      if (carve::geom2d::pointInPolySimple(vertices, p2_adapt_project(project), proj_p)) {
        intersection = p;
        return true;
      }
      return false;
    }

    // XXX: should try to return a pre-existing vertex in the case of a
    // line-vertex intersection.  as it stands, this code isn't used,
    // so... meh.
    IntersectionClass Face::lineSegmentIntersection(const carve::geom3d::LineSegment &line,
                                                    carve::geom3d::Vector &intersection) const {
      if (!line.OK()) return INTERSECT_NONE;

  
      carve::geom3d::Vector p;
      IntersectionClass intersects = carve::geom3d::lineSegmentPlaneIntersection(plane_eqn,
                                                                  line,
                                                                  p);
      if (intersects == INTERSECT_NONE || intersects == INTERSECT_BAD) {
        return intersects;
      }

      carve::geom2d::P2 proj_p(face::project(this, p));
      carve::geom2d::PolyInclusionInfo pi = carve::geom2d::pointInPoly(vertices, p2_adapt_project(project), proj_p);
      switch (pi.iclass) {
      case POINT_VERTEX:
        intersection = p;
        return INTERSECT_VERTEX;

      case POINT_EDGE:
        intersection = p;
        return INTERSECT_EDGE;

      case POINT_IN:
        intersection = p;
        return INTERSECT_FACE;
      
      case POINT_OUT:
        return INTERSECT_NONE;

      default:
        break;
      }
      return INTERSECT_NONE;
    }

  }
}

