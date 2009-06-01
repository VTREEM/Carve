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

#include <assert.h>
#include <list>

namespace carve {
  namespace poly {



    inline void Polyhedron::invert(int m_id) {
      std::vector<bool> selected_manifolds(manifold_is_closed.size(), false);
      if (m_id >=0 && (unsigned)m_id < selected_manifolds.size()) selected_manifolds[m_id] = true;
      invert(selected_manifolds);
    }


    
    inline void Polyhedron::invert() {
      invertAll();
    }



    inline const Face *Polyhedron::connectedFace(const Face *f, const Edge *e) const {
      const std::vector<const Face *> &edge_faces = connectivity.edge_to_face[edgeToIndex_fast(e)];
      for (size_t i = 0; i < (edge_faces.size() & ~1U); i++) {
        if (edge_faces[i] == f) return edge_faces[i^1];
      }
      return NULL;
    }



    inline Polyhedron *Polyhedron::makeCopy(int m_id) const {
      std::vector<bool> selected_manifolds(manifold_is_closed.size(), false);
      if (m_id >=0 && (unsigned)m_id < selected_manifolds.size()) selected_manifolds[m_id] = true;
      return makeCopy(selected_manifolds);
    }



    inline Polyhedron *Polyhedron::make(std::vector<Face> &_faces, bool _recalc) {
      return new Polyhedron(_faces, _recalc);
    }



    inline Polyhedron *Polyhedron::make(std::list<Face> &_faces, bool _recalc) {
      return new Polyhedron(_faces, _recalc);
    }



    template<typename T>
    int Polyhedron::vertexToEdges(const Vertex *v, T result) const {
      std::vector<const Edge *> &e = connectivity.vertex_to_edge[vertexToIndex_fast(v)];
      std::copy(e.begin(), e.end(), result);
      return e.size();
    }



    template<typename T>
    int Polyhedron::edgeManifolds(const Edge *e, T result) const {
      const std::vector<const Face *> &edge_faces = connectivity.edge_to_face[edgeToIndex_fast(e)];

      for (size_t i = 0; i < (edge_faces.size() & ~1U); ++i) {
        const Face *f1 = edge_faces[i];
        const Face *f2 = edge_faces[i+1];
        assert (f1 || f2);
        if (f1)
          *result++ = f1->manifold_id;
        else if (f2)
          *result++ = f2->manifold_id;
      }
      return edge_faces.size() >> 1;
    }



    template<typename T>
    int Polyhedron::vertexManifolds(const Vertex *v, T result) const {
      const std::vector<const Face *> &f = connectivity.vertex_to_face[vertexToIndex_fast(v)];
      std::set<int> em;

      for (size_t i = 0; i < f.size(); ++i) {
        em.insert(f[i]->manifold_id);
      }

      std::copy(em.begin(), em.end(), result);
      return em.size();
    }



    template<typename T>
    int Polyhedron::_faceNeighbourhood(const Face *f, int depth, T *result) const {
      if (depth < 0 || f->is_tagged()) return 0;

      f->tag();
      *(*result)++ = f;

      int r = 1;
      for (size_t i = 0; i < f->edges.size(); ++i) {
        const std::vector<const Face *> &edge_faces = connectivity.edge_to_face[edgeToIndex_fast(f->edges[i])];
        const Face *f2 = connectedFace(f, f->edges[i]);
        if (f2) {
          r += _faceNeighbourhood(f2, depth - 1, (*result));
        }
      }
      return r;
    }



    template<typename T>
    int Polyhedron::faceNeighbourhood(const Face *f, int depth, T result) const {
      tagable::tag_begin();

      return _faceNeighbourhood(f, depth, &result);
    }



    template<typename T>
    int Polyhedron::faceNeighbourhood(const Edge *e, int m_id, int depth, T result) const {
      tagable::tag_begin();

      int r = 0;
      const std::vector<const Face *> &edge_faces = connectivity.edge_to_face[edgeToIndex_fast(e)];
      for (size_t i = 0; i < edge_faces.size(); ++i) {
        Face *f = edge_faces[i];
        if (f && f->manifold_id == m_id) { r += _faceNeighbourhood(f, depth, &result); }
      }
      return r;
    }



    template<typename T>
    int Polyhedron::faceNeighbourhood(const Vertex *v, int m_id, int depth, T result) const {
      tagable::tag_begin();

      int r = 0;
      const std::vector<const Face *> &vertex_faces = connectivity.vertex_to_face[vertexToIndex_fast(v)];
      for (size_t i = 0; i < vertex_faces.size(); ++i) {
        Face *f = vertex_faces[i];
        if (f && f->manifold_id == m_id) { r += _faceNeighbourhood(f, depth, &result); }
      }
      return r;
    }



    template<typename T>
    int Polyhedron::edgeToFaces(const Edge *e, T result) const {
      const std::vector<const Face *> &edge_faces = connectivity.edge_to_face[edgeToIndex_fast(e)];
      int c = 0;
      for (size_t i = 0; i < edge_faces.size(); ++i) {
        if (edge_faces[i] != NULL) { *result++ = edge_faces[i]; ++c; }
      }
      return c;
    }



    template<typename T>
    int Polyhedron::vertexToFaces(const Vertex *v, T result) const {
      const std::vector<const Face *> &vertex_faces = connectivity.vertex_to_face[vertexToIndex_fast(v)];
      int c = 0;
      for (size_t i = 0; i < vertex_faces.size(); ++i) {
        *result++ = vertex_faces[i]; ++c;
      }
      return c;
    }



    template<typename T>
    void Polyhedron::transform(const T &xform) {
      for (size_t i = 0; i < poly_vertices.size(); i++) {
        poly_vertices[i].v = xform(poly_vertices[i].v);
      }
      faceRecalc();
      init();
    }



    inline ptrdiff_t Polyhedron::vertexToIndex_fast(const Vertex *v) const {
      return v - &poly_vertices[0];
    }



    inline ptrdiff_t Polyhedron::vertexToIndex(const Vertex *v) const {
      if (v < &poly_vertices.front() || v > &poly_vertices.back()) return -1;
      return v - &poly_vertices[0];
    }



    inline ptrdiff_t Polyhedron::edgeToIndex_fast(const Edge *e) const {
      return e - &edges[0];
    }



    inline ptrdiff_t Polyhedron::edgeToIndex(const Edge *e) const {
      if (e < &edges.front() || e > &edges.back()) return -1;
      return e - &edges[0];
    }



    inline size_t Polyhedron::manifoldCount() const {
      return manifold_is_closed.size();
    }



    inline bool Polyhedron::hasOpenManifolds() const {
      for (size_t i = 0; i < manifold_is_closed.size(); ++i) {
        if (!manifold_is_closed[i]) return true;
      }
      return false;
    }



    inline std::ostream &operator<<(std::ostream &o, const Polyhedron &p) {
      p.print(o);
      return o;
    }



  }
}
