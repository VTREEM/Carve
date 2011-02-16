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

#include <carve/geom2d.hpp>
#include <carve/geom3d.hpp>
#include <carve/djset.hpp>

#include <iostream>
#include <deque>

namespace carve {
  namespace mesh {



    namespace detail {
      template<typename list_t>
      struct circlist_iter_t {
        typedef std::bidirectional_iterator_tag iterator_category;
        typedef list_t                          value_type;
        typedef void                            difference_type;
        typedef value_type &                    reference;
        typedef value_type *                    pointer;

        list_t *curr;

        circlist_iter_t(list_t *pos) : curr(pos) { }
        circlist_iter_t operator++(int) { circlist_iter_t result(*this); curr = curr->next; return result; }
        circlist_iter_t operator--(int) { circlist_iter_t result(*this); curr = curr->prev; return result; }
        circlist_iter_t operator++() { curr = curr->next; return *this; }
        circlist_iter_t operator--() { curr = curr->prev; return *this; }
        bool operator==(const circlist_iter_t &other) const { return curr == other.curr; }
        bool operator!=(const circlist_iter_t &other) const { return curr != other.curr; }

        reference operator*() { return *curr; }
        pointer operator->() { return curr; }
      };



      template<typename list_t>
      struct fwd_circlist_iter_t {
        typedef std::forward_iterator_tag       iterator_category;
        typedef list_t                          value_type;
        typedef void                            difference_type;
        typedef value_type &                    reference;
        typedef value_type *                    pointer;

        list_t *curr, *start;

        void advance() {
          if (!start) start = curr;
          if (curr) curr = curr->next;
          if (curr == start) curr = start = NULL;
        }

        fwd_circlist_iter_t(list_t *pos) : curr(pos), start(NULL) { }
        fwd_circlist_iter_t() : curr(NULL), start(NULL) { } // end iterator
        fwd_circlist_iter_t operator++(int) { fwd_circlist_iter_t result(*this); advance(); return result; }
        fwd_circlist_iter_t operator++() { advance(); return *this; }
        bool operator==(const fwd_circlist_iter_t &other) const { return curr == other.curr; }
        bool operator!=(const fwd_circlist_iter_t &other) const { return curr != other.curr; }

        reference operator*() { return *curr; }
        pointer operator->() { return curr; }
      };



      template<typename iter_t, typename mapping_t>
      struct mapped_iter_t {
        typedef typename std::iterator_traits<iter_t>::iterator_category iterator_category;
        typedef typename std::iterator_traits<iter_t>::difference_type   difference_type;
        typedef typename mapping_t::output_type                          value_type;
        typedef value_type &                                             reference;
        typedef value_type *                                             pointer;

        iter_t base;

        mapped_iter_t(iter_t _base) : base(_base) { }

        bool operator==(const mapped_iter_t &other) const { return base == other.base; }
        bool operator!=(const mapped_iter_t &other) const { return base != other.base; }

        mapped_iter_t operator++(int) { mapped_iter_t result(*this); ++base; return result; }
        mapped_iter_t operator--(int) { mapped_iter_t result(*this); --base; return result; }
        mapped_iter_t operator++() { ++base; return *this; }
        mapped_iter_t operator--() { --base; return *this; }

        value_type operator*() { return mapping_t()(*base); }
      };


      template<unsigned ndim>
      struct edge_vertex_mapping {
        typedef typename Edge<ndim>::vertex_t output_type;
        output_type &operator()(Edge<ndim> &edge) { return *(edge.vert); }
      };



      template<unsigned ndim>
      struct const_edge_vertex_mapping {
        typedef const typename Edge<ndim>::vertex_t output_type;
        output_type &operator()(const Edge<ndim> &edge) { return *(edge.vert); }
      };
    }



    template<unsigned ndim>
    void Edge<ndim>::remove() {
      if (rev) { rev->rev = NULL; rev = NULL; }
      if (prev->rev) { prev->rev->rev = NULL; prev->rev = NULL; }

      next->prev = prev;
      prev->next = next;

      prev = next = this;
    }



    template<unsigned ndim>
    void Edge<ndim>::insertBefore(Edge<ndim> *other) {
      if (prev != this) remove();
      prev = other->prev;
      next = other;
      next->prev = this;
      prev->next = this;

      if (prev->rev) { prev->rev->rev = NULL;  prev->rev = NULL; }
    }



    template<unsigned ndim> typename Edge<ndim>::iter_t Edge<ndim>::iter() { return iter_t(this); }
    template<unsigned ndim> typename Edge<ndim>::const_iter_t Edge<ndim>::iter() const { return const_iter_t(this); }

    template<unsigned ndim> typename Edge<ndim>::fwd_iter_t Edge<ndim>::begin() { return fwd_iter_t(this); }
    template<unsigned ndim> typename Edge<ndim>::const_fwd_iter_t Edge<ndim>::begin() const { return const_fwd_iter_t(this); }

    template<unsigned ndim> typename Edge<ndim>::fwd_iter_t Edge<ndim>::end() { return fwd_iter_t(NULL); }
    template<unsigned ndim> typename Edge<ndim>::const_fwd_iter_t Edge<ndim>::end() const { return const_fwd_iter_t(NULL); }

    template<unsigned ndim> typename Edge<ndim>::vert_iter_t Edge<ndim>::viter() { return vert_iter_t(this); }
    template<unsigned ndim> typename Edge<ndim>::const_vert_iter_t Edge<ndim>::viter() const { return const_vert_iter_t(this); }

    template<unsigned ndim> typename Edge<ndim>::fwd_vert_iter_t Edge<ndim>::vbegin() { return fwd_vert_iter_t(this); }
    template<unsigned ndim> typename Edge<ndim>::const_fwd_vert_iter_t Edge<ndim>::vbegin() const { return const_fwd_vert_iter_t(this); }

    template<unsigned ndim> typename Edge<ndim>::fwd_vert_iter_t Edge<ndim>::vend() { return fwd_vert_iter_t(NULL); }
    template<unsigned ndim> typename Edge<ndim>::const_fwd_vert_iter_t Edge<ndim>::vend() const { return const_fwd_vert_iter_t(NULL); }



    template<unsigned ndim>
    void Edge<ndim>::insertAfter(Edge<ndim> *other) {
      if (prev != this) remove();
      next = other->next;
      prev = other;
      next->prev = this;
      prev->next = this;

      if (prev->rev) { prev->rev->rev = NULL;  prev->rev = NULL; }
    }



    template<unsigned ndim>
    size_t Edge<ndim>::loopSize() const {
      const Edge *e = this;
      size_t n = 0;
      do { e = e->next; ++n; } while (e != this);
      return n;
    }



    Edge<ndim>::Edge(vertex_t *_vert, face_t *_face) :
        vert(_vert), face(_face), prev(this), next(this), rev(NULL) {
      CARVE_ASSERT(face != NULL);
    }



    template<unsigned ndim>
    Edge<ndim>::~Edge() {
      remove();
    }

    template<unsigned ndim>
    typename Face<ndim>::aabb_t Face<ndim>::getAABB() const {
      aabb_t aabb;
      aabb.fit(edge->begin(), edge->end(), vector_mapping());
      return aabb;
    }

    template<unsigned ndim>
    bool Face<ndim>::recalc() {
      if (!carve::geom3d::fitPlane(edge->begin(), edge->end(), vector_mapping(), plane)) {
        return false;
      }

      int da = carve::geom::largestAxis(plane.N);
      double A = carve::geom2d::signedArea(edge->begin(), edge->end(), projection_mapping(getProjector(false, da)));

      if ((A < 0.0) ^ (plane.N.v[da] < 0.0)) {
        plane.negate();
      }

      project = getProjector(plane.N.v[da] > 0, da);
      unproject = getUnprojector(plane.N.v[da] > 0, da);

      return true;
    }

    template<unsigned ndim>
    void Face<ndim>::clearEdges() {
      if (!edge) return;

      edge_t *curr = edge;
      while (1) {
        edge_t *next = curr->next;
        delete curr;
        if (curr == next) break;
        curr = next;
      }
      edge = NULL;

      n_edges = 0;
    }

    template<unsigned ndim>
    template<typename iter_t>
    void Face<ndim>::loopFwd(iter_t begin, iter_t end) {
      clearEdges();
      if (begin == end) return;
      edge = new edge_t(*begin, this); ++n_edges; ++begin;
      while (begin != end) {
        edge_t *e = new edge_t(*begin, this);
        e->insertAfter(edge->prev);
        ++n_edges;
        ++begin;
      }
    }

    template<unsigned ndim>
    template<typename iter_t>
    void Face<ndim>::loopRev(iter_t begin, iter_t end) {
      clearEdges();
      if (begin == end) return;
      edge = new edge_t(*begin, this); ++n_edges; ++begin;
      while (begin != end) {
        edge_t *e = new edge_t(*begin, this);
        e->insertBefore(edge->next);
        ++n_edges;
        ++begin;
      }
    }

    template<unsigned ndim>
    template<typename iter_t>
    void Face<ndim>::init(iter_t begin, iter_t end) {
      loopFwd(begin, end);
    }

    template<unsigned ndim>
    void Face<ndim>::init(vertex_t *a, vertex_t *b, vertex_t *c) {
      clearEdges();
      edge_t *ea = new edge_t(a, this);
      edge_t *eb = new edge_t(b, this);
      edge_t *ec = new edge_t(c, this);
      eb->insertAfter(ea);
      ec->insertAfter(eb);
      edge = ea;
      n_edges = 3;
    }

    template<unsigned ndim>
    void Face<ndim>::init(vertex_t *a, vertex_t *b, vertex_t *c, vertex_t *d) {
      clearEdges();
      edge_t *ea = new edge_t(a, this);
      edge_t *eb = new edge_t(b, this);
      edge_t *ec = new edge_t(c, this);
      edge_t *ed = new edge_t(d, this);
      eb->insertAfter(ea);
      ec->insertAfter(eb);
      ed->insertAfter(ec);
      edge = ea;
      n_edges = 4;
    }

    template<unsigned ndim>
    void Face<ndim>::getVertices(std::vector<const vertex_t *> &verts) const {
      verts.clear();
      verts.reserve(n_edges);
      const edge_t *e = edge;
      do { verts.push_back(e->vert); e = e->next; } while (e != edge);
    }

    template<unsigned ndim>
    void Face<ndim>::getProjectedVertices(std::vector<carve::geom::vector<2> > &verts) const {
      verts.clear();
      verts.reserve(n_edges);
      const edge_t *e = edge;
      do { verts.push_back(project(e->vert->v)); e = e->next; } while (e != edge);
    }



    template<unsigned ndim>
    Mesh<ndim>::Mesh(std::vector<face_t *> &_faces) : faces(), open_edges(), closed_edges(), meshset(NULL) {
      faces.swap(_faces);

      edge_t *emin = faces[0]->edge;

      for (size_t i = 0; i < faces.size(); ++i) {
        face_t *face = faces[i];
        edge_t *e = face->edge;
        do {
          if (e->rev == NULL) {
            open_edges.push_back(e);
          } else if (e < e->rev) {
            closed_edges.push_back(e);
          }
          if (e->v1()->v < emin->v1()->v) emin = e;
          e = e->next;
        } while (e != face->edge);
      }

      if (open_edges.size()) {
        is_negative = false;
      } else {
        std::vector<face_t *> min_faces;
        edge_t *e = emin;
        do {
          min_faces.push_back(e->face);
          CARVE_ASSERT(e->rev != NULL);
          e = e->rev->next;
          CARVE_ASSERT(e->v1() == emin->v1());
        } while (e != emin);

        double max_abs_x = 0.0;
        for (size_t f = 0; f < min_faces.size(); ++f) {
          if (fabs(min_faces[f]->plane.N.x) > fabs(max_abs_x)) max_abs_x = min_faces[f]->plane.N.x;
        }

        is_negative =  max_abs_x > 0.0;
      }
    }

    template<unsigned ndim>
    Mesh<ndim>::~Mesh() {
      for (size_t i = 0; i < faces.size(); ++i) {
        delete faces[i];
      }
    }



  }
}
