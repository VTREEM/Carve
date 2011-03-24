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
    Edge<ndim> *Edge<ndim>::collapse() {
      next->v1()->v = (next->v1()->v + v1()->v) / 2;

      if (rev) {
        rev->next->prev = rev->prev;
        rev->prev->next = rev->next;
        if (rev->face) {
          rev->face->n_edges--;
          if (rev->face->edge == rev) {
            if (rev == rev->next) {
              rev->face->edge = NULL;
            } else {
              rev->face->edge = rev->next;
            }
          }
        }
        delete rev;
      }

      next->prev = prev;
      prev->next = next;
      if (face) {
        face->n_edges--;
        if (face->edge == this) {
          if (next == this) {
            face->edge = NULL;
            next = NULL;
          } else {
            face->edge = next;
          }
        }
      }

      Edge *n = next;
      delete this;
      return n;
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



    template<unsigned ndim>
    Edge<ndim> *Edge<ndim>::perimNext() const {
      if (rev) return NULL;
      Edge *e = next;
      while(e->rev) {
        e = e->rev->next;
      }
      return e;
    }



    template<unsigned ndim>
    Edge<ndim> *Edge<ndim>::perimPrev() const {
      if (rev) return NULL;
      Edge *e = prev;
      while(e->rev) {
        e = e->rev->prev;
      }
      return e;
    }



    template<unsigned ndim>
    Edge<ndim>::Edge(vertex_t *_vert, face_t *_face) :
        vert(_vert), face(_face), prev(this), next(this), rev(NULL) {
      CARVE_ASSERT(face != NULL);
    }



    template<unsigned ndim>
    Edge<ndim>::~Edge() {
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
      do {
        edge_t *next = curr->next;
        delete curr;
        curr = next;
      } while (curr != edge);

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
    typename Face<ndim>::vector_t Face<ndim>::centroid() const {
      vector_t v;
      edge_t *e = edge;
      do {
        v += e->vert->v;
        e = e->next;
      } while(e != edge);
      v /= n_edges;
      return v;
    }



    template<unsigned ndim>
    Face<ndim> *Face<ndim>::clone(const Face<ndim>::vertex_t *old_base,
                                  Face<ndim>::vertex_t *new_base,
                                  std::unordered_map<const Face<ndim>::edge_t *, Face<ndim>::edge_t *> &edge_map) const {
      Face *r = new Face(*this);

      edge_t *e = edge;
      edge_t *r_p = NULL;
      edge_t *r_e;
      do {
        r_e = new edge_t(e->vert - old_base + new_base, r);
        edge_map[e] = r_e;
        if (r_p) {
          r_p->next = r_e;
          r_e->prev = r_p;
        } else {
          r->edge = r_e;
        }
        r_p = r_e;
        e = e->next;
      } while (e != edge);
      r_e->next = r->edge;
      r->edge->prev = r_e;
      return r;
    }



    template<unsigned ndim>
    Mesh<ndim>::Mesh(std::vector<typename Mesh<ndim>::face_t *> &_faces,
                     std::vector<typename Mesh<ndim>::edge_t *> &_open_edges,
                     std::vector<typename Mesh<ndim>::edge_t *> &_closed_edges,
                     bool _is_negative) {
      std::swap(faces, _faces);
      std::swap(open_edges, _open_edges);
      std::swap(closed_edges, _closed_edges);
      is_negative = _is_negative;
      meshset = NULL;
      
      for (size_t i = 0; i < faces.size(); ++i) {
        faces[i]->mesh = this;
      }
    }



    namespace detail {
      template<typename iter_t>
      void FaceStitcher::initEdges(iter_t begin,
                                   iter_t end) {
        size_t c = 0;
        for (iter_t i = begin; i != end; ++i) {
          face_t *face = *i;
          CARVE_ASSERT(face->mesh == NULL); // for the moment, can only insert a face into a mesh once.

          face->id = c++;
          edge_t *e = face->edge;
          do {
            edges[vpair_t(e->v1(), e->v2())].push_back(e);
            e = e->next;
            if (e->rev) { e->rev->rev = NULL; e->rev = NULL; }
          } while (e != face->edge);
        }
        face_groups.init(c);
        is_open.clear();
        is_open.resize(c, false);
      }

      template<typename iter_t>
      void FaceStitcher::build(iter_t begin,
                               iter_t end,
                               std::vector<Mesh<3> *> &meshes) {
        // work out what set each face belongs to, and then construct
        // mesh instances for each set of faces.
        std::vector<size_t> index_set;
        std::vector<size_t> set_size;
        face_groups.get_index_to_set(index_set, set_size);

        std::vector<std::vector<face_t *> > mesh_faces;
        mesh_faces.resize(set_size.size());
        for (size_t i = 0; i < set_size.size(); ++i) {
          mesh_faces[i].reserve(set_size[i]);
        }
      
        for (iter_t i = begin; i != end; ++i) {
          face_t *face = *i;
          mesh_faces[index_set[face->id]].push_back(face);
        }

        meshes.clear();
        meshes.reserve(mesh_faces.size());
        for (size_t i = 0; i < mesh_faces.size(); ++i) {
          meshes.push_back(new Mesh<3>(mesh_faces[i]));
        }
      }

      template<typename iter_t>
      void FaceStitcher::create(iter_t begin,
                                iter_t end,
                                std::vector<Mesh<3> *> &meshes) {
        initEdges(begin, end);
        construct();
        build(begin, end, meshes);
      }
    }



    template<unsigned ndim>
    void Mesh<ndim>::cacheEdges() {
      closed_edges.clear();
      open_edges.clear();

      for (size_t i = 0; i < faces.size(); ++i) {
        face_t *face = faces[i];
        edge_t *e = face->edge;
        do {
          if (e->rev == NULL) {
            open_edges.push_back(e);
          } else if (e < e->rev) {
            closed_edges.push_back(e);
          }
          e = e->next;
        } while (e != face->edge);
      }
    }



    template<unsigned ndim>
    Mesh<ndim>::Mesh(std::vector<face_t *> &_faces) : faces(), open_edges(), closed_edges(), meshset(NULL) {
      faces.swap(_faces);

      cacheEdges();
      calcOrientation();
    }



    template<unsigned ndim>
    void Mesh<ndim>::calcOrientation() {
      if (open_edges.size() || !closed_edges.size()) {
        is_negative = false;
      } else {
        edge_t *emin = closed_edges[0];
        if (emin->rev->v1()->v < emin->v1()->v) emin = emin->rev;
        for (size_t i = 1; i < closed_edges.size(); ++i) {
          if (emin->v1()->v < closed_edges[i]->v1()->v) emin = closed_edges[i];
          if (emin->v1()->v < closed_edges[i]->rev->v1()->v) emin = closed_edges[i]->rev;
        }

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
    Mesh<ndim> *Mesh<ndim>::clone(const typename Mesh<ndim>::vertex_t *old_base,
                                  typename Mesh<ndim>::vertex_t *new_base) const {
      std::vector<face_t *> r_faces;
      std::vector<edge_t *> r_open_edges;
      std::vector<edge_t *> r_closed_edges;
      std::unordered_map<const edge_t *, edge_t *> edge_map;

      r_faces.reserve(faces.size());
      r_open_edges.reserve(r_open_edges.size());
      r_closed_edges.reserve(r_closed_edges.size());

      for (size_t i = 0; i < faces.size(); ++i) {
        r_faces.push_back(faces[i]->clone(old_base, new_base, edge_map));
      }
      for (size_t i = 0; i < closed_edges.size(); ++i) {
        r_closed_edges.push_back(edge_map[closed_edges[i]]);
        r_closed_edges.back()->rev = edge_map[closed_edges[i]->rev];
      }
      for (size_t i = 0; i < open_edges.size(); ++i) {
        r_open_edges.push_back(edge_map[open_edges[i]]);
      }

      return new Mesh(r_faces, r_open_edges, r_closed_edges, is_negative);
    }



    template<unsigned ndim>
    Mesh<ndim>::~Mesh() {
      for (size_t i = 0; i < faces.size(); ++i) {
        delete faces[i];
      }
    }



    template<unsigned ndim>
    template<typename iter_t>
    void Mesh<ndim>::create(iter_t begin, iter_t end, std::vector<Mesh<ndim> *> &meshes) {
      meshes.clear();
    }



    template<>
    template<typename iter_t>
    void Mesh<3>::create(iter_t begin, iter_t end, std::vector<Mesh<3> *> &meshes) {
      detail::FaceStitcher().create(begin, end, meshes);
    }



    template<unsigned ndim>
    MeshSet<ndim>::MeshSet(const std::vector<typename MeshSet<ndim>::vertex_t::vector_t> &points,
                           size_t n_faces,
                           const std::vector<int> &face_indices) {
      vertex_storage.reserve(points.size());
      std::vector<face_t *> faces;
      faces.reserve(n_faces);
      for (size_t i = 0; i < points.size(); ++i) {
        vertex_storage.push_back(vertex_t(points[i]));
      }

      std::vector<vertex_t *> v;
      size_t p = 0;
      for (size_t i = 0; i < n_faces; ++i) {
        const size_t N = face_indices[p++];
        v.clear();
        v.reserve(N);
        for (size_t j = 0; j < N; ++j) {
          v.push_back(&vertex_storage[face_indices[p++]]);
        }
        faces.push_back(new face_t(v.begin(), v.end()));
      }
      CARVE_ASSERT(p == face_indices.size());
      mesh_t::create(faces.begin(), faces.end(), meshes);
    }



    template<unsigned ndim>
    MeshSet<ndim>::MeshSet(std::vector<typename MeshSet<ndim>::vertex_t> &_vertex_storage,
                           std::vector<MeshSet<ndim>::mesh_t *> &_meshes) {
      vertex_storage.swap(_vertex_storage);
      meshes.swap(_meshes);

      for (size_t i = 0; i < meshes.size(); ++i) {
        meshes[i]->meshset = this;
      }
    }



    template<unsigned ndim>
    MeshSet<ndim>::MeshSet(std::vector<typename MeshSet<ndim>::mesh_t *> &_meshes) {
      meshes.swap(_meshes);
      std::unordered_map<vertex_t *, size_t> vert_idx;

      for (size_t m = 0; m < meshes.size(); ++m) {
        mesh_t *mesh = meshes[m];
        CARVE_ASSERT(mesh->meshset == NULL);
        mesh->meshset = this;
        for (size_t f = 0; f < mesh->faces.size(); ++f) {
          face_t *face = mesh->faces[f];
          edge_t *edge = face->edge;
          do {
            vert_idx[edge->vert] = 0;
            edge = edge->next;
          } while (edge != face->edge);
        }
      }

      vertex_storage.reserve(vert_idx.size());
      for (typename std::unordered_map<vertex_t *, size_t>::iterator i = vert_idx.begin(); i != vert_idx.end(); ++i) {
        (*i).second = vertex_storage.size();
        vertex_storage.push_back(*(*i).first);
      }

      for (size_t m = 0; m < meshes.size(); ++m) {
        mesh_t *mesh = meshes[m];
        for (size_t f = 0; f < mesh->faces.size(); ++f) {
          face_t *face = mesh->faces[f];
          edge_t *edge = face->edge;
          do {
            size_t i = vert_idx[edge->vert];
            edge->vert = &vertex_storage[i];
            edge = edge->next;
          } while (edge != face->edge);
        }
      }
    }



    template<unsigned ndim>
    MeshSet<ndim> *MeshSet<ndim>::clone() const {
      std::vector<vertex_t> r_vertex_storage = vertex_storage;
      std::vector<mesh_t *> r_meshes;
      for (size_t i = 0; i < meshes.size(); ++i) {
        r_meshes.push_back(meshes[i]->clone(&vertex_storage[0], &r_vertex_storage[0]));
      }
      return new MeshSet(r_vertex_storage, r_meshes);
    }



    template<unsigned ndim>
    MeshSet<ndim>::~MeshSet() {
      for (size_t i = 0; i < meshes.size(); ++i) {
        delete meshes[i];
      }
    }



    template<unsigned ndim>
    MeshSet<ndim>::FaceIter::FaceIter(const MeshSet<ndim> *_obj, size_t _mesh, size_t _face) : obj(_obj), mesh(_mesh), face(_face) {
    }



    template<unsigned ndim>
    void MeshSet<ndim>::FaceIter::fwd(size_t n) {
      if (mesh < obj->meshes.size()) {
        face += n;
        while (face >= obj->meshes[mesh]->faces.size()) {
          face -= obj->meshes[mesh++]->faces.size();
          if (mesh == obj->meshes.size()) { face = 0; break; }
        }
      }
    }



    template<unsigned ndim>
    void MeshSet<ndim>::FaceIter::rev(size_t n) {
      while (n > face) {
        n -= face;
        if (mesh == 0) { face = 0; return; }
        face = obj->meshes[--mesh]->faces.size() - 1;
      }
      face -= n;
    }



    template<unsigned ndim>
    void MeshSet<ndim>::FaceIter::adv(int n) {
      if (n > 0) {
        fwd((size_t)n);
      } else if (n < 0) {
        rev((size_t)-n);
      }
    }



    template<unsigned ndim>
    typename MeshSet<ndim>::FaceIter::super::difference_type MeshSet<ndim>::FaceIter::operator-(const FaceIter &other) const {
      CARVE_ASSERT(obj == other.obj);
      if (mesh == other.mesh) return face - other.face;

      size_t m = 0;
      for (size_t i = std::min(mesh, other.mesh) + 1; i < std::max(mesh, other.mesh); ++i) {
        m += obj->meshes[i]->faces.size();
      }

      if (mesh < other.mesh) {
        return -(obj->meshes[mesh]->faces.size() - face +
                 m +
                 other.face);
      } else {
        return +(obj->meshes[other.mesh]->faces.size() - other.face +
                 m +
                 face);
      }
    }
  }
}
