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

#if defined(DEBUG)
#define DEBUG_CONTAINS_VERTEX
#endif

#include <carve/geom.hpp>
#include <carve/poly.hpp>

#include <carve/octree_impl.hpp>

#include <carve/timing.hpp>

#include <algorithm>

namespace {



  struct VPtrSort {
    bool operator()(carve::poly::Vertex const *a, carve::poly::Vertex const *b) const {
      return a->v < b->v;
    }
  };



  struct FV {
    carve::poly::Face *face;
    size_t vertex;
    FV(carve::poly::Face *f, size_t v) : face(f), vertex(v) { }
  };



  struct EdgeFaces {
    std::list<FV> fwd, rev;
    carve::poly::Edge *edge;
  };



  // Interstingly, for one set of inserts and a number of complete
  // traversals, a map seems to be faster than an
  // unordered_map. This may apply in other places.

  struct EdgeFaceMap :
    public std::map<std::pair<const carve::poly::Vertex *, const carve::poly::Vertex *>,
                    EdgeFaces> {
    void record(const carve::poly::Vertex *v1,
                const carve::poly::Vertex *v2,
                carve::poly::Face *f,
                size_t i) {
      if (v1 < v2) {
        (*this)[std::make_pair(v1, v2)].fwd.push_back(FV(f, i));
      } else {
        (*this)[std::make_pair(v2, v1)].rev.push_back(FV(f, i));
      }
    }
  };



  struct FaceOrder {
    double ang;
    const FV *fv;
    bool fwd;

    FaceOrder(double _ang, const FV *_fv, bool _fwd) : ang(_ang), fv(_fv), fwd(_fwd) { }
  };



  bool operator<(const FaceOrder &a, const FaceOrder &b) { 
    return a.ang < b.ang || (a.ang == b.ang && a.fwd && !b.fwd); 
  }



  static inline std::ostream &operator<<(std::ostream &o, const FaceOrder &a) {
    o << (a.fwd ? "+" : "-") << " " << a.ang << " " << a.fv;
    return o;
  }



  bool makeFacePairs(EdgeFaces &ef,
                     carve::poly::Edge *e,
                     std::vector<const carve::poly::Face *> &edge_face_pairs) {
    edge_face_pairs.clear();

    static carve::TimingName FUNC_NAME("static Polyhedron makeFacePairs()");
    carve::TimingBlock block(FUNC_NAME);
  
    carve::geom3d::Vector evec = (e->v2->v - e->v1->v).normalized();
    std::vector<FaceOrder> sorted_faces;

    if (ef.fwd.size() == 0) {
      for (std::list<FV>::const_iterator
             f_i = ef.rev.begin(), f_e = ef.rev.end(); f_i != f_e; ++f_i) {
        const FV &fv2 = ((*f_i));

        edge_face_pairs.push_back(NULL);
        edge_face_pairs.push_back(fv2.face);
      }
      return true;
    } else if (ef.rev.size() == 0) {
      for (std::list<FV>::const_iterator
             f_i = ef.fwd.begin(), f_e = ef.fwd.end(); f_i != f_e; ++f_i) {
        const FV &fv1 = ((*f_i));

        edge_face_pairs.push_back(fv1.face);
        edge_face_pairs.push_back(NULL);
      }
      return true;
    }

    carve::geom3d::Vector base;

    base = ef.fwd.front().face->plane_eqn.N;

    for (std::list<FV>::const_iterator
           f_i = ef.fwd.begin(), f_e = ef.fwd.end(); f_i != f_e; ++f_i) {
      double ang = carve::geom3d::antiClockwiseAngle((*f_i).face->plane_eqn.N, base, evec);
      if (ang == 0.0 && f_i != ef.fwd.begin()) ang = M_TWOPI + carve::EPSILON;
      sorted_faces.push_back(FaceOrder(ang, &(*f_i), true));
    }
    for (std::list<FV>::const_iterator
           f_i = ef.rev.begin(), f_e = ef.rev.end(); f_i != f_e; ++f_i) {
      double ang = carve::geom3d::antiClockwiseAngle(-(*f_i).face->plane_eqn.N, base, evec);
      if (ang == 0.0) ang = M_TWOPI + carve::EPSILON;
      sorted_faces.push_back(FaceOrder(ang, &(*f_i), false));
    }
    std::sort(sorted_faces.begin(), sorted_faces.end());

    for (unsigned i = 0; i < sorted_faces.size();) {
      if (!sorted_faces[i].fwd) {
        const FV &fv2 = (*(sorted_faces[i++].fv));

        edge_face_pairs.push_back(NULL);
        edge_face_pairs.push_back(fv2.face);
      } else if (i == sorted_faces.size() - 1 || sorted_faces[i + 1].fwd) {
        const FV &fv1 = (*(sorted_faces[i++].fv));

        edge_face_pairs.push_back(fv1.face);
        edge_face_pairs.push_back(NULL);
      } else {
        const FV &fv1 = (*(sorted_faces[i++].fv));
        const FV &fv2 = (*(sorted_faces[i++].fv));

        edge_face_pairs.push_back(fv1.face);
        edge_face_pairs.push_back(fv2.face);
      }
    }

    return true;
  }



  bool emb_test(carve::poly::Polyhedron *poly,
                 std::map<int, std::set<int> > &embedding,
                 carve::geom3d::Vector v,
                 int m_id) {

    std::map<int, carve::PointClass> result;
#if defined(DEBUG)
    std::cerr << "test " << v << " (m_id:" << m_id << ")" << std::endl;
#endif
    poly->testVertexAgainstClosedManifolds(v, result, true);
    std::set<int> inside;
    for (std::map<int, carve::PointClass>::iterator j = result.begin();
         j != result.end();
         ++j) {
      if ((*j).first == m_id) continue;
      if ((*j).second == carve::POINT_IN) inside.insert((*j).first);
      else if ((*j).second == carve::POINT_ON) {
#if defined(DEBUG)
        std::cerr << " FAIL" << std::endl;
#endif
        return false;
      }
    }
#if defined(DEBUG)
    std::cerr << " OK (inside.size()==" << inside.size() << ")" << std::endl;
#endif
    embedding[m_id] = inside;
    return true;
  }



  struct order_faces {
    bool operator()(const carve::poly::Face * const &a, const carve::poly::Face * const &b) const {
      return std::lexicographical_compare(a->vertices.begin(), a->vertices.end(),
                                          b->vertices.begin(), b->vertices.end());
    }
  };



}



namespace carve {
  namespace poly {



    bool Polyhedron::initSpatialIndex() {
      static carve::TimingName FUNC_NAME("Polyhedron::initSpatialIndex()");
      carve::TimingBlock block(FUNC_NAME);

      octree.setBounds(aabb);
      octree.addFaces(faces);
      octree.addEdges(edges);
      octree.splitTree();

      return true;
    }



    void Polyhedron::invertAll() {
      for (size_t i = 0; i < faces.size(); ++i) {
        faces[i].invert();
      }

      for (size_t i = 0; i < edges.size(); ++i) {
        std::vector<const Face *> &f = connectivity.edge_to_face[i];
        for (size_t j = 0; j < (f.size() & ~1U); j += 2) {
          std::swap(f[j], f[j+1]);
        }
      }

      for (size_t i = 0; i < manifold_is_negative.size(); ++i) {
        manifold_is_negative[i] = !manifold_is_negative[i];
      }
    }



    void Polyhedron::invert(const std::vector<bool> &selected_manifolds) {
      bool altered = false;
      for (size_t i = 0; i < faces.size(); ++i) {
        if (faces[i].manifold_id >= 0 &&
            (unsigned)faces[i].manifold_id < selected_manifolds.size() &&
            selected_manifolds[faces[i].manifold_id]) {
          altered = true;
          faces[i].invert();
        }
      }

      if (altered) {
        for (size_t i = 0; i < edges.size(); ++i) {
          std::vector<const Face *> &f = connectivity.edge_to_face[i];
          for (size_t j = 0; j < (f.size() & ~1U); j += 2) {
            int m_id = -1;
            if (f[j]) m_id = f[j]->manifold_id;
            if (f[j+1]) m_id = f[j+1]->manifold_id;
            if (m_id >= 0 && (unsigned)m_id < selected_manifolds.size() && selected_manifolds[m_id]) {
              std::swap(f[j], f[j+1]);
            }
          }
        }

        for (size_t i = 0; i < std::min(selected_manifolds.size(), manifold_is_negative.size()); ++i) {
          manifold_is_negative[i] = !manifold_is_negative[i];
        }
      }
    }



    Polyhedron *Polyhedron::makeCopy(const std::vector<bool> &selected_manifolds) const {
      std::vector<Face> newFaces;
      size_t n_faces = 0;

      for (size_t i = 0; i < faces.size(); ++i) {
        if (faces[i].manifold_id >= 0 &&
            (unsigned)faces[i].manifold_id < selected_manifolds.size() &&
            selected_manifolds[faces[i].manifold_id]) n_faces++;
      }
      newFaces.reserve(n_faces);

      for (size_t i = 0; i < faces.size(); ++i) {
        const Face &src = faces[i];
        if (src.manifold_id < 0 ||
            (unsigned)src.manifold_id >= selected_manifolds.size() ||
            !selected_manifolds[src.manifold_id]) continue;
        newFaces.push_back(src);
      }

      return new Polyhedron(newFaces);
    }



    Polyhedron *Polyhedron::makeCopy() const {
      std::vector<Face> newFaces;

      newFaces.reserve(faces.size());

      for (size_t i = 0; i < faces.size(); ++i) {
        const Face &src = faces[i];
        newFaces.push_back(src);
      }

      return new Polyhedron(newFaces);
    }



    bool Polyhedron::buildEdges() {
      static carve::TimingName COUNT_VERTEX_FACES("Polyhedron::buildEdges() - count vertex faces");
      static carve::TimingName MAKE_HASHMAP("Polyhedron::buildEdges() - construct hashmap");
      static carve::TimingName MAKE_PAIRS("Polyhedron::buildEdges() - construct pairs");

      static carve::TimingName FUNC_NAME("Polyhedron::buildEdges()");
      carve::TimingBlock block(FUNC_NAME);

      EdgeFaceMap ef_map;
      bool is_ok = true;

      edges.clear();

      std::vector<size_t> vertex_face_count;

      carve::Timing::start(COUNT_VERTEX_FACES);

      vertex_face_count.resize(poly_vertices.size());

      // work out how many faces/edges each vertex is connected to, in
      // order to save on array reallocs.
      for (unsigned i = 0; i < faces.size(); ++i) {
        std::vector<const Vertex *> &v = faces[i].vertices;
        for (unsigned j = 0; j < v.size(); j++) {
          vertex_face_count[vertexToIndex_fast(v[j])]++;
        }
      }

      carve::Timing::stop();

      carve::Timing::start(MAKE_HASHMAP);

      // make a mapping from pairs of vertices denoting edges to pairs
      // of <face,vertex index> that incorporate this edge in the
      // forward and reverse directions.
      for (unsigned i = 0; i < faces.size(); ++i) {
        Face &f = faces[i];
        std::vector<const Vertex *> &v = f.vertices;

        for (unsigned j = 0; j < v.size() - 1; j++) {
          ef_map.record(v[j], v[j+1], &f, j);
        }
        ef_map.record(v.back(), v.front(), &f, v.size() - 1);

        f.edges.clear();
        f.edges.resize(v.size(), NULL);
        f.manifold_id = -1;
      }

      // now we know how many edges this polyhedron has.
      edges.reserve(ef_map.size());

      carve::Timing::stop();

      carve::Timing::start(MAKE_PAIRS);

      // make an edge object for each entry in ef_map.
      for (EdgeFaceMap::iterator i = ef_map.begin(), e = ef_map.end();
           i != e;
           ++i) {
        Vertex *v1 = const_cast<Vertex *>((*i).first.first);
        Vertex *v2 = const_cast<Vertex *>((*i).first.second);

        edges.push_back(Edge(v1, v2, this));
        (*i).second.edge = &edges.back();
      }

      // allocate space for connectivity info.
      connectivity.vertex_to_edge.resize(poly_vertices.size());
      connectivity.vertex_to_face.resize(poly_vertices.size());
      connectivity.edge_to_face.resize(edges.size());

      for (size_t i = 0; i < poly_vertices.size(); ++i) {
        connectivity.vertex_to_edge[i].reserve(vertex_face_count[i]);
        connectivity.vertex_to_face[i].reserve(vertex_face_count[i]);
      }

      // record connectivity from vertex to edges.
      for (size_t i = 0; i < edges.size(); ++i) {
        size_t v1i = vertexToIndex_fast(edges[i].v1);
        size_t v2i = vertexToIndex_fast(edges[i].v2);

        connectivity.vertex_to_edge[v1i].push_back(&edges[i]);
        connectivity.vertex_to_edge[v2i].push_back(&edges[i]);
      }

      // record connectivity from vertex to faces.
      for (size_t i = 0; i < faces.size(); ++i) {
        Face &f = faces[i];
        std::vector<const Vertex *> &v = f.vertices;

        for (unsigned j = 0; j < v.size(); j++) {
          size_t vi = vertexToIndex_fast(v[j]);
          connectivity.vertex_to_face[vi].push_back(&f);
        }
      }

      // set up face->edge lists.
      for (EdgeFaceMap::iterator i = ef_map.begin(); i != ef_map.end(); ++i) {
        Edge *edge = (*i).second.edge;
        std::list<FV> &fwd_faces = ((*i).second.fwd);
        std::list<FV> &rev_faces = ((*i).second.rev);
        for (std::list<FV>::iterator j = fwd_faces.begin(); j != fwd_faces.end(); ++j) {
          ASSERT((*j).vertex < (*j).face->edges.size());

          (*j).face->edges[(*j).vertex] = edge;
        }

        for (std::list<FV>::iterator j = rev_faces.begin(); j != rev_faces.end(); ++j) {
          ASSERT((*j).vertex < (*j).face->edges.size());

          (*j).face->edges[(*j).vertex] = edge;
        }
      }

      // pair up incident faces for each edge.
      std::vector<const Face *> edge_face_pairs;
      for (EdgeFaceMap::iterator i = ef_map.begin(); i != ef_map.end(); ++i) {
        std::list<FV> &fwd_faces = ((*i).second.fwd);
        std::list<FV> &rev_faces = ((*i).second.rev);

        Edge *edge = (*i).second.edge;
        size_t edge_index = edgeToIndex_fast(edge);

        if (fwd_faces.size() == 1 && rev_faces.size() == 1) {
          Face *f1 = fwd_faces.front().face;
          Face *f2 = rev_faces.front().face;
          edge_face_pairs.clear();
          edge_face_pairs.push_back(f1);
          edge_face_pairs.push_back(f2);

          connectivity.edge_to_face[edge_index] = edge_face_pairs;
        } else {
          if (makeFacePairs((*i).second, edge, edge_face_pairs)) {
            connectivity.edge_to_face[edge_index] = edge_face_pairs;
          } else {
            is_ok = false;
          }
        }
      }

      carve::Timing::stop();
  
      return is_ok;
    }



    bool Polyhedron::sortVertices() {
      static carve::TimingName FUNC_NAME("Polyhedron::sortVertices()");
      carve::TimingBlock block(FUNC_NAME);

      std::vector<Vertex *> vptr;
      std::map<const Vertex *, const Vertex *> vmap;
      std::vector<Vertex> vout;
      const std::vector<const Vertex *>::size_type l = poly_vertices.size();

      vptr.reserve(poly_vertices.size());
      vout.reserve(poly_vertices.size());

      for (size_t i = 0; i != l; ++i) {
        vptr.push_back(&poly_vertices[i]);
      }
      std::sort(vptr.begin(), vptr.end(), VPtrSort());

      for (size_t i = 0; i != l; ++i) {
        vout.push_back(*vptr[i]);
        vout.back().owner = this;
        vmap[vptr[i]] = &vout.back();
      }

      vout.swap(poly_vertices);

      for (size_t i = 0; i < faces.size(); ++i) {
        Face &f = faces[i];
        for (size_t j = 0; j < f.vertices.size(); ++j) {
          f.vertices[j] = vmap[f.vertices[j]];
        }
        f.owner = this;
      }
      for (size_t i = 0; i < edges.size(); ++i) {
        edges[i].v1 = vmap[edges[i].v1];
        edges[i].v2 = vmap[edges[i].v2];
      }

      return true;
    }



    bool Polyhedron::calcManifoldEmbedding() {
      static carve::TimingName FUNC_NAME("Polyhedron::calcManifoldEmbedding()");
      static carve::TimingName CME_V("Polyhedron::calcManifoldEmbedding() (vertices)");
      static carve::TimingName CME_E("Polyhedron::calcManifoldEmbedding() (edges)");
      static carve::TimingName CME_F("Polyhedron::calcManifoldEmbedding() (faces)");

      carve::TimingBlock block(FUNC_NAME);

      const unsigned MCOUNT = manifoldCount();
      if (MCOUNT < 2) return true;

      std::set<int> vertex_manifolds;
      std::map<int, std::set<int> > embedding;

      carve::Timing::start(CME_V);
      for (size_t i = 0; i < poly_vertices.size(); ++i) {
        vertex_manifolds.clear();
        if (vertexManifolds(&poly_vertices[i], set_inserter(vertex_manifolds)) != 1) continue;
        int m_id = *vertex_manifolds.begin();
        if (embedding.find(m_id) == embedding.end()) {
          if (emb_test(this, embedding, poly_vertices[i].v, m_id) && embedding.size() == MCOUNT) {
            carve::Timing::stop();
            goto done;
          }
        }
      }
      carve::Timing::stop();

      carve::Timing::start(CME_E);
      for (size_t i = 0; i < edges.size(); ++i) {
        if (connectivity.edge_to_face[i].size() == 2) {
          int m_id;
          const Face *f1 = connectivity.edge_to_face[i][0];
          const Face *f2 = connectivity.edge_to_face[i][1];
          if (f1) m_id = f1->manifold_id;
          if (f2) m_id = f2->manifold_id;
          if (embedding.find(m_id) == embedding.end()) {
            if (emb_test(this, embedding, (edges[i].v1->v + edges[i].v2->v) / 2, m_id) && embedding.size() == MCOUNT) {
              carve::Timing::stop();
              goto done;
            }
          }
        }
      }
      carve::Timing::stop();

      carve::Timing::start(CME_F);
      for (size_t i = 0; i < faces.size(); ++i) {
        int m_id = faces[i].manifold_id;
        if (embedding.find(m_id) == embedding.end()) {
          carve::geom2d::P2 pv;
          if (!carve::geom2d::pickContainedPoint(faces[i].vertices, p2_adapt_project(faces[i].project), pv)) continue;
          carve::geom3d::Vector v = carve::poly::face::unproject(faces[i], pv);
          if (emb_test(this, embedding, v, m_id) && embedding.size() == MCOUNT) {
            carve::Timing::stop();
            goto done;
          }
        }
      }
      carve::Timing::stop();

      std::cerr << "could not find test points!!!" << std::endl;
      return true;

      ASSERT(!!!"could not find test points");

    done:;
      for (std::map<int, std::set<int> >::iterator i = embedding.begin(); i != embedding.end(); ++i) {
#if defined(DEBUG)
        std::cerr << (*i).first << " : ";
        std::copy((*i).second.begin(), (*i).second.end(), std::ostream_iterator<int>(std::cerr, ","));
        std::cerr << std::endl;
#endif
        (*i).second.insert(-1);
      }
      std::set<int> parents, new_parents;
      parents.insert(-1);

      while (embedding.size()) {
        new_parents.clear();
        for (std::map<int, std::set<int> >::iterator i = embedding.begin(); i != embedding.end(); ++i) {
          if ((*i).second.size() == 1) {
            if (parents.find(*(*i).second.begin()) != parents.end()) {
              new_parents.insert((*i).first);
#if defined(DEBUG)
              std::cerr << "parent(" << (*i).first << "): " << *(*i).second.begin() << std::endl;
#endif
            } else {
#if defined(DEBUG)
              std::cerr << "no parent: " << (*i).first << " (looking for: " << *(*i).second.begin() << ")" << std::endl;
#endif
            }
          }
        }
        for (std::set<int>::const_iterator i = new_parents.begin(); i != new_parents.end(); ++i) {
          embedding.erase(*i);
        }
        for (std::map<int, std::set<int> >::iterator i = embedding.begin(); i != embedding.end(); ++i) {
          size_t n = 0;
          for (std::set<int>::const_iterator j = parents.begin(); j != parents.end(); ++j) {
            n += (*i).second.erase((*j));
          }
          ASSERT(n != 0);
        }
        parents.swap(new_parents);
      }

      return true;
    }



    bool Polyhedron::markManifolds() {
      static carve::TimingName FUNC_NAME("Polyhedron::markManifolds()");
      carve::TimingBlock block(FUNC_NAME);

      std::vector<Face *> to_mark;
      size_t i = 0;
      int m_id = 0;
      int closed_manifold_count = 0;

      const Vertex *min_vertex = NULL;
      std::set<const Face *> min_faces;

      manifold_is_closed.clear();
      manifold_is_negative.clear();

      while (1) {

        while (i < faces.size() && faces[i].manifold_id != -1) ++i;
        if (i == faces.size()) break;

        to_mark.push_back(&faces[i]);
        min_vertex = faces[i].vertices[0];

        bool is_closed = true;

        while (to_mark.size()) {
          Face *f = to_mark.back();
          to_mark.pop_back();

          if (f->manifold_id == -1) {
            f->manifold_id = m_id;

            const Vertex *v = f->vertices[0];
            for (size_t j = 1; j < f->vertices.size(); ++j) {
              if (f->vertices[j]->v < v->v) {
                v = f->vertices[j];
              }
            }
            if (v->v < min_vertex->v) {
              min_vertex = v;
            }

            for (size_t j = 0; j < f->edges.size(); ++j) {
              Face *g = const_cast<Face *>(connectedFace(f, f->edges[j]));

              if (g) {
                if (g->manifold_id == -1) to_mark.push_back(g);
              } else {
                is_closed = false;
              }
            }
          }
        }

        vertexToFaces(min_vertex, set_inserter(min_faces));

        double max_abs_x = 0.0;
        for (std::set<const Face *>::iterator i = min_faces.begin(); i != min_faces.end(); ++i) {
          if (fabs((*i)->plane_eqn.N.x) > fabs(max_abs_x)) max_abs_x = (*i)->plane_eqn.N.x;
        }

        manifold_is_closed.push_back(is_closed);
        manifold_is_negative.push_back(is_closed && max_abs_x > 0.0);
#if defined(DEBUG)
        std::cerr << "{manifold: " << m_id << (manifold_is_negative.back() ? " is" : " is not") << " negative}" << std::endl;
#endif
        if (is_closed) closed_manifold_count++;
        ++m_id;
      }

#if defined(DEBUG)
      std::cerr << "polyhedron " << this << " has " << m_id << " manifolds (" << closed_manifold_count << " closed)" << std::endl;
#endif

      return true;
    }



    bool Polyhedron::init() {
      static carve::TimingName FUNC_NAME("Polyhedron::init()");
      carve::TimingBlock block(FUNC_NAME);
  
      aabb.fit(poly_vertices.begin(), poly_vertices.end(), vec_adapt_vertex_ref());

      for (size_t i = 0; i < poly_vertices.size(); ++i) poly_vertices[i].owner = this;
      
      return sortVertices() && buildEdges() && initSpatialIndex() && markManifolds(); //  && calcManifoldEmbedding();
    }



    void Polyhedron::faceRecalc() {
      for (size_t i = 0; i < faces.size(); ++i) {
        if (!faces[i].recalc()) {
          std::ostringstream out;
          out << "face " << i << " recalc failed";
          throw carve::exception(out.str());
        }
      }
    }



    Polyhedron::Polyhedron(const std::vector<carve::geom3d::Vector> &vertices, int n_faces, const std::vector<int> &face_indices) {

      // okay, our polyhedron has a vector of vertices, which we want to copy, 
      // and we need to generate a set of Face*'s from it's face index list.
      poly_vertices.clear();
      poly_vertices.resize(vertices.size());
      for (size_t i = 0; i < vertices.size(); ++i) {
        poly_vertices[i].v = vertices[i];
      }

      faces.reserve(n_faces);
  
      std::vector<int>::const_iterator iter = face_indices.begin();
      std::vector<const Vertex *> v;
      for (int i = 0; i < n_faces; ++i) {
        int vertexCount = *iter++;
    
        v.clear();
    
        while (vertexCount--) {
          ASSERT(*iter >= 0);
          ASSERT((unsigned)*iter < poly_vertices.size());
          v.push_back(&poly_vertices[*iter++]);
        }
        faces.push_back(Face(v));
      }
      if (!init()) {
        // std::cerr << "polyhedron creation failed" << std::endl;
        throw carve::exception("polyhedron creation failed");
      }
    }



    Polyhedron::Polyhedron(std::vector<Face> &_faces, std::vector<Vertex> &_vertices, bool _recalc) {
      faces.swap(_faces);
      poly_vertices.swap(_vertices);

      if (_recalc) faceRecalc();

      if (!init()) {
        throw carve::exception("polyhedron creation failed");
      }
    }



    Polyhedron::Polyhedron(std::vector<Face> &_faces, bool _recalc) {
      faces.swap(_faces);
      commonFaceInit(_recalc);
    }



    Polyhedron::Polyhedron(std::list<Face> &_faces, bool _recalc) {
      faces.reserve(_faces.size());
      std::copy(_faces.begin(), _faces.end(), std::back_inserter(faces));
      commonFaceInit(_recalc);
    }



    void Polyhedron::collectFaceVertices(std::vector<Face> &faces,
                                         std::vector<Vertex> &vertices,
                                         carve::csg::VVMap &vmap) {
      vertices.clear();
      vmap.clear();

      for (size_t i = 0, il = faces.size(); i != il; ++i) {
        Face &f = faces[i];

        for (size_t j = 0, jl = f.vertices.size(); j != jl; ++j) {
          vmap[f.vertices[j]] = NULL;
        }
      }

      vertices.reserve(vmap.size());

      for (carve::csg::VVMap::iterator i = vmap.begin(),
             e = vmap.end();
           i != e;
           ++i) {
        vertices.push_back(*(*i).first);
        (*i).second = &vertices.back();
      }

      for (FacePtrVector::size_type i = 0, il = faces.size(); i != il; ++i) {
        Face &f = faces[i];

        for (size_t j = 0, jl = f.vertices.size(); j != jl; ++j) {
          f.vertices[j] = vmap[f.vertices[j]];
        }
      }
    }



    void Polyhedron::collectFaceVertices(std::vector<Face> &faces,
                                         std::vector<Vertex> &vertices) {
      std::unordered_map<const Vertex *,
        const Vertex *,
        hash_vertex_ptr> vmap;
      collectFaceVertices(faces, vertices, vmap);
    }



    void Polyhedron::commonFaceInit(bool _recalc) {
      collectFaceVertices(faces, poly_vertices);

      if (_recalc) faceRecalc();

      if (!init()) {
        throw carve::exception("polyhedron creation failed");
      }
    }



    Polyhedron::~Polyhedron() {
    }



    void Polyhedron::testVertexAgainstClosedManifolds(const carve::geom3d::Vector &v,
                                                      std::map<int, PointClass> &result,
                                                      bool ignore_orientation) const {

      for (FacePtrVector::size_type i = 0; i < faces.size(); i++) {
        if (!manifold_is_closed[faces[i].manifold_id]) continue; // skip open manifolds
        if (faces[i].containsPoint(v)) {
          result[faces[i].manifold_id] = POINT_ON;
        }
      }

      double ray_len = aabb.extent.length() * 2;

      std::vector<const Face *> possible_faces;

      std::vector<std::pair<const Face *, carve::geom3d::Vector> > manifold_intersections;

      while (1) {
        double a1 = random() / double(RAND_MAX) * M_TWOPI;
        double a2 = random() / double(RAND_MAX) * M_TWOPI;

        carve::geom3d::Vector ray_dir = carve::geom::VECTOR(sin(a1) * sin(a2), cos(a1) * sin(a2), cos(a2));

        carve::geom3d::Vector v2 = v + ray_dir * ray_len;

        bool failed = false;
        carve::geom3d::LineSegment line(v, v2);
        carve::geom3d::Vector intersection;

        possible_faces.clear();
        manifold_intersections.clear();
        octree.findFacesNear(line, possible_faces);

        for (unsigned i = 0; !failed && i < possible_faces.size(); i++) {
          if (!manifold_is_closed[possible_faces[i]->manifold_id]) continue; // skip open manifolds
          if (result.find(possible_faces[i]->manifold_id) != result.end()) continue; // already ON

          switch (possible_faces[i]->lineSegmentIntersection(line, intersection)) {
          case INTERSECT_FACE: {
            manifold_intersections.push_back(std::make_pair(possible_faces[i], intersection));
            break;
          }
          case INTERSECT_NONE: {
            break;
          }
          default: {
            failed = true;
            break;
          }
          }
        }

        if (!failed) break;
      }

      std::vector<int> crossings(manifold_is_closed.size(), 0);

      for (size_t i = 0; i < manifold_intersections.size(); ++i) {
        const Face *f = manifold_intersections[i].first;
        crossings[f->manifold_id]++;
      }

      for (size_t i = 0; i < crossings.size(); ++i) {
#if defined(DEBUG)
        std::cerr << "crossing: " << i << " = " << crossings[i] << " is_negative = " << manifold_is_negative[i] << std::endl;
#endif
        if (!manifold_is_closed[i]) continue;
        if (result.find(i) != result.end()) continue;
        PointClass pc = (crossings[i] & 1) ? POINT_IN : POINT_OUT;
        if (!ignore_orientation && manifold_is_negative[i]) pc = (PointClass)-pc;
        result[i] = pc;
      }
    }



    PointClass Polyhedron::containsVertex(const carve::geom3d::Vector &v,
                                          const Face **hit_face,
                                          bool even_odd,
                                          int manifold_id) const {
      if (hit_face) *hit_face = NULL;

#if defined(DEBUG_CONTAINS_VERTEX)
      std::cerr << "{containsVertex " << v << "}" << std::endl;
#endif

      if (!aabb.containsPoint(v)) {
#if defined(DEBUG_CONTAINS_VERTEX)
        std::cerr << "{final:OUT(aabb short circuit)}" << std::endl;
#endif
        // XXX: if the top level manifolds are negative, this should be POINT_IN.
        // for the moment, this only works for a single manifold.
        if (manifold_is_negative.size() == 1 && manifold_is_negative[0]) return POINT_IN;
        return POINT_OUT;
      }

      for (FacePtrVector::size_type i = 0; i < faces.size(); i++) {
        if (manifold_id != -1 && manifold_id != faces[i].manifold_id) continue;

        // XXX: Do allow the tested vertex to be ON an open
        // manifold. This was here originally because of the
        // possibility of an open manifold contained within a closed
        // manifold.

        // if (!manifold_is_closed[faces[i].manifold_id]) continue;

        if (faces[i].containsPoint(v)) {
#if defined(DEBUG_CONTAINS_VERTEX)
          std::cerr << "{final:ON(hits face " << faces[i] << ")}" << std::endl;
#endif
          if (hit_face) *hit_face = &faces[i];
          return POINT_ON;
        }
      }

      double ray_len = aabb.extent.length() * 2;

      std::vector<const Face *> possible_faces;

      std::vector<std::pair<const Face *, carve::geom3d::Vector> > manifold_intersections;

      while (1) {
        double a1 = random() / double(RAND_MAX) * M_TWOPI;
        double a2 = random() / double(RAND_MAX) * M_TWOPI;

        carve::geom3d::Vector ray_dir = carve::geom::VECTOR(sin(a1) * sin(a2), cos(a1) * sin(a2), cos(a2));

#if defined(DEBUG_CONTAINS_VERTEX)
        std::cerr << "{testing ray: " << ray_dir << "}" << std::endl;
#endif

        carve::geom3d::Vector v2 = v + ray_dir * ray_len;

        bool failed = false;
        carve::geom3d::LineSegment line(v, v2);
        carve::geom3d::Vector intersection;

        possible_faces.clear();
        manifold_intersections.clear();
        octree.findFacesNear(line, possible_faces);

        for (unsigned i = 0; !failed && i < possible_faces.size(); i++) {
          if (manifold_id != -1 && manifold_id != faces[i].manifold_id) continue;

          if (!manifold_is_closed[possible_faces[i]->manifold_id]) continue;

          switch (possible_faces[i]->lineSegmentIntersection(line, intersection)) {
          case INTERSECT_FACE: {

#if defined(DEBUG_CONTAINS_VERTEX)
            std::cerr << "{intersects face: " << possible_faces[i]
                      << " dp: " << dot(ray_dir, possible_faces[i]->plane_eqn.N) << "}" << std::endl;
#endif

            if (!even_odd && fabs(dot(ray_dir, possible_faces[i]->plane_eqn.N)) < EPSILON) {

#if defined(DEBUG_CONTAINS_VERTEX)
              std::cerr << "{failing(small dot product)}" << std::endl;
#endif

              failed = true;
              break;
            }
            manifold_intersections.push_back(std::make_pair(possible_faces[i], intersection));
            break;
          }
          case INTERSECT_NONE: {
            break;
          }
          default: {

#if defined(DEBUG_CONTAINS_VERTEX)
            std::cerr << "{failing(degenerate intersection)}" << std::endl;
#endif
            failed = true;
            break;
          }
          }
        }

        if (!failed) {
          if (even_odd) {
            return (manifold_intersections.size() & 1) ? POINT_IN : POINT_OUT;
          }

#if defined(DEBUG_CONTAINS_VERTEX)
          std::cerr << "{intersections ok [count:"
                    << manifold_intersections.size()
                    << "], sorting}"
                    << std::endl;
#endif

          carve::geom3d::sortInDirectionOfRay(ray_dir,
                                              manifold_intersections.begin(),
                                              manifold_intersections.end(),
                                              carve::geom3d::vec_adapt_pair_second());

          std::vector<int> crossings(manifold_is_closed.size(), 0);

          for (size_t i = 0; i < manifold_intersections.size(); ++i) {
            const Face *f = manifold_intersections[i].first;
            if (dot(ray_dir, f->plane_eqn.N) < 0.0) {
              crossings[f->manifold_id]++;
            } else {
              crossings[f->manifold_id]--;
            }
          }

#if defined(DEBUG_CONTAINS_VERTEX)
          for (size_t i = 0; i < crossings.size(); ++i) {
            std::cerr << "{manifold " << i << " crossing count: " << crossings[i] << "}" << std::endl;
          }
#endif

          for (size_t i = 0; i < manifold_intersections.size(); ++i) {
            const Face *f = manifold_intersections[i].first;

#if defined(DEBUG_CONTAINS_VERTEX)
            std::cerr << "{intersection at "
                      << manifold_intersections[i].second
                      << " id: "
                      << f->manifold_id
                      << " count: "
                      << crossings[f->manifold_id]
                      << "}"
                      << std::endl;
#endif

            if (crossings[f->manifold_id] < 0) {
              // inside this manifold.

#if defined(DEBUG_CONTAINS_VERTEX)
              std::cerr << "{final:IN}" << std::endl;
#endif

              return POINT_IN;
            } else if (crossings[f->manifold_id] > 0) {
              // outside this manifold, but it's an infinite manifold. (for instance, an inverted cube)

#if defined(DEBUG_CONTAINS_VERTEX)
              std::cerr << "{final:OUT}" << std::endl;
#endif

              return POINT_OUT;
            }
          }

#if defined(DEBUG_CONTAINS_VERTEX)
          std::cerr << "{final:OUT(default)}" << std::endl;
#endif

          return POINT_OUT;
        }
      }
    }



    void Polyhedron::findEdgesNear(const carve::geom3d::LineSegment &line,
                                   std::vector<const Edge*> &outEdges) const {
      outEdges.clear();
      octree.findEdgesNear(line, outEdges);
    }



    void Polyhedron::findEdgesNear(const carve::geom3d::Vector &v,
                                   std::vector<const Edge*> &outEdges) const {
      outEdges.clear();
      octree.findEdgesNear(v, outEdges);
    }



    void Polyhedron::findEdgesNear(const Face &face,
                                   std::vector<const Edge*> &edges) const {
      edges.clear();
      octree.findEdgesNear(face, edges);
    }



    void Polyhedron::findEdgesNear(const Edge &edge,
                                   std::vector<const Edge*> &outEdges) const {
      outEdges.clear();
      octree.findEdgesNear(edge, outEdges);
    }



    void Polyhedron::findFacesNear(const carve::geom3d::LineSegment &line,
                                   std::vector<const Face*> &outFaces) const {
      outFaces.clear();
      octree.findFacesNear(line, outFaces);
    }



    void Polyhedron::findFacesNear(const Edge &edge,
                                   std::vector<const Face*> &outFaces) const {
      outFaces.clear();
      octree.findFacesNear(edge, outFaces);
    }



    void Polyhedron::transform(const carve::math::Matrix &xform) {
      for (size_t i = 0; i < poly_vertices.size(); i++) {
        poly_vertices[i].v = xform * poly_vertices[i].v;
      }
      for (size_t i = 0; i < faces.size(); i++) {
        faces[i].recalc();
      }
      init();
    }



    void Polyhedron::print(std::ostream &o) const {
      o << "Polyhedron@" << this << " {" << std::endl;
      for (std::vector<Vertex>::const_iterator
             i = poly_vertices.begin(), e = poly_vertices.end(); i != e; ++i) {
        o << "  V@" << &(*i) << " " << (*i).v << std::endl;
      }
      for (std::vector<Edge>::const_iterator
             i = edges.begin(), e = edges.end(); i != e; ++i) {
        o << "  E@" << &(*i) << " {" << std::endl;
        o << "    V@" << (*i).v1 << " - " << "V@" << (*i).v2 << std::endl;
        const std::vector<const Face *> &faces = connectivity.edge_to_face[edgeToIndex_fast(&(*i))];
        for (size_t j = 0; j < (faces.size() & ~1U); j += 2) {
          o << "      fp: F@" << faces[j] << ", F@" << faces[j+1] << std::endl;
        }
        o << "  }" << std::endl;
      }
      for (std::vector<Face>::const_iterator
             i = faces.begin(), e = faces.end(); i != e; ++i) {
        o << "  F@" << &(*i) << " {" << std::endl;
        o << "    vertices {" << std::endl;
        for (std::vector<const Vertex *>::const_iterator
               j = (*i).vertices.begin(), je = (*i).vertices.end(); j != je; ++j) {
          o << "      V@" << (*j) << std::endl;
        }
        o << "    }" << std::endl;
        o << "    edges {" << std::endl;
        for (std::vector<const Edge *>::const_iterator
               j = (*i).edges.begin(), je = (*i).edges.end(); j != je; ++j) {
          o << "      E@" << (*j) << std::endl;
        }
        carve::geom::plane<3> p = (*i).plane_eqn;
        o << "    }" << std::endl;
        o << "    normal " << (*i).plane_eqn.N << std::endl;
        o << "    centroid " << (*i).centroid << std::endl;
        o << "    aabb " << (*i).aabb << std::endl;
        o << "    plane_eqn ";
        carve::geom::operator<< <3>(o, p);
        o << std::endl;
        o << "  }" << std::endl;
      }

      o << "}" << std::endl;
    }



    void Polyhedron::canonicalize() {
      sortVertices();
      for (size_t i = 0; i < faces.size(); i++) {
        Face &f = faces[i];
        size_t j = std::distance(f.vertices.begin(),
                                 std::min_element(f.vertices.begin(),
                                                  f.vertices.end()));
        if (j) {
          {
            std::vector<const Vertex *> temp;
            temp.reserve(f.vertices.size());
            std::copy(f.vertices.begin() + j, f.vertices.end(),       std::back_inserter(temp));
            std::copy(f.vertices.begin(),     f.vertices.begin() + j, std::back_inserter(temp));
            temp.swap(f.vertices);
          }
          {
            std::vector<const Edge *> temp;
            temp.reserve(f.vertices.size());
            std::copy(f.edges.begin() + j, f.edges.end(),       std::back_inserter(temp));
            std::copy(f.edges.begin(),     f.edges.begin() + j, std::back_inserter(temp));
            temp.swap(f.edges);
          }
        }
      }

      std::vector<Face *> face_ptrs;
      face_ptrs.reserve(faces.size());
      for (size_t i = 0; i < faces.size(); ++i) face_ptrs.push_back(&faces[i]);
      std::sort(face_ptrs.begin(), face_ptrs.end(), order_faces());
      std::vector<Face> sorted_faces;
      sorted_faces.reserve(faces.size());
      for (size_t i = 0; i < faces.size(); ++i) sorted_faces.push_back(*face_ptrs[i]);
      std::swap(faces, sorted_faces);
    }



  }
}

