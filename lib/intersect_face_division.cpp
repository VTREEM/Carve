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

#include <carve/csg.hpp>
#include <carve/polyline.hpp>
#include <carve/debug_hooks.hpp>
#include <carve/timing.hpp>
#include <carve/triangulator.hpp>

#include <list>
#include <set>
#include <iostream>

#include <algorithm>
#include <assert.h>

#include "intersect_common.hpp"

typedef carve::poly::Polyhedron poly_t;

namespace {



  struct GraphEdge {
    GraphEdge *next;
    GraphEdge *prev;
    GraphEdge *loop_next;
    const poly_t::vertex_t *src;
    const poly_t::vertex_t *tgt;
    double ang;
    int visited;

    GraphEdge(const poly_t::vertex_t *_src, const poly_t::vertex_t *_tgt) :
      next(NULL), prev(NULL), loop_next(NULL),
      src(_src), tgt(_tgt),
      ang(0.0), visited(-1) {
    }
  };



  struct GraphEdges {
    GraphEdge *edges;
    carve::geom2d::P2 proj;

    GraphEdges() : edges(NULL), proj() {
    }
  };



  struct Graph {
    typedef std::unordered_map<const poly_t::vertex_t *, GraphEdges, carve::poly::hash_vertex_ptr> graph_t;

    graph_t graph;

    Graph() : graph() {
    }

    ~Graph() {
      int c = 0;

      GraphEdge *edge;
      for (graph_t::iterator i = graph.begin(), e =  graph.end(); i != e; ++i) {
        edge = (*i).second.edges;
        while (edge) {
          GraphEdge *temp = edge;
          ++c;
          edge = edge->next;
          delete temp;
        }
      }

      if (c) {
        std::cerr << "warning: "
                  << c
                  << " edges should have already been removed at graph destruction time"
                  << std::endl;
      }
    }

    const carve::geom2d::P2 &projection(const poly_t::vertex_t *v) const {
      graph_t::const_iterator i = graph.find(v);
      ASSERT(i != graph.end());
      return (*i).second.proj;
    }

    void computeProjection(const poly_t::face_t *face) {
      for (graph_t::iterator i = graph.begin(), e =  graph.end(); i != e; ++i) {
        (*i).second.proj = carve::poly::face::project(face, (*i).first->v);
      }
      for (graph_t::iterator i = graph.begin(), e =  graph.end(); i != e; ++i) {
        for (GraphEdge *e = (*i).second.edges; e; e = e->next) {
          e->ang = carve::math::ANG(carve::geom2d::atan2(projection(e->tgt) - projection(e->src)));
        }
      }
    }

    void print(const carve::csg::VertexIntersections *vi) const {
      for (graph_t::const_iterator i = graph.begin(), e =  graph.end(); i != e; ++i) {
        std::cerr << (*i).first << (*i).first->v << '(' << projection((*i).first).x << ',' << projection((*i).first).y << ") :";
        for (const GraphEdge *e = (*i).second.edges; e; e = e->next) {
          std::cerr << ' ' << e->tgt << e->tgt->v << '(' << projection(e->tgt).x << ',' << projection(e->tgt).y << ')';
        }
        std::cerr << std::endl;
        if (vi) {
          carve::csg::VertexIntersections::const_iterator j = vi->find((*i).first);
          if (j != vi->end()) {
            std::cerr << "   (int) ";
            for (carve::csg::IObjPairSet::const_iterator
                   k = (*j).second.begin(), ke = (*j).second.end(); k != ke; ++k) {
              if ((*k).first < (*k).second) {
                std::cerr << (*k).first << ".." << (*k).second << "; ";
              }
            }
            std::cerr << std::endl;
          }
        }
      }
    }

    void addEdge(const poly_t::vertex_t *v1, const poly_t::vertex_t *v2) {
      GraphEdges &edges = graph[v1];
      GraphEdge *edge = new GraphEdge(v1, v2);
      if (edges.edges) edges.edges->prev = edge;
      edge->next = edges.edges;
      edges.edges = edge;
    }

    void removeEdge(GraphEdge *edge) {
      if (edge->prev != NULL) {
        edge->prev->next = edge->next;
      } else {
        if (edge->next != NULL) {
          GraphEdges &edges = (graph[edge->src]);
          edges.edges = edge->next;
        } else {
          graph.erase(edge->src);
        }
      }
      if (edge->next != NULL) {
        edge->next->prev = edge->prev;
      }
      delete edge;
    }

    bool empty() const {
      return graph.size() == 0;
    }

    GraphEdge *pickStartEdge() {
      // Try and find a vertex from which there is only one outbound edge. Won't always succeed.
      for (graph_t::iterator i = graph.begin(); i != graph.end(); ++i) {
        GraphEdges &ge = i->second;
        if (ge.edges->next == NULL) {
          return ge.edges;
        }
      }
      return (*graph.begin()).second.edges;
    }

    GraphEdge *outboundEdges(const poly_t::vertex_t *v) {
      return graph[v].edges;
    }
  };



  /** 
   * \brief Take a set of new edges and split a face based upon those edges.
   * 
   * @param[in] face The face to be split.
   * @param[in] edges 
   * @param[out] face_loops Output list of face loops
   * @param[out] hole_loops Output list of hole loops
   * @param vi 
   */
  static void splitFace(const poly_t::face_t *face,
                        const carve::csg::V2Set &edges,
                        std::list<std::vector<const poly_t::vertex_t *> > &face_loops,
                        std::list<std::vector<const poly_t::vertex_t *> > &hole_loops,
                        const carve::csg::VertexIntersections &vi) {
    Graph graph;

#if defined(DEBUG)
    std::cerr << "splitFace()" << " face=" << face << " face->vertices.size()=" << face->vertices.size() << " edges.size()=" << edges.size() << std::endl;
#endif

    for (carve::csg::V2Set::const_iterator
           i = edges.begin(), e = edges.end();
         i != e;
         ++i) {
      const poly_t::vertex_t *v1 = ((*i).first), *v2 = ((*i).second);
      if (carve::geom::equal(v1->v, v2->v)) std::cerr << "WARNING! " << v1->v << "==" << v2->v << std::endl;
      graph.addEdge(v1, v2);
    }

    graph.computeProjection(face);
#if defined(DEBUG)
    graph.print(&vi);
#endif

    while (!graph.empty()) {
      GraphEdge *edge;
      GraphEdge *start;
      start = edge = graph.pickStartEdge();

      edge->visited = 0;

      int len = 0;

      while (1) {
        double in_ang = M_PI + edge->ang;
        if (in_ang > M_TWOPI) in_ang -= M_TWOPI;

        GraphEdge *opts;
        GraphEdge *out = NULL;
        double best = M_TWOPI + 1.0;

        for (opts = graph.outboundEdges(edge->tgt); opts; opts = opts->next) {
          if (opts->tgt == edge->src) {
            if (out == NULL && opts->next == NULL) out = opts;
          } else {
            double out_ang = carve::math::ANG(in_ang - opts->ang);

            if (out == NULL || out_ang < best) {
              out = opts;
              best = out_ang;
            }
          }
        }

        ASSERT(out != NULL);

        edge->loop_next = out;

        if (out->visited >= 0) {
          while (start != out) {
            GraphEdge *e = start;
            start = start->loop_next;
            e->loop_next = NULL;
            e->visited = -1;
          }
          len = edge->visited - out->visited + 1;
          break;
        }

        out->visited = edge->visited + 1;
        edge = out;
      }

      std::vector<const poly_t::vertex_t *> loop(len);
      std::vector<carve::geom2d::P2> projected(len);

      edge = start;
      for (int i = 0; i < len; ++i) {
        GraphEdge *next = edge->loop_next;
        loop[i] = edge->src;
        projected[i] = graph.projection(edge->src);
        graph.removeEdge(edge);
        edge = next;
      }
#if defined(DEBUG)
      std::cerr << "===============================================" << std::endl;
      graph.print(&vi);
#endif
      ASSERT(edge == start);
#if defined(DEBUG)
      std::cerr << "signed area of loop: " << carve::geom2d::signedArea(projected) << std::endl;
#endif
      if (carve::geom2d::signedArea(projected) < 0) {
#if defined(DEBUG)
        std::cerr << "output face loop size: " << loop.size() << " : ";
        for (size_t i = 0; i < loop.size(); ++i) std::cerr << " " << loop[i];
        std::cerr << std::endl;
#endif
        face_loops.push_back(std::vector<const poly_t::vertex_t *>());
        face_loops.back().swap(loop);
      } else {
#if defined(DEBUG)
        std::cerr << "output hole loop size: " << loop.size() << " : ";
        for (size_t i = 0; i < loop.size(); ++i) std::cerr << " " << loop[i];
        std::cerr << std::endl;
#endif
        hole_loops.push_back(std::vector<const poly_t::vertex_t *>());
        hole_loops.back().swap(loop);
      }
    }
#if defined(DEBUG)
    std::cerr << "===============================================" << std::endl;

    std::cerr << "result: " << face_loops.size() << " face loops (";
    for (std::list<std::vector<const poly_t::vertex_t *> >::const_iterator i = face_loops.begin(); i != face_loops.end(); ++i) {
      std::cerr << ((i != face_loops.begin()) ? " " : "") << (*i).size();
      for (unsigned j = 0; j < (*i).size(); ++j) {
        if (std::find((*i).begin() + j + 1, (*i).end(), (*i)[j]) != (*i).end()) {
          std::cerr << "[!]";
          break;
        }
      }
    }
    std::cerr << ") " << hole_loops.size() << " hole loops (";
    for (std::list<std::vector<const poly_t::vertex_t *> >::const_iterator i = hole_loops.begin(); i != hole_loops.end(); ++i) {
      std::cerr << ((i != hole_loops.begin()) ? " " : "") << (*i).size();
      for (unsigned j = 0; j < (*i).size(); ++j) {
        if (std::find((*i).begin() + j + 1, (*i).end(), (*i)[j]) != (*i).end()) {
          std::cerr << "[!]";
          break;
        }
      }
    }
    std::cerr << ")" << std::endl;
#endif
  }



  /** 
   * \brief Determine the relationship between a face loop and a hole loop.
   * 
   * Determine whether a face and hole share an edge, or a vertex,
   * or do not touch. Find a hole vertex that is not part of the
   * face, and a hole,face vertex pair that are coincident, if such
   * a pair exists.
   *
   * @param[in] f A face loop.
   * @param[in] f_sort A vector indexing \a f in address order
   * @param[in] h A hole loop.
   * @param[in] h_sort A vector indexing \a h in address order
   * @param[out] f_idx Index of a face vertex that is shared with the hole.
   * @param[out] h_idx Index of the hole vertex corresponding to \a f_idx.
   * @param[out] unmatched_h_idx Index of a hole vertex that is not part of the face.
   * @param[out] shares_vertex Boolean indicating that the face and the hole share a vertex.
   * @param[out] shares_edge Boolean indicating that the face and the hole share an edge.
   */
  static void compareFaceLoopAndHoleLoop(const std::vector<const poly_t::vertex_t *> &f,
                                         const std::vector<unsigned> &f_sort,
                                         const std::vector<const poly_t::vertex_t *> &h,
                                         const std::vector<unsigned> &h_sort,
                                         unsigned &f_idx,
                                         unsigned &h_idx,
                                         int &unmatched_h_idx,
                                         bool &shares_vertex,
                                         bool &shares_edge) {
    const size_t F = f.size();
    const size_t H = h.size();

    shares_vertex = shares_edge = false;
    unmatched_h_idx = -1;

    unsigned I, J;
    for (I = J = 0; I < F && J < H;) {
      unsigned i = f_sort[I], j = h_sort[J];
      if (f[i] == h[j]) {
        shares_vertex = true;
        f_idx = i;
        h_idx = j;
        if (f[(i + F - 1) % F] == h[(j + 1) % H]) {
          shares_edge = true;
          return;
        }
        ++I; ++J;
      } else if (f[i] < h[j]) {
        ++I;
      } else {
        unmatched_h_idx = j;
        ++J;
      }
    }
  }



  /** 
   * \brief Compute an embedding for a set of face loops and hole loops.
   *
   * Because face and hole loops may be contained within each other,
   * it must be determined which hole loops are directly contained
   * within a face loop.
   * 
   * @param[in] face The face from which these face and hole loops derive.
   * @param[in] face_loops 
   * @param[in] hole_loops 
   * @param[out] containing_faces     A vector which for each hole loop
   *                                  lists the indices of the face
   *                                  loops it is containined in.
   * @param[out] hole_shared_vertices A map from a face,hole pair to
   *                                  a shared vertex pair.
   */
  static void computeContainment(const poly_t::face_t *face,
                                 std::vector<std::vector<const poly_t::vertex_t *> > &face_loops,
                                 std::vector<std::vector<const poly_t::vertex_t *> > &hole_loops,
                                 std::vector<std::vector<int> > &containing_faces,
                                 std::map<int, std::map<int, std::pair<unsigned, unsigned> > > &hole_shared_vertices) {
    unsigned m, n;

    std::vector<std::vector<carve::geom2d::P2> > face_loops_projected, hole_loops_projected;
    std::vector<std::vector<unsigned> > face_loops_sorted, hole_loops_sorted;

    std::vector<double> face_loop_areas, hole_loop_areas;

    face_loops_projected.resize(face_loops.size());
    face_loops_sorted.resize(face_loops.size());
    face_loop_areas.resize(face_loops.size());

    hole_loops.resize(hole_loops.size());
    hole_loops_projected.resize(hole_loops.size());
    hole_loops_sorted.resize(hole_loops.size());
    hole_loop_areas.resize(hole_loops.size());

    // produce a projection of each face loop onto a 2D plane, and a
    // index vector which sorts vertices by address.
    for (size_t m = 0; m < face_loops.size(); ++m) {
      const std::vector<const poly_t::vertex_t *> &f_loop = (face_loops[m]);
      face_loops_projected[m].reserve(f_loop.size());
      face_loops_sorted[m].reserve(f_loop.size());
      for (n = 0; n < f_loop.size(); ++n) {
        face_loops_projected[m].push_back(carve::poly::face::project(face, f_loop[n]->v));
        face_loops_sorted[m].push_back(n);
      }
      face_loop_areas.push_back(carve::geom2d::signedArea(face_loops_projected[m]));
      std::sort(face_loops_sorted[m].begin(), face_loops_sorted[m].end(), 
                carve::index_sort<std::vector<const poly_t::vertex_t *> >(face_loops[m]));
    }

    // produce a projection of each hole loop onto a 2D plane, and a
    // index vector which sorts vertices by address.
    for (size_t m = 0; m < hole_loops.size(); ++m) {
      const std::vector<const poly_t::vertex_t *> &h_loop = (hole_loops[m]);
      hole_loops_projected[m].reserve(h_loop.size());
      hole_loops_projected[m].reserve(h_loop.size());
      for (n = 0; n < h_loop.size(); ++n) {
        hole_loops_projected[m].push_back(carve::poly::face::project(face, h_loop[n]->v));
        hole_loops_sorted[m].push_back(n);
      }
      hole_loop_areas.push_back(carve::geom2d::signedArea(hole_loops_projected[m]));
      std::sort(hole_loops_sorted[m].begin(), hole_loops_sorted[m].end(), 
                carve::index_sort<std::vector<const poly_t::vertex_t *> >(hole_loops[m]));
    }

    containing_faces.resize(hole_loops.size());

    for (unsigned i = 0; i < hole_loops.size(); ++i) {

      for (unsigned j = 0; j < face_loops.size(); ++j) {
        unsigned f_idx, h_idx;
        int unmatched_h_idx;
        bool shares_vertex, shares_edge;
        compareFaceLoopAndHoleLoop(face_loops[j],
                                   face_loops_sorted[j],
                                   hole_loops[i],
                                   hole_loops_sorted[i],
                                   f_idx, h_idx,
                                   unmatched_h_idx,
                                   shares_vertex,
                                   shares_edge);

#if defined(DEBUG)
        std::cerr << "face: " << j
                  << " hole: " << i
                  << " shares_vertex: " << shares_vertex
                  << " shares_edge: " << shares_edge
                  << std::endl;
#endif

        carve::geom3d::Vector test = hole_loops[i][0]->v;
        carve::geom2d::P2 test_p = carve::poly::face::project(face, test);

        if (shares_vertex) {
          hole_shared_vertices[i][j] = std::make_pair(h_idx, f_idx);
          // Hole touches face. Should be able to connect it up
          // trivially. Still need to record its containment, so that
          // the assignment below works.
          if (unmatched_h_idx != -1) {
#if defined(DEBUG)
            std::cerr << "using unmatched vertex: " << unmatched_h_idx << std::endl;
#endif
            test = hole_loops[i][unmatched_h_idx]->v;
            test_p = carve::poly::face::project(face, test);
          } else {
            // XXX: hole shares ALL vertices with face. Pick a point
            // internal to the projected poly.
            if (shares_edge) {
              // Hole shares edge with face => face can't contain hole.
              continue;
            }

            // XXX: how is this possible? Doesn't share an edge, but
            // also doesn't have any vertices that are not in
            // common. Degenerate hole?

            // XXX: come up with a test case for this.
            ASSERT(!!!"implement me");
          }
        }


        // XXX: use loop area to avoid some point-in-poly tests? Loop
        // area is faster, but not sure which is more robust.
        if (carve::geom2d::pointInPolySimple(face_loops_projected[j], test_p)) {
#if defined(DEBUG)
          std::cerr << "contains: " << i << " - " << j << std::endl;
#endif
          containing_faces[i].push_back(j);
        } else {
#if defined(DEBUG)
          std::cerr << "does not contain: " << i << " - " << j << std::endl;
#endif
        }
      }

#if defined(DEBUG)
      if (containing_faces[i].size() == 0) {
        //HOOK(drawFaceLoopWireframe(hole_loops[i], face->normal, 1.0, 0.0, 0.0, 1.0););
        std::cerr << "hole loop: ";
        for (unsigned j = 0; j < hole_loops[i].size(); ++j) {
          std::cerr << " " << hole_loops[i][j] << ":" << hole_loops[i][j]->v;
        }
        std::cerr << std::endl;
        for (unsigned j = 0; j < face_loops.size(); ++j) {
          //HOOK(drawFaceLoopWireframe(face_loops[j], face->normal, 0.0, 1.0, 0.0, 1.0););
        }
      }
#endif

      // ASSERT(containing_faces[i].size() >= 1);
    }
  }



  /** 
   * \brief Merge face loops and hole loops to produce a set of face loops without holes.
   * 
   * @param[in] face The face from which these face loops derive.
   * @param[in,out] f_loops A list of face loops.
   * @param[in] h_loops A list of hole loops to be incorporated into face loops.
   */
  static void mergeFacesAndHoles(const poly_t::face_t *face,
                                 std::list<std::vector<const poly_t::vertex_t *> > &f_loops,
                                 std::list<std::vector<const poly_t::vertex_t *> > &h_loops,
                                 carve::csg::CSG::Hooks &hooks) {
    std::vector<std::vector<const poly_t::vertex_t *> > face_loops;
    std::vector<std::vector<const poly_t::vertex_t *> > hole_loops;

    std::vector<std::vector<int> > containing_faces;
    std::map<int, std::map<int, std::pair<unsigned, unsigned> > > hole_shared_vertices;

    {
      // move input face and hole loops to temp vectors.
      size_t m;
      face_loops.resize(f_loops.size());
      m = 0;
      for (std::list<std::vector<const poly_t::vertex_t *> >::iterator
             i = f_loops.begin(), ie = f_loops.end();
           i != ie;
           ++i, ++m) {
        face_loops[m].swap((*i));
      }

      hole_loops.resize(h_loops.size());
      m = 0;
      for (std::list<std::vector<const poly_t::vertex_t *> >::iterator
             i = h_loops.begin(), ie = h_loops.end();
           i != ie;
           ++i, ++m) {
        hole_loops[m].swap((*i));
      }
      f_loops.clear();
      h_loops.clear();
    }

    // work out the embedding of holes and faces.
    computeContainment(face, face_loops, hole_loops, containing_faces, hole_shared_vertices);

    int unassigned = (int)hole_loops.size();

    std::vector<std::vector<int> > face_holes;
    face_holes.resize(face_loops.size());

    for (unsigned i = 0; i < containing_faces.size(); ++i) {
      if (containing_faces[i].size() == 0) {
        std::map<int, std::map<int, std::pair<unsigned, unsigned> > >::iterator it = hole_shared_vertices.find(i);
        if (it != hole_shared_vertices.end()) {
          std::map<int, std::pair<unsigned, unsigned> >::iterator it2 = (*it).second.begin();
          int f = (*it2).first;
          unsigned h_idx = (*it2).second.first;
          unsigned f_idx = (*it2).second.second;

          // patch the hole into the face directly. because
          // f_loop[f_idx] == h_loop[h_idx], we don't need to
          // duplicate the f_loop vertex.

          std::vector<const poly_t::vertex_t *> &f_loop = face_loops[f];
          std::vector<const poly_t::vertex_t *> &h_loop = hole_loops[i];

          f_loop.insert(f_loop.begin() + f_idx + 1, h_loop.size(), NULL);

          unsigned p = f_idx + 1;
          for (unsigned a = h_idx + 1; a < h_loop.size(); ++a, ++p) {
            f_loop[p] = h_loop[a];
          }
          for (unsigned a = 0; a <= h_idx; ++a, ++p) {
            f_loop[p] = h_loop[a];
          }

#if defined(DEBUG)
          std::cerr << "hook face " << f << " to hole " << i << "(vertex)" << std::endl;
#endif
        } else {
          std::cerr << "uncontained hole loop does not share vertices with any face loop!" << std::endl;
        }
        unassigned--;
      }
    }


    // work out which holes are directly contained within which faces.
    while (unassigned) {
      std::set<int> removed;

      for (unsigned i = 0; i < containing_faces.size(); ++i) {
        if (containing_faces[i].size() == 1) {
          int f = containing_faces[i][0];
          face_holes[f].push_back(i);
#if defined(DEBUG)
          std::cerr << "hook face " << f << " to hole " << i << std::endl;
#endif
          removed.insert(f);
          unassigned--;
        }
      }
      for (std::set<int>::iterator f = removed.begin(); f != removed.end(); ++f) {
        for (unsigned i = 0; i < containing_faces.size(); ++i) {
          containing_faces[i].erase(std::remove(containing_faces[i].begin(),
                                                containing_faces[i].end(),
                                                *f),
                                    containing_faces[i].end());
        }
      }
    }

    // patch holes into faces.
    for (unsigned i = 0; i < face_loops.size(); ++i) {
      std::vector<std::vector<const poly_t::vertex_t *> > face_hole_loops;
      face_hole_loops.resize(face_holes[i].size());
      for (unsigned j = 0; j < face_holes[i].size(); ++j) {
        face_hole_loops[j].swap(hole_loops[face_holes[i][j]]);
      }
      if (face_hole_loops.size()) {

        f_loops.push_back(carve::triangulate::incorporateHolesIntoPolygon(carve::poly::p2_adapt_project<3>(face->project), face_loops[i], face_hole_loops));
      } else {
        f_loops.push_back(face_loops[i]);
      }
    }
  }



  /** 
   * \brief Assemble the base loop for a face.
   *
   * The base loop is the original face loop, including vertices
   * created by intersections crossing any of its edges.
   * 
   * @param[in] face The face to process.
   * @param[in] vmap 
   * @param[in] face_split_edges 
   * @param[in] divided_edges A mapping from edge pointer to sets of
   *            ordered vertices corrsponding to the intersection points
   *            on that edge.
   * @param[out] base_loop A vector of the vertices of the base loop.
   */
  static void assembleBaseLoop(const poly_t::face_t *face,
                               const carve::csg::VVMap &vmap,
                               const carve::csg::FV2SMap &face_split_edges,
                               const carve::csg::EVVMap &divided_edges,
                               std::vector<const poly_t::vertex_t *> &base_loop) {
    base_loop.clear();

    // XXX: assumes that face->edges is in the same order as
    // face->vertices. (Which it is)
    for (size_t j = 0, je = face->vertices.size(); j < je; ++j) {
      base_loop.push_back(carve::csg::map_vertex(vmap, face->vertices[j]));

      const poly_t::edge_t *e = face->edges[j];
      carve::csg::EVVMap::const_iterator ev = divided_edges.find(e);

      if (ev != divided_edges.end()) {
        const std::vector<const poly_t::vertex_t *> &ev_vec = ((*ev).second);

        if (e->v1 == face->vertices[j]) {
          // edge is forward;
          for (size_t k = 0, ke = ev_vec.size(); k < ke;) {
            base_loop.push_back(ev_vec[k++]);
          }
        } else {
          // edge is backward;
          for (size_t k = ev_vec.size(); k;) {
            base_loop.push_back(ev_vec[--k]);
          }
        }
      }
    }
  }

}

/** 
 * \brief Build a set of face loops for all (split) faces of a Polyhedron.
 * 
 * @param[in] poly The polyhedron to process
 * @param vmap 
 * @param face_split_edges 
 * @param divided_edges 
 * @param[out] face_loops_out The resulting face loops
 * 
 * @return The number of edges generated.
 */
size_t carve::csg::CSG::generateFaceLoops(const poly_t *poly,
                                          const VVMap &vmap,
                                          const FV2SMap &face_split_edges,
                                          const EVVMap &divided_edges,
                                          FaceLoopList &face_loops_out) {
  static carve::TimingName FUNC_NAME("CSG::generateFaceLoops()");
  carve::TimingBlock block(FUNC_NAME);
  size_t generated_edges = 0;
  std::vector<const poly_t::vertex_t *> base_loop;
  std::list<std::vector<const poly_t::vertex_t *> > face_loops, hole_loops;
  
  for (std::vector<poly_t::face_t >::const_iterator
         i = poly->faces.begin(), e = poly->faces.end();
       i != e;
       ++i) {
    const poly_t::face_t *face = &(*i);

    assembleBaseLoop(face, vmap, face_split_edges, divided_edges, base_loop);

    FV2SMap::const_iterator fse_iter = face_split_edges.find(face);
    if (fse_iter != face_split_edges.end()) {
      // complex case: input face is split into multiple output faces.
      carve::csg::V2Set face_edges;

      // collect the perimeter edges
#if defined(DEBUG)
      std::cerr << "base loop is:";
      for (size_t j = 0, je = base_loop.size(); j < je; ++j) {
        std::cerr << " " << base_loop[j];
      }
      std::cerr << std::endl;
#endif

      for (size_t j = 0, je = base_loop.size() - 1; j < je; ++j) {
        face_edges.insert(std::make_pair(base_loop[j], base_loop[j + 1]));
      }
      face_edges.insert(std::make_pair(base_loop[base_loop.size() - 1], base_loop[0]));

      // collect the split edges (as long as they're not on the perimeter)
      const FV2SMap::mapped_type &fse = ((*fse_iter).second);
      bool actually_split = false;
      size_t added_edges = 0;
      for (FV2SMap::mapped_type::const_iterator
             j = fse.begin(), je =  fse.end();
           j != je;
           ++j) {
        const poly_t::vertex_t *v1 = ((*j).first), *v2 = ((*j).second);
#if defined(DEBUG)
        std::cerr << "testing edge: " << v1 << " - " << v2 << std::endl;
#endif

        if (face_edges.find(std::make_pair(v1, v2)) == face_edges.end() &&
            face_edges.find(std::make_pair(v2, v1)) == face_edges.end()) {
#if defined(DEBUG)
          std::cerr << "  adding edge: " << v1 << " - " << v2 << std::endl;
#endif
          actually_split = true;
          added_edges++;
          face_edges.insert(std::make_pair(v1, v2));
          face_edges.insert(std::make_pair(v2, v1));
        }
      }

      // if there's at least one non-perimeter edge, then the face is split.
      if (actually_split) {
        hole_loops.clear();
        face_loops.clear();

#if defined(DEBUG)
        std::cerr << std::endl << "splitting face " << face << " (" << added_edges << " added edges)" << std::endl;
#endif

        // trace the edge graph to segregate edges into face and hole loops.
        splitFace(face, face_edges, face_loops, hole_loops, vertex_intersections);

        // if there are holes, then they need to be merged with faces.
        if (hole_loops.size()) {
#if defined(DEBUG)
          std::cerr << "before split: "
                    << face_loops.size() << " face loops "
                    << hole_loops.size() << " hole loops"
                    << std::endl;
#endif

#if defined(DEBUG_WRITE_PLY_DATA)
          {
            std::map<const poly_t::vertex_t *, size_t> v_included;

            for (std::list<std::vector<const poly_t::vertex_t *> >::iterator
                   i = face_loops.begin(); i != face_loops.end(); ++i) {
              for (size_t j = 0; j < (*i).size(); ++j) {
                if (v_included.find((*i)[j]) == v_included.end()) {
                  size_t &p = v_included[(*i)[j]];
                  p = v_included.size() - 1;
                }
              }
            }

            for (std::list<std::vector<const poly_t::vertex_t *> >::iterator
                   i = hole_loops.begin(); i != hole_loops.end(); ++i) {
              for (size_t j = 0; j < (*i).size(); ++j) {
                if (v_included.find((*i)[j]) == v_included.end()) {
                  size_t &p = v_included[(*i)[j]];
                  p = v_included.size() - 1;
                }
              }
            }

            carve::line::PolylineSet fh;
            fh.vertices.resize(v_included.size());
            for (std::map<const poly_t::vertex_t *, size_t>::const_iterator
                   i = v_included.begin(); i != v_included.end(); ++i) {
              fh.vertices[(*i).second].v = (*i).first->v;
            }

            {
              std::vector<size_t> connected;
              for (std::list<std::vector<const poly_t::vertex_t *> >::iterator
                     i = face_loops.begin(); i != face_loops.end(); ++i) {
                connected.clear();
                for (size_t j = 0; j < (*i).size(); ++j) {
                  connected.push_back(v_included[(*i)[j]]);
                }
                fh.addPolyline(true, connected.begin(), connected.end());
              }
              for (std::list<std::vector<const poly_t::vertex_t *> >::iterator
                     i = hole_loops.begin(); i != hole_loops.end(); ++i) {
                connected.clear();
                for (size_t j = 0; j < (*i).size(); ++j) {
                  connected.push_back(v_included[(*i)[j]]);
                }
                fh.addPolyline(true, connected.begin(), connected.end());
              }
            }

            void writePLY(std::string &out_file, const carve::line::PolylineSet *lines, bool ascii);
            std::string out("/tmp/hole_merge.ply");
            writePLY(out, &fh, true);
          }
#endif

          mergeFacesAndHoles(face, face_loops, hole_loops, hooks);

#if defined(DEBUG)
          std::cerr << "after split: "
                    << face_loops.size() << " face loops "
                    << hole_loops.size() << " hole loops"
                    << std::endl;
#endif
        }
      } else {
        // in this case, the intersections just add vertices to
        // the edge of the face.
        face_loops.clear();
        face_loops.push_back(base_loop);
      }
    } else {
      // simple case: input face is output face (possibly with the
      // addition of vertices at intersections).

      face_loops.clear();
      face_loops.push_back(base_loop);
    }

    // now record all the resulting face loops.
    for (std::list<std::vector<const poly_t::vertex_t *> >::const_iterator
           f = face_loops.begin(), fe = face_loops.end();
         f != fe;
         ++f) {
      face_loops_out.append(new FaceLoop(face, *f));
      generated_edges += (*f).size();
    }
  }
  return generated_edges;
}
