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

#include <carve/mesh.hpp>
#include <carve/mesh_impl.hpp>

#include <carve/poly.hpp>

namespace {
  inline double CALC_X(const carve::geom::plane<3> &p, double y, double z) { return -(p.d + p.N.y * y + p.N.z * z) / p.N.x; }
  inline double CALC_Y(const carve::geom::plane<3> &p, double x, double z) { return -(p.d + p.N.x * x + p.N.z * z) / p.N.y; }
  inline double CALC_Z(const carve::geom::plane<3> &p, double x, double y) { return -(p.d + p.N.x * x + p.N.y * y) / p.N.z; }

  carve::geom::vector<2> _project_1(const carve::geom::vector<3> &v) {
    return carve::geom::VECTOR(v.z, v.y);
  }

  carve::geom::vector<2> _project_2(const carve::geom::vector<3> &v) {
    return carve::geom::VECTOR(v.x, v.z);
  }

  carve::geom::vector<2> _project_3(const carve::geom::vector<3> &v) {
    return carve::geom::VECTOR(v.y, v.x);
  }

  carve::geom::vector<2> _project_4(const carve::geom::vector<3> &v) {
    return carve::geom::VECTOR(v.y, v.z);
  }

  carve::geom::vector<2> _project_5(const carve::geom::vector<3> &v) {
    return carve::geom::VECTOR(v.z, v.x);
  }

  carve::geom::vector<2> _project_6(const carve::geom::vector<3> &v) {
    return carve::geom::VECTOR(v.x, v.y);
  }
  
  carve::geom::vector<3> _unproject_1(const carve::geom::vector<2> &p, const carve::geom3d::Plane &plane) {
    return carve::geom::VECTOR(CALC_X(plane, p.y, p.x), p.y, p.x);
  }

  carve::geom::vector<3> _unproject_2(const carve::geom::vector<2> &p, const carve::geom3d::Plane &plane) {
    return carve::geom::VECTOR(p.x, CALC_Y(plane, p.x, p.y), p.y);
  }

  carve::geom::vector<3> _unproject_3(const carve::geom::vector<2> &p, const carve::geom3d::Plane &plane) {
    return carve::geom::VECTOR(p.y, p.x, CALC_Z(plane, p.y, p.x));
  }

  carve::geom::vector<3> _unproject_4(const carve::geom::vector<2> &p, const carve::geom3d::Plane &plane) {
    return carve::geom::VECTOR(CALC_X(plane, p.x, p.y), p.x, p.y);
  }

  carve::geom::vector<3> _unproject_5(const carve::geom::vector<2> &p, const carve::geom3d::Plane &plane) {
    return carve::geom::VECTOR(p.y, CALC_Y(plane, p.y, p.x), p.x);
  }

  carve::geom::vector<3> _unproject_6(const carve::geom::vector<2> &p, const carve::geom3d::Plane &plane) {
    return carve::geom::VECTOR(p.x, p.y, CALC_Z(plane, p.x, p.y));
  }

  static carve::geom::vector<2> (*project_tab[2][3])(const carve::geom::vector<3> &) = {
    { &_project_1, &_project_2, &_project_3 },
    { &_project_4, &_project_5, &_project_6 }
  };

  static carve::geom::vector<3> (*unproject_tab[2][3])(const carve::geom::vector<2> &, const carve::geom3d::Plane &) = {
    { &_unproject_1, &_unproject_2, &_unproject_3 },
    { &_unproject_4, &_unproject_5, &_unproject_6 }
  };

}

namespace carve {
  namespace mesh {

    namespace detail {

      void FaceStitcher::matchSimpleEdges() {
        // join faces that share an edge, where no other faces are incident.
        for (edge_map_t::iterator i = edges.begin(); i != edges.end(); ++i) {
          const vpair_t &ev = (*i).first;
          edge_map_t::iterator j = edges.find(vpair_t(ev.second, ev.first));
          if (j == edges.end()) {
            for (edgelist_t::iterator k = (*i).second.begin(); k != (*i).second.end(); ++k) {
              is_open[ (*k)->face->id] = true;
            }
          } else if ((*i).second.size() != 1 || (*j).second.size() != 1) {
            std::swap(complex_edges[(*i).first], (*i).second);
          } else {
            // simple edge.
            edge_t *a = (*i).second.front();
            edge_t *b = (*j).second.front();
            if (a < b) {
              // every simple edge pair is encountered twice. only merge once.
              a->rev = b;
              b->rev = a;
              face_groups.merge_sets(a->face->id, b->face->id);
            }
          }
        }
      }

      size_t FaceStitcher::faceGroupID(const Face<3> *face) {
        return face_groups.find_set_head(face->id);
      }

      size_t FaceStitcher::faceGroupID(const Edge<3> *edge) {
        return face_groups.find_set_head(edge->face->id);
      }

      void FaceStitcher::orderForwardAndReverseEdges(std::vector<std::vector<Edge<3> *> > &efwd,
                                                     std::vector<std::vector<Edge<3> *> > &erev,
                                                     std::vector<std::vector<EdgeOrderData> > &result) {
        const size_t Nfwd = efwd.size();
        const size_t Nrev = erev.size();
        const size_t N = efwd[0].size();

        result.resize(N);

        for (size_t i = 0; i < N; ++i) {
          Edge<3> *base = efwd[0][i];

          result[i].reserve(Nfwd + Nrev);
          for (size_t j = 0; j < Nfwd; ++j) {
            result[i].push_back(EdgeOrderData(efwd[j][i], j, false));
            CARVE_ASSERT(efwd[0][i]->v1() == efwd[j][i]->v1());
            CARVE_ASSERT(efwd[0][i]->v2() == efwd[j][i]->v2());
          }
          for (size_t j = 0; j < Nrev; ++j) {
            result[i].push_back(EdgeOrderData(erev[j][i], j, true));
            CARVE_ASSERT(erev[0][i]->v1() == erev[j][i]->v1());
            CARVE_ASSERT(erev[0][i]->v2() == erev[j][i]->v2());
          }

          std::sort(result[i].begin(),
                    result[i].end(),
                    EdgeOrderData::Cmp(base->v2()->v - base->v1()->v, result[i][0].face_dir));
        }
      }

      void FaceStitcher::edgeIncidentGroups(const vpair_t &e,
                                            const edge_map_t &all_edges,
                                            std::pair<std::set<size_t>, std::set<size_t> > &groups) {
        groups.first.clear();
        groups.second.clear();
        edge_map_t::const_iterator i;

        i = all_edges.find(e);
        if (i != all_edges.end()) {
          for (edgelist_t::const_iterator j = (*i).second.begin(); j != (*i).second.end(); ++j) {
            groups.first.insert(faceGroupID(*j));
          }
        }

        i = all_edges.find(vpair_t(e.second, e.first));
        if (i != all_edges.end()) {
          for (edgelist_t::const_iterator j = (*i).second.begin(); j != (*i).second.end(); ++j) {
            groups.second.insert(faceGroupID(*j));
          }
        }
      }

      void FaceStitcher::buildEdgeGraph(const edge_map_t &all_edges) {
        for (edge_map_t::const_iterator i = all_edges.begin();
             i != all_edges.end();
             ++i) {
          edge_graph[(*i).first.first].insert((*i).first.second);
        }
      }

      void FaceStitcher::extractPath(std::vector<const vertex_t *> &path) {
        path.clear();

        edge_graph_t::iterator iter = edge_graph.begin();

        
        const vertex_t *init = (*iter).first;
        const vertex_t *next = *(*iter).second.begin();
        const vertex_t *prev = NULL;
        const vertex_t *vert = init;

        while ((*iter).second.size() == 2) {
          prev = *std::find_if((*iter).second.begin(),
                               (*iter).second.end(),
                               std::bind2nd(std::not_equal_to<const vertex_t *>(), next));
          next = vert;
          vert = prev;
          iter = edge_graph.find(vert);
          CARVE_ASSERT(iter != edge_graph.end());
          if (vert == init) break;
        }
        init = vert;

        std::vector<const edge_t *> efwd;
        std::vector<const edge_t *> erev;

        edge_map_t::iterator edgeiter;
        edgeiter = complex_edges.find(vpair_t(vert, next));
        std::copy((*edgeiter).second.begin(), (*edgeiter).second.end(), std::back_inserter(efwd));

        edgeiter = complex_edges.find(vpair_t(next, vert));
        std::copy((*edgeiter).second.begin(), (*edgeiter).second.end(), std::back_inserter(erev));

        path.push_back(vert);

        prev = vert;
        vert = next;
        path.push_back(vert);
        iter = edge_graph.find(vert);
        CARVE_ASSERT(iter != edge_graph.end());

        while (vert != init && (*iter).second.size() == 2) {
          next = *std::find_if((*iter).second.begin(),
                               (*iter).second.end(),
                               std::bind2nd(std::not_equal_to<const vertex_t *>(), prev));

          edgeiter = complex_edges.find(vpair_t(vert, next));
          if ((*edgeiter).second.size() != efwd.size()) goto done;

          for (size_t i = 0; i < efwd.size(); ++i) {
            Edge<3> *e_next = efwd[i]->perimNext();
            if (e_next->v2() != next) goto done;
            efwd[i] = e_next;
          }

          edgeiter = complex_edges.find(vpair_t(next, vert));
          if ((*edgeiter).second.size() != erev.size()) goto done;

          for (size_t i = 0; i < erev.size(); ++i) {
            Edge<3> *e_prev = erev[i]->perimPrev();
            if (e_prev->v1() != next) goto done;
            erev[i] = e_prev;
          }

          prev = vert;
          vert = next;
          path.push_back(vert);
          iter = edge_graph.find(vert);
          CARVE_ASSERT(iter != edge_graph.end());
        }
      done:;
      }

      void FaceStitcher::removePath(const std::vector<const vertex_t *> &path) {
        for (size_t i = 1; i < path.size() - 1; ++i) {
          edge_graph.erase(path[i]);
        }

        edge_graph[path[0]].erase(path[1]);
        if (edge_graph[path[0]].size() == 0) {
          edge_graph.erase(path[0]);
        }

        edge_graph[path[path.size()-1]].erase(path[path.size()-2]);
        if (edge_graph[path[path.size()-1]].size() == 0) {
          edge_graph.erase(path[path.size()-1]);
        }
      }

      void FaceStitcher::reorder(std::vector<EdgeOrderData> &ordering,
                                 size_t grp) {
        if (!ordering[0].is_reversed && ordering[0].group_id == grp) return;
        for (size_t i = 1; i < ordering.size(); ++i) {
          if (!ordering[i].is_reversed && ordering[i].group_id == grp) {
            std::vector<EdgeOrderData> temp;
            temp.reserve(ordering.size());
            std::copy(ordering.begin() + i, ordering.end(), std::back_inserter(temp));
            std::copy(ordering.begin(), ordering.begin() + i, std::back_inserter(temp));
            std::copy(temp.begin(), temp.end(), ordering.begin());
            return;
          }
        }
      }

      struct lt_second {
        template<typename pair_t>
        bool operator()(const pair_t &a, const pair_t &b) const {
          return a.second < b.second;
        }
      };

      void FaceStitcher::fuseEdges(std::vector<Edge<3> *> &fwd,
                                   std::vector<Edge<3> *> &rev) {
        for (size_t i = 0; i < fwd.size(); ++i) {
          fwd[i]->rev = rev[i];
          rev[i]->rev = fwd[i];
          face_groups.merge_sets(fwd[i]->face->id, rev[i]->face->id);
        }
      }

      void FaceStitcher::joinGroups(std::vector<std::vector<Edge<3> *> > &efwd,
                                    std::vector<std::vector<Edge<3> *> > &erev,
                                    size_t fwd_grp,
                                    size_t rev_grp) {
        fuseEdges(efwd[fwd_grp], erev[rev_grp]);
      }

      void FaceStitcher::matchOrderedEdges(const std::vector<std::vector<EdgeOrderData> >::iterator begin,
                                           const std::vector<std::vector<EdgeOrderData> >::iterator end,
                                           std::vector<std::vector<Edge<3> *> > &efwd,
                                           std::vector<std::vector<Edge<3> *> > &erev) {
        typedef std::unordered_map<std::pair<size_t, size_t>, size_t> pair_counts_t;
        while (1) {
          pair_counts_t pair_counts;

          for (std::vector<std::vector<EdgeOrderData> >::iterator i = begin; i != end; ++i) {
            std::vector<EdgeOrderData> &e = *i;
            for (size_t j = 0; j < e.size(); ++j) {
              if (!e[j].is_reversed && e[(j+1)%e.size()].is_reversed) {
                pair_counts[std::make_pair(e[j].group_id,
                                           e[(j+1)%e.size()].group_id)]++;
              }
            }
          }

          if (!pair_counts.size()) break;

          std::vector<std::pair<size_t, std::pair<size_t, size_t> > > counts;
          counts.reserve(pair_counts.size());
          for (pair_counts_t::iterator iter = pair_counts.begin(); iter != pair_counts.end(); ++iter) {
            counts.push_back(std::make_pair((*iter).second, (*iter).first));
          }
          std::make_heap(counts.begin(), counts.end());

          std::set<size_t> rem_fwd, rem_rev;

          while (counts.size()) {
            std::pair<size_t, size_t> join = counts.front().second;
            std::pop_heap(counts.begin(), counts.end());
            counts.pop_back();
            if (rem_fwd.find(join.first) != rem_fwd.end()) continue;
            if (rem_rev.find(join.second) != rem_rev.end()) continue;

            size_t g1 = join.first;
            size_t g2 = join.second;

            joinGroups(efwd, erev, g1, g2);

            for (std::vector<std::vector<EdgeOrderData> >::iterator i = begin; i != end; ++i) {
              (*i).erase(std::remove_if((*i).begin(), (*i).end(), EdgeOrderData::TestGroups(g1, g2)), (*i).end());
            }

            rem_fwd.insert(g1);
            rem_rev.insert(g2);
          }
        }
      }

      void FaceStitcher::resolveOpenEdges() {
        // Remove open regions of mesh. Doing this may make additional
        // edges simple (for example, removing a fin from the edge of
        // a cube), and may also expose more open mesh regions. In the
        // latter case, the process must be repeated to deal with the
        // newly uncovered regions.
        std::unordered_set<size_t> open_groups;

        for (size_t i = 0; i < is_open.size(); ++i) {
          if (is_open[i]) open_groups.insert(face_groups.find_set_head(i));
        }

        while (!open_groups.empty()) {
          std::list<vpair_t> edge_0, edge_1;

          for (edge_map_t::iterator i = complex_edges.begin(); i != complex_edges.end(); ++i) {
            bool was_modified = false;
            for(edgelist_t::iterator j = (*i).second.begin(); j != (*i).second.end(); ) {
              if (open_groups.find(faceGroupID(*j)) != open_groups.end()) {
                j = (*i).second.erase(j);
                was_modified = true;
              } else {
                ++j;
              }
            }
            if (was_modified) {
              if ((*i).second.empty()) {
                edge_0.push_back((*i).first);
              } else if ((*i).second.size() == 1) {
                edge_1.push_back((*i).first);
              }
            }
          }

          for (std::list<vpair_t>::iterator i = edge_1.begin(); i != edge_1.end(); ++i) {
            vpair_t e1 = *i;
            edge_map_t::iterator e1i = complex_edges.find(e1);
            if (e1i == complex_edges.end()) continue;
            vpair_t e2 = vpair_t(e1.second, e1.first);
            edge_map_t::iterator e2i = complex_edges.find(e2);
            CARVE_ASSERT(e2i != complex_edges.end()); // each complex edge should have a mate.

            if ((*e2i).second.size() == 1) {
              // merge newly simple edges, delete both from complex_edges.
              edge_t *a = (*e1i).second.front();
              edge_t *b = (*e2i).second.front();
              a->rev = b;
              b->rev = a;
              face_groups.merge_sets(a->face->id, b->face->id);
              complex_edges.erase(e1i);
              complex_edges.erase(e2i);
            }
          }

          open_groups.clear();

          for (std::list<vpair_t>::iterator i = edge_0.begin(); i != edge_0.end(); ++i) {
            vpair_t e1 = *i;
            edge_map_t::iterator e1i = complex_edges.find(e1);
            vpair_t e2 = vpair_t(e1.second, e1.first);
            edge_map_t::iterator e2i = complex_edges.find(e2);
            CARVE_ASSERT(e2i != complex_edges.end()); // each complex edge should have a mate.

            for (edgelist_t::iterator j = (*e2i).second.begin(); j != (*e2i).second.end(); ++j) {
              open_groups.insert(faceGroupID(*j));
            }
            complex_edges.erase(e1i);
            complex_edges.erase(e2i);
          }
        }
      }

      void FaceStitcher::extractConnectedEdges(std::vector<const vertex_t *>::iterator begin,
                                               std::vector<const vertex_t *>::iterator end,
                                               std::vector<std::vector<Edge<3> *> > &efwd,
                                               std::vector<std::vector<Edge<3> *> > &erev) {
        const size_t N = std::distance(begin, end) - 1;

        std::vector<const vertex_t *>::iterator e1, e2;
        e1 = e2 = begin; ++e2;
        vpair_t start_f = vpair_t(*e1, *e2);
        vpair_t start_r = vpair_t(*e2, *e1);

        const size_t Nfwd = complex_edges[start_f].size();
        const size_t Nrev = complex_edges[start_r].size();

        size_t j;
        edgelist_t::iterator ji;

        efwd.clear(); efwd.resize(Nfwd);
        erev.clear(); erev.resize(Nrev);

        for (j = 0, ji = complex_edges[start_f].begin();
             ji != complex_edges[start_f].end();
             ++j, ++ji) {
          efwd[j].reserve(N);
          efwd[j].push_back(*ji);
        }

        for (j = 0, ji = complex_edges[start_r].begin();
             ji != complex_edges[start_r].end();
             ++j, ++ji) {
          erev[j].reserve(N);
          erev[j].push_back(*ji);
        }

        std::vector<Edge<3> *> temp_f, temp_r;
        temp_f.resize(Nfwd);
        temp_r.resize(Nrev);

        for (j = 1; j < N; ++j) {
          ++e1; ++e2;
          vpair_t ef = vpair_t(*e1, *e2);
          vpair_t er = vpair_t(*e2, *e1);

          if (complex_edges[ef].size() != Nfwd || complex_edges[ef].size() != Nrev) break;

          for (size_t k = 0; k < Nfwd; ++k) {
            Edge<3> *e_next = efwd[k].back()->perimNext();
            CARVE_ASSERT(e_next == NULL || e_next->rev == NULL);
            if (e_next == NULL || e_next->v2() != *e2) goto done;
            CARVE_ASSERT(e_next->v1() == *e1);
            CARVE_ASSERT(std::find(complex_edges[ef].begin(), complex_edges[ef].end(), e_next) != complex_edges[ef].end());
            temp_f[k] = e_next;
          }

          for (size_t k = 0; k < Nrev; ++k) {
            Edge<3> *e_next = erev[k].back()->perimPrev();
            if (e_next == NULL || e_next->v1() != *e2) goto done;
            CARVE_ASSERT(e_next->v2() == *e1);
            CARVE_ASSERT(std::find(complex_edges[er].begin(), complex_edges[er].end(), e_next) != complex_edges[er].end());
            temp_r[k] = e_next;
          }

          for (size_t k = 0; k < Nfwd; ++k) {
            efwd[k].push_back(temp_f[k]);
          }

          for (size_t k = 0; k < Nrev; ++k) {
            erev[k].push_back(temp_r[k]);
          }
        }
      done:;
      }

      void FaceStitcher::construct() {
        matchSimpleEdges();
        if (!complex_edges.size()) return;

        resolveOpenEdges();
        if (!complex_edges.size()) return;

        buildEdgeGraph(complex_edges);

        std::list<std::vector<const vertex_t *> > paths;

        while (edge_graph.size()) {
          paths.push_back(std::vector<const vertex_t *>());
          extractPath(paths.back());
          removePath(paths.back());
        };


        for (std::list<std::vector<const vertex_t *> >::iterator path = paths.begin(); path != paths.end(); ++path) {
          for (size_t i = 0; i < (*path).size() - 1;) {
            std::vector<std::vector<Edge<3> *> > efwd, erev;

            extractConnectedEdges((*path).begin() + i, (*path).end(), efwd, erev);

            std::vector<std::vector<EdgeOrderData> > orderings;
            orderForwardAndReverseEdges(efwd, erev, orderings);

            matchOrderedEdges(orderings.begin(), orderings.end(), efwd, erev);
            i += efwd[0].size();
          }
        }
      }
    }




    template<unsigned ndim>
    typename Face<ndim>::project_t Face<ndim>::getProjector(bool positive_facing, int axis) {
      return NULL;
    }



    template<>
    Face<3>::project_t Face<3>::getProjector(bool positive_facing, int axis) {
      return project_tab[positive_facing ? 1 : 0][axis];
    }



    template<unsigned ndim>
    typename Face<ndim>::unproject_t Face<ndim>::getUnprojector(bool positive_facing, int axis) {
      return NULL;
    }



    template<>
    Face<3>::unproject_t Face<3>::getUnprojector(bool positive_facing, int axis) {
      return unproject_tab[positive_facing ? 1 : 0][axis];
    }

    template class Vertex<3>;
    template class Edge<3>;
    template class Face<3>;
    template class Mesh<3>;
    template class MeshSet<3>;
  }



  // construct a MeshSet from a Polyhedron, maintaining on the
  // connectivity information in the Polyhedron.
  mesh::MeshSet<3> *meshFromPolyhedron(const poly::Polyhedron *poly, int manifold_id) {
    typedef mesh::Vertex<3> vertex_t;
    typedef mesh::Vertex<3>::vector_t vector_t;
    typedef mesh::Edge<3> edge_t;
    typedef mesh::Face<3> face_t;
    typedef mesh::Mesh<3> mesh_t;
    typedef mesh::MeshSet<3> meshset_t;

    std::vector<vertex_t> vertex_storage;
    vertex_storage.reserve(poly->vertices.size());
    for (size_t i = 0; i < poly->vertices.size(); ++i) {
      vertex_storage.push_back(vertex_t(poly->vertices[i].v));
    }

    std::vector<std::vector<face_t *> > faces;
    faces.resize(poly->manifold_is_closed.size());

    std::unordered_map<std::pair<size_t, size_t>, std::list<edge_t *> > vertex_to_edge;

    std::vector<vertex_t *> vert_ptrs;
    for (size_t i = 0; i < poly->faces.size(); ++i) {
      const poly::Polyhedron::face_t &src = poly->faces[i];
      if (manifold_id != -1 && src.manifold_id != manifold_id) continue;
      vert_ptrs.clear();
      vert_ptrs.reserve(src.nVertices());
      for (size_t j = 0; j < src.nVertices(); ++j) {
        size_t vi = poly->vertexToIndex_fast(src.vertex(j));
        vert_ptrs.push_back(&vertex_storage[vi]);
      }
      face_t *face = new face_t(vert_ptrs.begin(), vert_ptrs.end());
      face->id = src.manifold_id;
      faces[src.manifold_id].push_back(face);

      edge_t *edge = face->edge;
      do {
        vertex_to_edge[std::make_pair(size_t(edge->v1() - &vertex_storage[0]),
                                      size_t(edge->v2() - &vertex_storage[0]))].push_back(edge);
        edge = edge->next;
      } while (edge != face->edge);
    }

    // copy connectivity from Polyhedron.
    for (size_t i = 0; i < poly->edges.size(); ++i) {
      const poly::Polyhedron::edge_t &src = poly->edges[i];
      size_t v1i = poly->vertexToIndex_fast(src.v1);
      size_t v2i = poly->vertexToIndex_fast(src.v2);

      std::list<edge_t *> &efwd = vertex_to_edge[std::make_pair(v1i, v2i)];
      std::list<edge_t *> &erev = vertex_to_edge[std::make_pair(v2i, v1i)];

      const std::vector<const poly::Polyhedron::face_t *> &facepairs = poly->connectivity.edge_to_face[i];
      for (size_t j = 0; j < facepairs.size(); j += 2) {
        const poly::Polyhedron::face_t *fa, *fb;
        fa = facepairs[j];
        fb = facepairs[j+1];
        if (!fa || !fb) continue;
        CARVE_ASSERT(fa->manifold_id == fb->manifold_id);
        if (manifold_id != -1 && fa->manifold_id != manifold_id) continue;

        std::list<edge_t *>::iterator efwdi, erevi;
        for (efwdi = efwd.begin(); efwdi != efwd.end() && (*efwdi)->face->id != fa->manifold_id; ++efwdi);
        for (erevi = erev.begin(); erevi != erev.end() && (*erevi)->face->id != fa->manifold_id; ++erevi);
        CARVE_ASSERT(efwdi != efwd.end() && erevi != erev.end());

        (*efwdi)->rev = (*erevi);
        (*erevi)->rev = (*efwdi);
      }
    }

    std::vector<mesh_t *> meshes;
    meshes.reserve(faces.size());
    for (size_t i = 0; i < faces.size(); ++i) {
      if (faces[i].size()) {
        meshes.push_back(new mesh_t(faces[i]));
      }
    }

    return new meshset_t(vertex_storage, meshes);
  }



  static void copyMeshFaces(const mesh::Mesh<3> *mesh,
                            size_t manifold_id,
                            const mesh::Vertex<3> *Vbase,
                            poly::Polyhedron *poly,
                            std::unordered_map<std::pair<size_t, size_t>, std::list<mesh::Edge<3> *> > &edges,
                            std::unordered_map<const mesh::Face<3> *, size_t> &face_map) {
    std::vector<const poly::Polyhedron::vertex_t *> vert_ptr;
    for (size_t f = 0; f < mesh->faces.size(); ++f) {
      mesh::Face<3> *src = mesh->faces[f];
      vert_ptr.clear();
      vert_ptr.reserve(src->nVertices());
      mesh::Edge<3> *e = src->edge;
      do {
        vert_ptr.push_back(&poly->vertices[e->vert - Vbase]);
        edges[std::make_pair(e->v1() - Vbase, e->v2() - Vbase)].push_back(e);
        e = e->next;
      } while (e != src->edge);
      
      face_map[src] = poly->faces.size();;
      
      poly->faces.push_back(poly::Polyhedron::face_t(vert_ptr));
      poly->faces.back().manifold_id = manifold_id;
      poly->faces.back().owner = poly;
    }
  }



  // construct a Polyhedron from a MeshSet
  poly::Polyhedron *polyhedronFromMesh(const mesh::MeshSet<3> *mesh, int manifold_id) {
    typedef poly::Polyhedron poly_t;
    typedef poly::Polyhedron::vertex_t vertex_t;
    typedef poly::Polyhedron::edge_t edge_t;
    typedef poly::Polyhedron::face_t face_t;

    poly::Polyhedron *poly = new poly::Polyhedron();
    const mesh::Vertex<3> *Vbase = &mesh->vertex_storage[0];

    poly->vertices.reserve(mesh->vertex_storage.size());
    for (size_t i = 0; i < mesh->vertex_storage.size(); ++i) {
      poly->vertices.push_back(vertex_t(mesh->vertex_storage[i].v));
      poly->vertices.back().owner = poly;
    }

    size_t n_faces = 0;
    if (manifold_id == -1) {
      poly->manifold_is_closed.resize(mesh->meshes.size());
      poly->manifold_is_negative.resize(mesh->meshes.size());
      for (size_t m = 0; m < mesh->meshes.size(); ++m) {
        n_faces += mesh->meshes[m]->faces.size();
        poly->manifold_is_closed[m] = mesh->meshes[m]->isClosed();
        poly->manifold_is_negative[m] = mesh->meshes[m]->isNegative();
      }
    } else {
      poly->manifold_is_closed.resize(1);
      poly->manifold_is_negative.resize(1);
      n_faces = mesh->meshes[manifold_id]->faces.size();
      poly->manifold_is_closed[manifold_id] = mesh->meshes[manifold_id]->isClosed();
      poly->manifold_is_negative[manifold_id] = mesh->meshes[manifold_id]->isNegative();
    }

    std::unordered_map<std::pair<size_t, size_t>, std::list<mesh::Edge<3> *> > edges;
    std::unordered_map<const mesh::Face<3> *, size_t> face_map;
    poly->faces.reserve(n_faces);

    if (manifold_id == -1) {
      for (size_t m = 0; m < mesh->meshes.size(); ++m) {
        copyMeshFaces(mesh->meshes[m], m, Vbase, poly, edges, face_map);
      }
    } else {
      copyMeshFaces(mesh->meshes[manifold_id], 0, Vbase, poly, edges, face_map);
    }

    size_t n_edges = 0;
    for (std::unordered_map<std::pair<size_t, size_t>, std::list<mesh::Edge<3> *> >::iterator i = edges.begin(); i != edges.end(); ++i) {
      if ((*i).first.first < (*i).first.second || edges.find(std::make_pair((*i).first.second, (*i).first.first)) == edges.end()) {
        n_edges++;
      }
    }

    poly->edges.reserve(n_edges);
    for (std::unordered_map<std::pair<size_t, size_t>, std::list<mesh::Edge<3> *> >::iterator i = edges.begin(); i != edges.end(); ++i) {
      if ((*i).first.first < (*i).first.second || edges.find(std::make_pair((*i).first.second, (*i).first.first)) == edges.end()) {
        poly->edges.push_back(edge_t(&poly->vertices[(*i).first.first],
                                     &poly->vertices[(*i).first.second],
                                     poly));
      }
    }

    poly->initVertexConnectivity();

    // build edge entries for face.
    for (size_t f = 0; f < poly->faces.size(); ++f) {
      face_t &face = poly->faces[f];
      size_t N = face.nVertices();
      for (size_t v = 0; v < N; ++v) {
        size_t v1i = poly->vertexToIndex_fast(face.vertex(v));
        size_t v2i = poly->vertexToIndex_fast(face.vertex((v+1)%N));
        std::vector<const edge_t *> found_edge;
        std::set_intersection(poly->connectivity.vertex_to_edge[v1i].begin(), poly->connectivity.vertex_to_edge[v1i].end(),
                              poly->connectivity.vertex_to_edge[v2i].begin(), poly->connectivity.vertex_to_edge[v2i].end(),
                              std::back_inserter(found_edge));
        CARVE_ASSERT(found_edge.size() == 1);
        face.edge(v) = found_edge[0];
      }
    }

    poly->connectivity.edge_to_face.resize(poly->edges.size());

    for (size_t i = 0; i < poly->edges.size(); ++i) {
      size_t v1i = poly->vertexToIndex_fast(poly->edges[i].v1);
      size_t v2i = poly->vertexToIndex_fast(poly->edges[i].v2);
      std::list<mesh::Edge<3> *> &efwd = edges[std::make_pair(v1i, v2i)];
      std::list<mesh::Edge<3> *> &erev = edges[std::make_pair(v1i, v2i)];

      for (std::list<mesh::Edge<3> *>::iterator j = efwd.begin(); j != efwd.end(); ++j) {
        mesh::Edge<3> *edge = *j;
        if (face_map.find(edge->face) != face_map.end()) {
          poly->connectivity.edge_to_face[i].push_back(&poly->faces[face_map[edge->face]]);
          if (edge->rev == NULL) {
            poly->connectivity.edge_to_face[i].push_back(NULL);
          } else {
            poly->connectivity.edge_to_face[i].push_back(&poly->faces[face_map[edge->rev->face]]);
          }
        }
      }
      for (std::list<mesh::Edge<3> *>::iterator j = erev.begin(); j != erev.end(); ++j) {
        mesh::Edge<3> *edge = *j;
        if (face_map.find(edge->face) != face_map.end()) {
          if (edge->rev == NULL) {
            poly->connectivity.edge_to_face[i].push_back(NULL);
            poly->connectivity.edge_to_face[i].push_back(&poly->faces[face_map[edge->face]]);
          }
        }
      }

    }

    poly->initSpatialIndex();

    // XXX: at this point, manifold_is_negative is not set up. This
    // info should be computed/stored in Mesh instances.

    return poly;
  }
}



