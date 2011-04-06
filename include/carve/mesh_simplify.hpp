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
#include <carve/mesh.hpp>
#include <carve/mesh_ops.hpp>
#include <carve/geom2d.hpp>
#include <carve/heap.hpp>

#include <fstream>
#include <string>
#include <utility>
#include <set>
#include <algorithm>
#include <vector>



namespace carve {
  namespace mesh {


    class MeshSimplifier {
      typedef carve::mesh::MeshSet<3> meshset_t;
      typedef carve::mesh::Mesh<3> mesh_t;
      typedef mesh_t::vertex_t vertex_t;
      typedef mesh_t::edge_t edge_t;
      typedef mesh_t::face_t face_t;



      struct EdgeInfo {
        edge_t *edge;
        double delta_v;

        double c[4];
        double l[2], t1[2], t2[2];
        size_t heap_idx;

        void update() {
          const vertex_t *v1 = edge->vert;
          const vertex_t *v2 = edge->next->vert;
          const vertex_t *v3 = edge->next->next->vert;
          const vertex_t *v4 = edge->rev ? edge->rev->next->next->vert : NULL;

          l[0] = (v1->v - v2->v).length();

          t1[0] = (v3->v - v1->v).length();
          t1[1] = (v3->v - v2->v).length();

          c[0] = std::max((t1[0] + t1[1]) / l[0] - 1.0, 0.0);

          if (v4) {
            l[1] = (v3->v - v4->v).length();
            t2[0] = (v4->v - v1->v).length();
            t2[1] = (v4->v - v2->v).length();
            c[1] = std::max((t2[0] + t2[1]) / l[0] - 1.0, 0.0);
            c[2] = std::max((t1[0] + t2[0]) / l[1] - 1.0, 0.0);
            c[3] = std::max((t1[1] + t2[1]) / l[1] - 1.0, 0.0);
            delta_v = carve::geom3d::tetrahedronVolume(v1->v, v2->v, v3->v, v4->v);
          } else {
            l[1] = 0.0;
            t2[0] = t2[1] = 0.0;
            c[1] = c[2] = c[3] = 0.0;
            delta_v = 0.0;
          }
        }

        EdgeInfo(edge_t *e) : edge(e) {
          update();
        }

        EdgeInfo() : edge(NULL) {
          delta_v = 0.0;
          c[0] = c[1] = c[2] = c[3] = 0.0;
          l[0] = l[1] = 0.0;
          t1[0] = t1[1] = 0.0;
          t2[0] = t2[1] = 0.0;
        }

        struct NotifyPos {
          void operator()(EdgeInfo *edge, size_t pos) const { edge->heap_idx = pos; }
          void operator()(EdgeInfo &edge, size_t pos) const { edge.heap_idx = pos; }
        };
      };



      struct FlippableBase {
        virtual bool canFlip(const EdgeInfo *) const =0;

        bool open(const EdgeInfo *e) const {
          return e->edge->rev == NULL;
        }

        bool flippable_DotProd(const EdgeInfo *e) const { 
          if (open(e)) return false;

          edge_t *edge = e->edge;

          const vertex_t *v1 = edge->vert;
          const vertex_t *v2 = edge->next->vert;
          const vertex_t *v3 = edge->next->next->vert;
          const vertex_t *v4 = edge->rev->next->next->vert;

          if (carve::geom::dot(carve::geom::cross(v3->v - v2->v, v1->v - v2->v),
                               carve::geom::cross(v4->v - v1->v, v2->v - v1->v)) < 0.95) return false;

          if (carve::geom::dot(carve::geom::cross(v3->v - v4->v, v1->v - v4->v),
                               carve::geom::cross(v4->v - v3->v, v2->v - v3->v)) < 0.95) return false;

          return true;
        }

        double score(const EdgeInfo *e) const {
          return std::min(e->c[2], e->c[3]) - std::min(e->c[0], e->c[1]);
        }

        struct Priority {
          const FlippableBase &flip;
          Priority(const FlippableBase &_flip) : flip(_flip) {}
          bool operator()(const EdgeInfo *a, const EdgeInfo *b) const { return flip.score(a) > flip.score(b); }
        };

        Priority priority() const {
          return Priority(*this);
        }
      };



      struct FlippableConservative : public FlippableBase {
        bool connectsAlmostCoplanarFaces(const EdgeInfo *e) const {
          // XXX: remove hard coded constants.
          if (e->c[0] < 1e-10 || e->c[1] < 1e-10) return true;
          return fabs(carve::geom::dot(e->edge->face->plane.N, e->edge->rev->face->plane.N) - 1.0) < 1e-10;
        }

        bool connectsExactlyCoplanarFaces(const EdgeInfo *e) const {
          edge_t *edge = e->edge;
          return
            carve::geom3d::orient3d(edge->vert->v,
                                    edge->next->vert->v,
                                    edge->next->next->vert->v,
                                    edge->rev->next->next->vert->v) == 0.0 &&
            carve::geom3d::orient3d(edge->rev->vert->v,
                                    edge->rev->next->vert->v,
                                    edge->rev->next->next->vert->v,
                                    edge->next->next->vert->v) == 0.0;
        }

        bool canFlip(const EdgeInfo *e) const {
          return !open(e) && score(e) > 0.0 && connectsExactlyCoplanarFaces(e) && flippable_DotProd(e);
        }
      };



      struct Flippable : public FlippableBase {
        double min_colinearity;
        double min_delta_v;

        Flippable(double _min_colinearity,
                  double _min_delta_v) :
            FlippableBase(),
            min_colinearity(_min_colinearity),
            min_delta_v(_min_delta_v) {
        }


        bool canFlip(const EdgeInfo *e) const {
          if (open(e)) return false;
          if (score(e) <= 0.0) return false;

          if (fabs(e->delta_v) > min_delta_v) return false;

          if (std::min(e->c[0], e->c[1]) > min_colinearity) return false;

          return flippable_DotProd(e);
        }
      };



      struct EdgeMerger {
        double min_edgelen;

        virtual bool canMerge(const EdgeInfo *e) const {
          return e->l[0] <= min_edgelen;
        }

        EdgeMerger(double _min_edgelen) : min_edgelen(_min_edgelen) {
        }

        double score(const EdgeInfo *e) const {
          return min_edgelen - e->l[0];
        }

        struct Priority {
          const EdgeMerger &merger;
          Priority(const EdgeMerger &_merger) : merger(_merger) {
          }
          bool operator()(const EdgeInfo *a, const EdgeInfo *b) const {
            // collapse edges in order from shortest to longest.
            return merger.score(a) < merger.score(b);
          }
        };

        Priority priority() const {
          return Priority(*this);
        }
      };



      typedef std::unordered_map<edge_t *, EdgeInfo *> edge_info_map_t;
      std::unordered_map<edge_t *, EdgeInfo *> edge_info;



      void initEdgeInfo(mesh_t *mesh) {
        for (size_t i = 0; i < mesh->faces.size(); ++i) {
          edge_t *e = mesh->faces[i]->edge;
          do {
            edge_info[e] = new EdgeInfo(e);
            e = e->next;
          } while (e != mesh->faces[i]->edge);
        }
      }



      void initEdgeInfo(meshset_t *meshset) {
        for (size_t m = 0; m < meshset->meshes.size(); ++m) {
          mesh_t *mesh = meshset->meshes[m];
          initEdgeInfo(mesh);
        }
      }



      void clearEdgeInfo() {
        for (edge_info_map_t::iterator i = edge_info.begin(); i != edge_info.end(); ++i) {
          delete (*i).second;
        }
      }



      void updateEdgeFlipHeap(std::vector<EdgeInfo *> &edge_heap,
                              edge_t *edge,
                              const FlippableBase &flipper) {
        std::unordered_map<edge_t *, EdgeInfo *>::const_iterator i = edge_info.find(edge);
        CARVE_ASSERT(i != edge_info.end());
        EdgeInfo *e = (*i).second;

        bool heap_pre = e->heap_idx != -1U;
        (*i).second->update();
        bool heap_post = edge->v1() < edge->v2() && flipper.canFlip(e);

        if (!heap_pre and heap_post) {
          edge_heap.push_back(e);
          carve::heap::push_heap(edge_heap.begin(),
                                 edge_heap.end(),
                                 flipper.priority(),
                                 EdgeInfo::NotifyPos());
        } else if (heap_pre && !heap_post) {
          CARVE_ASSERT(edge_heap[e->heap_idx] == e);
          carve::heap::remove_heap(edge_heap.begin(),
                                   edge_heap.end(),
                                   edge_heap.begin() + e->heap_idx,
                                   flipper.priority(),
                                   EdgeInfo::NotifyPos());
          CARVE_ASSERT(edge_heap.back() == e);
          edge_heap.pop_back();
          e->heap_idx = -1U;
        } else if (heap_pre && heap_post) {
          CARVE_ASSERT(edge_heap[e->heap_idx] == e);
          carve::heap::adjust_heap(edge_heap.begin(),
                                   edge_heap.end(),
                                   edge_heap.begin() + e->heap_idx,
                                   flipper.priority(),
                                   EdgeInfo::NotifyPos());
          CARVE_ASSERT(edge_heap[e->heap_idx] == e);
        }
      }



      size_t flipEdges(meshset_t *meshset,
                       const FlippableBase &flipper) {
        size_t n_mods = 0;

        std::vector<EdgeInfo *> edge_heap;

        edge_heap.reserve(edge_info.size());

        for (edge_info_map_t::iterator i = edge_info.begin();
             i != edge_info.end();
             ++i) {
          EdgeInfo *e = (*i).second;

          if (e->edge->v1() < e->edge->v2() && flipper.canFlip(e)) {
            edge_heap.push_back(e);
          } else {
            e->heap_idx = -1U;
          }
        }

        carve::heap::make_heap(edge_heap.begin(),
                               edge_heap.end(),
                               flipper.priority(),
                               EdgeInfo::NotifyPos());

        while (edge_heap.size()) {
          carve::heap::pop_heap(edge_heap.begin(),
                                edge_heap.end(),
                                flipper.priority(),
                                EdgeInfo::NotifyPos());
          EdgeInfo *e = edge_heap.back();
          edge_heap.pop_back();
          e->heap_idx = -1U;

          n_mods++;
          CARVE_ASSERT(flipper.canFlip(e));
          carve::mesh::flipTriEdge(e->edge);

          updateEdgeFlipHeap(edge_heap, e->edge, flipper);
          updateEdgeFlipHeap(edge_heap, e->edge->rev, flipper);

          CARVE_ASSERT(!flipper.canFlip(e));

          updateEdgeFlipHeap(edge_heap, e->edge->next, flipper);
          updateEdgeFlipHeap(edge_heap, e->edge->next->next, flipper);
          updateEdgeFlipHeap(edge_heap, e->edge->rev->next, flipper);
          updateEdgeFlipHeap(edge_heap, e->edge->rev->next->next, flipper);
          updateEdgeFlipHeap(edge_heap, e->edge->next->rev, flipper);
          updateEdgeFlipHeap(edge_heap, e->edge->next->next->rev, flipper);
          updateEdgeFlipHeap(edge_heap, e->edge->rev->next->rev, flipper);
          updateEdgeFlipHeap(edge_heap, e->edge->rev->next->next->rev, flipper);
        }

        return n_mods;
      }



      void removeFromEdgeMergeHeap(std::vector<EdgeInfo *> &edge_heap,
                                   EdgeInfo *edge,
                                   const EdgeMerger &merger) {
        if (edge->heap_idx != -1U) {
          CARVE_ASSERT(edge_heap[edge->heap_idx] == edge);
          carve::heap::remove_heap(edge_heap.begin(),
                                   edge_heap.end(),
                                   edge_heap.begin() + edge->heap_idx,
                                   merger.priority(),
                                   EdgeInfo::NotifyPos());
          CARVE_ASSERT(edge_heap.back() == edge);
          edge_heap.pop_back();
          edge->heap_idx = -1U;
        }
      }

      void updateEdgeMergeHeap(std::vector<EdgeInfo *> &edge_heap,
                               EdgeInfo *edge,
                               const EdgeMerger &merger) {
        bool heap_pre = edge->heap_idx != -1U;
        edge->update();
        bool heap_post = merger.canMerge(edge);

        if (!heap_pre and heap_post) {
          edge_heap.push_back(edge);
          carve::heap::push_heap(edge_heap.begin(),
                                 edge_heap.end(),
                                 merger.priority(),
                                 EdgeInfo::NotifyPos());
        } else if (heap_pre && !heap_post) {
          CARVE_ASSERT(edge_heap[edge->heap_idx] == edge);
          carve::heap::remove_heap(edge_heap.begin(),
                                   edge_heap.end(),
                                   edge_heap.begin() + edge->heap_idx,
                                   merger.priority(),
                                   EdgeInfo::NotifyPos());
          CARVE_ASSERT(edge_heap.back() == edge);
          edge_heap.pop_back();
          edge->heap_idx = -1U;
        } else if (heap_pre && heap_post) {
          CARVE_ASSERT(edge_heap[edge->heap_idx] == edge);
          carve::heap::adjust_heap(edge_heap.begin(),
                                   edge_heap.end(),
                                   edge_heap.begin() + edge->heap_idx,
                                   merger.priority(),
                                   EdgeInfo::NotifyPos());
          CARVE_ASSERT(edge_heap[edge->heap_idx] == edge);
        }
      }



      size_t mergeEdges(meshset_t *mesh,
                        const EdgeMerger &merger) {
        size_t n_mods = 0;

        std::vector<EdgeInfo *> edge_heap;
        std::unordered_map<vertex_t *, std::set<EdgeInfo *> > vert_to_edges;

        edge_heap.reserve(edge_info.size());

        for (edge_info_map_t::iterator i = edge_info.begin();
             i != edge_info.end();
             ++i) {
          EdgeInfo *e = (*i).second;

          vert_to_edges[e->edge->v1()].insert(e);
          vert_to_edges[e->edge->v2()].insert(e);

          if (merger.canMerge(e)) {
            edge_heap.push_back(e);
          } else {
            e->heap_idx = -1U;
          }
        }

        carve::heap::make_heap(edge_heap.begin(),
                               edge_heap.end(),
                               merger.priority(),
                               EdgeInfo::NotifyPos());

        while (edge_heap.size()) {
          carve::heap::pop_heap(edge_heap.begin(),
                                edge_heap.end(),
                                merger.priority(),
                                EdgeInfo::NotifyPos());
          EdgeInfo *e = edge_heap.back();
          edge_heap.pop_back();
          e->heap_idx = -1U;

          edge_t *edge = e->edge;
          vertex_t *v1 = edge->v1();
          vertex_t *v2 = edge->v2();

          std::vector<EdgeInfo *> edges_to_merge;
          std::vector<EdgeInfo *> v1_incident;
          std::vector<EdgeInfo *> v2_incident;

          std::set_intersection(vert_to_edges[v1].begin(), vert_to_edges[v1].end(),
                                vert_to_edges[v2].begin(), vert_to_edges[v2].end(),
                                std::back_inserter(edges_to_merge));

          CARVE_ASSERT(edges_to_merge.size() > 0);

          std::set_difference(vert_to_edges[v1].begin(), vert_to_edges[v1].end(),
                              edges_to_merge.begin(), edges_to_merge.end(),
                              std::back_inserter(v1_incident));
          std::set_difference(vert_to_edges[v2].begin(), vert_to_edges[v2].end(),
                              edges_to_merge.begin(), edges_to_merge.end(),
                              std::back_inserter(v2_incident));

          ++n_mods;

          double frac = 0.5; // compute this based upon v1_incident and v2_incident?
          v2->v = frac * v1->v + (1 - frac) * v2->v;

          for (size_t i = 0; i < v1_incident.size(); ++i) {
            if (v1_incident[i]->edge->vert == v1) {
              v1_incident[i]->edge->vert = v2;
            }
          }

          for (size_t i = 0; i < v1_incident.size(); ++i) {
            updateEdgeMergeHeap(edge_heap, v1_incident[i], merger);
          }

          for (size_t i = 0; i < v2_incident.size(); ++i) {
            updateEdgeMergeHeap(edge_heap, v2_incident[i], merger);
          }

          vert_to_edges[v2].insert(vert_to_edges[v1].begin(), vert_to_edges[v1].end());
          vert_to_edges.erase(v1);

          for (size_t i = 0; i < edges_to_merge.size(); ++i) {
            EdgeInfo *e = edges_to_merge[i];

            removeFromEdgeMergeHeap(edge_heap, e, merger);
            edge_info.erase(e->edge);

            vert_to_edges[v1].erase(e);
            vert_to_edges[v2].erase(e);

            face_t *f1 = e->edge->face;

            e->edge->removeHalfEdge();

            if (f1->n_edges == 2) {
              edge_t *e1 = f1->edge;
              edge_t *e2 = f1->edge->next;
              e1->rev->rev = e2->rev;
              e2->rev->rev = e1->rev;
              EdgeInfo *e1i = edge_info[e1];
              EdgeInfo *e2i = edge_info[e2];
              CARVE_ASSERT(e1i != NULL);
              CARVE_ASSERT(e2i != NULL);
              vert_to_edges[e1->v1()].erase(e1i);
              vert_to_edges[e1->v2()].erase(e1i);
              vert_to_edges[e2->v1()].erase(e2i);
              vert_to_edges[e2->v2()].erase(e2i);
              removeFromEdgeMergeHeap(edge_heap, e1i, merger);
              removeFromEdgeMergeHeap(edge_heap, e2i, merger);
              delete e1i;
              delete e2i;
              edge_info.erase(e1);
              edge_info.erase(e2);
              f1->clearEdges();
            }
          }
        }

        return n_mods;
      }



      void removeRemnantFaces(mesh_t *mesh) {
        size_t n = 0;
        for (size_t i = 0; i < mesh->faces.size(); ++i) {
          if (mesh->faces[i]->nEdges() == 0) {
            delete mesh->faces[i];
          } else {
            mesh->faces[n++] = mesh->faces[i];
          }
        }
        mesh->faces.resize(n);
      }



      void removeRemnantFaces(meshset_t *mesh) {
        for (size_t i = 0; i < mesh->meshes.size(); ++i) {
          removeRemnantFaces(mesh->meshes[i]);
        }
      }



    public:
      size_t improveMesh_conservative(meshset_t *meshset) {
        initEdgeInfo(meshset);
        size_t modifications = flipEdges(meshset, FlippableConservative());
        clearEdgeInfo();
        return modifications;
      }



      size_t improveMesh(meshset_t *meshset,
                         double min_colinearity,
                         double min_delta_v) {
        initEdgeInfo(meshset);
        size_t modifications = flipEdges(meshset, Flippable(min_colinearity, min_delta_v));
        clearEdgeInfo();
        return modifications;
      }



      size_t eliminateShortEdges(meshset_t *meshset,
                                 double min_length) {
        initEdgeInfo(meshset);
        size_t modifications = mergeEdges(meshset, EdgeMerger(min_length));
        removeRemnantFaces(meshset);
        clearEdgeInfo();
        return modifications;
      }



      size_t simplify(meshset_t *meshset,
                      double min_colinearity,
                      double min_delta_v,
                      double min_length) {
        size_t modifications = 0;
        size_t n_flip, n_merge;

        initEdgeInfo(meshset);

        std::cerr << "initial merge" << std::endl;
        modifications = mergeEdges(meshset, EdgeMerger(0.0));

        do {
          n_flip = n_merge = 0;
          std::cerr << "flip conservative" << std::endl;
          n_flip = flipEdges(meshset, FlippableConservative());
          // std::cerr << "flip" << std::endl;
          // n_flip += flipEdges(meshset, Flippable(min_colinearity, min_delta_v));
          std::cerr << "merge" << std::endl;
          n_merge = mergeEdges(meshset, EdgeMerger(min_length));
          modifications += n_flip + n_merge;
          std::cerr << "stats:" << n_flip << " " << n_merge << std::endl;
        } while (n_flip || n_merge);

        removeRemnantFaces(meshset);
        clearEdgeInfo();
        return modifications;
      }
    };
  }
}
