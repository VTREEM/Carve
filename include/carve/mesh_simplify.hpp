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
        double min_dp;

        FlippableBase(double _min_dp = 0.0) : min_dp(_min_dp) {
        }

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
                               carve::geom::cross(v4->v - v1->v, v2->v - v1->v)) < min_dp) return false;

          if (carve::geom::dot(carve::geom::cross(v3->v - v4->v, v1->v - v4->v),
                               carve::geom::cross(v4->v - v3->v, v2->v - v3->v)) < min_dp) return false;

          return true;
        }

        double score(const EdgeInfo *e) const {
          return std::min(e->c[2], e->c[3]) - std::min(e->c[0], e->c[1]);
        }

        class Priority {
          Priority &operator=(const Priority &);
          const FlippableBase &flip;

        public:
          Priority(const FlippableBase &_flip) : flip(_flip) {}
          bool operator()(const EdgeInfo *a, const EdgeInfo *b) const { return flip.score(a) > flip.score(b); }
        };

        Priority priority() const {
          return Priority(*this);
        }
      };



      struct FlippableConservative : public FlippableBase {
        FlippableConservative() : FlippableBase(0.0) {
        }

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
                  double _min_delta_v,
                  double _min_normal_angle) :
            FlippableBase(cos(_min_normal_angle)),
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

        class Priority {
          Priority &operator=(const Priority &);

        public:
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

        bool heap_pre = e->heap_idx != ~0U;
        (*i).second->update();
        bool heap_post = edge->v1() < edge->v2() && flipper.canFlip(e);

        if (!heap_pre && heap_post) {
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
          e->heap_idx = ~0U;
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



      size_t flipEdges(meshset_t * /* meshset */,
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
            e->heap_idx = ~0U;
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
          e->heap_idx = ~0U;

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
        if (edge->heap_idx != ~0U) {
          CARVE_ASSERT(edge_heap[edge->heap_idx] == edge);
          carve::heap::remove_heap(edge_heap.begin(),
                                   edge_heap.end(),
                                   edge_heap.begin() + edge->heap_idx,
                                   merger.priority(),
                                   EdgeInfo::NotifyPos());
          CARVE_ASSERT(edge_heap.back() == edge);
          edge_heap.pop_back();
          edge->heap_idx = ~0U;
        }
      }

      void updateEdgeMergeHeap(std::vector<EdgeInfo *> &edge_heap,
                               EdgeInfo *edge,
                               const EdgeMerger &merger) {
        bool heap_pre = edge->heap_idx != ~0U;
        edge->update();
        bool heap_post = merger.canMerge(edge);

        if (!heap_pre && heap_post) {
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
          edge->heap_idx = ~0U;
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



      // collapse edges edges based upon the predicate implemented by EdgeMerger.
      size_t collapseEdges(meshset_t * /* mesh */,
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
            e->heap_idx = ~0U;
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
          e->heap_idx = ~0U;

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



      size_t mergeCoplanarFaces(mesh_t *mesh, double min_normal_angle) {
        std::unordered_set<edge_t *> coplanar_face_edges;
        double min_dp = cos(min_normal_angle);
        size_t n_merge = 0;

        for (size_t i = 0; i < mesh->closed_edges.size(); ++i) {
          edge_t *e = mesh->closed_edges[i];
          face_t *f1 = e->face;
          face_t *f2 = e->rev->face;

          if (carve::geom::dot(f1->plane.N, f2->plane.N) < min_dp) {
            continue;
          }

          coplanar_face_edges.insert(std::min(e, e->rev));
        }

        while (coplanar_face_edges.size()) {
          edge_t *edge = *coplanar_face_edges.begin();
          if (edge->face == edge->rev->face) {
            coplanar_face_edges.erase(edge);
            continue;
          }

          edge_t *removed = edge->mergeFaces();
          if (removed == NULL) {
            coplanar_face_edges.erase(edge);
            ++n_merge;
          } else {
            edge_t *e = removed;
            do {
              edge_t *n = e->next;
              coplanar_face_edges.erase(std::min(e, e->rev));
              delete e->rev;
              delete e;
              e = n;
            } while (e != removed);
          }
        }
        return n_merge;
      }



      uint8_t affected_axes(const face_t *face) {
        uint8_t r = 0;
        if (fabs(carve::geom::dot(face->plane.N, carve::geom::VECTOR(1,0,0))) > 0.001) r |= 1;
        if (fabs(carve::geom::dot(face->plane.N, carve::geom::VECTOR(0,1,0))) > 0.001) r |= 2;
        if (fabs(carve::geom::dot(face->plane.N, carve::geom::VECTOR(0,0,1))) > 0.001) r |= 4;
        return r;
      }



      double median(std::vector<double> &v) {
        if (v.size() & 1) {
          size_t N = v.size() / 2 + 1;
          std::nth_element(v.begin(), v.begin() + N, v.end());
          return v[N];
        } else {
          size_t N = v.size() / 2;
          std::nth_element(v.begin(), v.begin() + N, v.end());
          return (v[N] + *std::min_element(v.begin() + N + 1, v.end())) / 2.0;
        }
      }



      double harmonicmean(const std::vector<double> &v) {
        double m = 0.0;
        for (size_t i = 0; i < v.size(); ++i) {
          m *= v[i];
        }
        return pow(m, 1.0 / v.size());
      }



      double mean(const std::vector<double> &v) {
        double m = 0.0;
        for (size_t i = 0; i < v.size(); ++i) {
          m += v[i];
        }
        return m / v.size();
      }



      template<typename iter_t>
      void snapFaces(iter_t begin, iter_t end, double grid, int axis) {
        std::set<vertex_t *> vertices;
        for (iter_t i = begin; i != end; ++i) {
          face_t *face = *i;
          edge_t *edge = face->edge;
          do {
            vertices.insert(edge->vert);
            edge = edge->next;
          } while (edge != face->edge);
        }

        std::vector<double> pos;
        pos.reserve(vertices.size());
        for (std::set<vertex_t *>::iterator i = vertices.begin(); i != vertices.end(); ++i) {
          pos.push_back((*i)->v.v[axis]);
        }

        double med = median(pos);

        double snap_pos = med;
        if (grid) snap_pos = round(snap_pos / grid) * grid;

        for (std::set<vertex_t *>::iterator i = vertices.begin(); i != vertices.end(); ++i) {
          (*i)->v.v[axis] = snap_pos;
        }

        for (iter_t i = begin; i != end; ++i) {
          face_t *face = *i;
          face->recalc();
          edge_t *edge = face->edge;
          do {
            if (edge->rev && edge->rev->face) edge->rev->face->recalc();
            edge = edge->next;
          } while (edge != face->edge);
        }
      }

      carve::geom::plane<3> quantizePlane(const face_t *face,
                                          int angle_xy_quantization,
                                          int angle_z_quantization) {
        if (!angle_xy_quantization && !angle_z_quantization) {
          return face->plane;
        }
        carve::geom::vector<3> normal = face->plane.N;

        if (angle_z_quantization) {
          if (normal.x || normal.y) {
            double a = asin(std::min(std::max(normal.z, 0.0), 1.0));
            a = round(a * angle_z_quantization / (M_PI * 2)) * (M_PI * 2) / angle_z_quantization;
            normal.z = sin(a);
            double s = sqrt((1 - normal.z * normal.z) / (normal.x * normal.x + normal.y * normal.y));
            normal.x = normal.x * s;
            normal.y = normal.y * s;
          }
        }
        if (angle_xy_quantization) {
          if (normal.x || normal.y) {
            double a = atan2(normal.y, normal.x);
            a = round(a * angle_xy_quantization / (M_PI * 2)) * (M_PI * 2) / angle_xy_quantization;
            double s = sqrt(1 - normal.z * normal.z);
            s = std::min(std::max(s, 0.0), 1.0);
            normal.x = cos(a) * s;
            normal.y = sin(a) * s;
          }
        }

        std::cerr << "normal = " << normal << std::endl;

        std::vector<double> d_vec;
        d_vec.reserve(face->nVertices());
        edge_t *e = face->edge;
        do {
          d_vec.push_back(-carve::geom::dot(normal, e->vert->v));
          e = e->next;
        } while (e != face->edge);

        return carve::geom::plane<3>(normal, mean(d_vec));
      }



      double summedError(const carve::geom::vector<3> &vert, const std::list<carve::geom::plane<3> > &planes) {
        double d = 0;
        for (std::list<carve::geom::plane<3> >::const_iterator i = planes.begin(); i != planes.end(); ++i) {
          d += fabs(carve::geom::distance2(*i, vert));
        }
        return d;
      }



      double minimize(carve::geom::vector<3> &vert, const std::list<carve::geom::plane<3> > &planes, int axis) {
        double num = 0.0;
        double den = 0.0;
        int a1 = (axis + 1) % 3;
        int a2 = (axis + 2) % 3;
        for (std::list<carve::geom::plane<3> >::const_iterator i = planes.begin(); i != planes.end(); ++i) {
          const carve::geom::vector<3> &N = (*i).N;
          const double d = (*i).d;
          den += N.v[axis] * N.v[axis];
          num -= N.v[axis] * (N.v[a1] * vert.v[a1] + N.v[a2] * vert.v[a2] + d);
        }
        if (fabs(den) < 1e-5) return vert.v[axis];
        return num / den;
      }



      size_t cleanFaceEdges(mesh_t *mesh) {
        size_t n_removed = 0;
        for (size_t i = 0; i < mesh->faces.size(); ++i) {
          face_t *face = mesh->faces[i];
          edge_t *start = face->edge;
          edge_t *edge = start;
          do {
            if (edge->next == edge->rev || edge->prev == edge->rev) {
              edge = edge->removeEdge();
              ++n_removed;
              start = edge->prev;
            } else {
              edge = edge->next;
            }
          } while (edge != start);
        }
        return n_removed;
      }



      size_t cleanFaceEdges(meshset_t *mesh) {
        size_t n_removed = 0;
        for (size_t i = 0; i < mesh->meshes.size(); ++i) {
          n_removed += cleanFaceEdges(mesh->meshes[i]);
        }
        return n_removed;
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
      // Merge adjacent coplanar faces (where coplanar is determined
      // by dot-product >= cos(min_normal_angle)).
      size_t mergeCoplanarFaces(meshset_t *meshset, double min_normal_angle) {
        size_t n_removed = 0;
        for (size_t i = 0; i < meshset->meshes.size(); ++i) {
          n_removed += mergeCoplanarFaces(meshset->meshes[i], min_normal_angle);
          removeRemnantFaces(meshset->meshes[i]);
          cleanFaceEdges(meshset->meshes[i]);
          meshset->meshes[i]->cacheEdges();
        }
        return n_removed;
      }

      size_t improveMesh_conservative(meshset_t *meshset) {
        initEdgeInfo(meshset);
        size_t modifications = flipEdges(meshset, FlippableConservative());
        clearEdgeInfo();
        return modifications;
      }



      size_t improveMesh(meshset_t *meshset,
                         double min_colinearity,
                         double min_delta_v,
                         double min_normal_angle) {
        initEdgeInfo(meshset);
        size_t modifications = flipEdges(meshset, Flippable(min_colinearity, min_delta_v, min_normal_angle));
        clearEdgeInfo();
        return modifications;
      }



      size_t eliminateShortEdges(meshset_t *meshset,
                                 double min_length) {
        initEdgeInfo(meshset);
        size_t modifications = collapseEdges(meshset, EdgeMerger(min_length));
        removeRemnantFaces(meshset);
        clearEdgeInfo();
        return modifications;
      }



      // Snap vertices to grid, aligning almost flat axis-aligned
      // faces to the axis, and flattening other faces as much as is
      // possible. Passing a number less than DBL_MIN_EXPONENT (-1021)
      // turns off snapping to grid (but face alignment is still
      // performed).
      void snap(meshset_t *meshset,
                int log2_grid,
                int angle_xy_quantization = 0,
                int angle_z_quantization = 0) {
        double grid = 0.0;
        if (log2_grid >= std::numeric_limits<double>::min_exponent) grid = pow(2.0, (double)log2_grid);

        typedef std::unordered_map<face_t *, uint8_t> axis_influence_map_t;
        axis_influence_map_t axis_influence;

        typedef std::unordered_map<face_t *, std::set<face_t *> > interaction_graph_t;
        interaction_graph_t interacting_faces;

        for (size_t m = 0; m < meshset->meshes.size(); ++m) {
          mesh_t *mesh = meshset->meshes[m];
          for (size_t f = 0; f < mesh->faces.size(); ++f) {
            face_t *face = mesh->faces[f];
            axis_influence[face] = affected_axes(face);
          }
        }

        std::map<vertex_t *, std::list<carve::geom::plane<3> > > non_axis_vertices;
        std::unordered_map<vertex_t *, uint8_t> vertex_constraints;

        for (axis_influence_map_t::iterator i = axis_influence.begin(); i != axis_influence.end(); ++i) {
          face_t *face = (*i).first;
          uint8_t face_axes = (*i).second;
          edge_t *edge = face->edge;
          if (face_axes != 1 && face_axes != 2 && face_axes != 4) {
            do {
              non_axis_vertices[edge->vert].push_back(quantizePlane(face,
                                                                    angle_xy_quantization,
                                                                    angle_z_quantization));
              edge = edge->next;
            } while (edge != face->edge);
          } else {
            interacting_faces[face].insert(face);
            do {
              vertex_constraints[edge->vert] |= face_axes;

              if (edge->rev && edge->rev->face) {
                face_t *face2 = edge->rev->face;
                uint8_t face2_axes = axis_influence[face2];
                if (face2_axes == face_axes) {
                  interacting_faces[face].insert(face2);
                }
              }
              edge = edge->next;
            } while (edge != face->edge);
          }
        }

        while (interacting_faces.size()) {
          std::set<face_t *> face_set;
          uint8_t axes = 0;

          std::set<face_t *> open;
          open.insert((*interacting_faces.begin()).first);
          while (open.size()) {
            face_t *curr = *open.begin();
            open.erase(open.begin());
            face_set.insert(curr);
            axes |= axis_influence[curr];
            for (interaction_graph_t::data_type::iterator i = interacting_faces[curr].begin(), e = interacting_faces[curr].end(); i != e; ++i) {
              face_t *f = *i;
              if (face_set.find(f) != face_set.end()) continue;
              open.insert(f);
            }
          }

          switch (axes) {
          case 1: snapFaces(face_set.begin(), face_set.end(), grid, 0); break;
          case 2: snapFaces(face_set.begin(), face_set.end(), grid, 1); break;
          case 4: snapFaces(face_set.begin(), face_set.end(), grid, 2); break;
          default: CARVE_FAIL("should not be reached");
          }

          for (std::set<face_t *>::iterator i = face_set.begin(); i != face_set.end(); ++i) {
            interacting_faces.erase((*i));
          }
        }

        for (std::map<vertex_t *, std::list<carve::geom::plane<3> > >::iterator i = non_axis_vertices.begin(); i != non_axis_vertices.end(); ++i) {
          vertex_t *vert = (*i).first;
          std::list<carve::geom::plane<3> > &planes = (*i).second;
          uint8_t constraint = vertex_constraints[vert];

          if (constraint == 7) continue;

          double d = summedError(vert->v, planes);
          for (size_t N = 0; ; N = (N+1) % 3) {
            if (constraint & (1 << N)) continue;
            vert->v[N] = minimize(vert->v, planes, N);
            double d_next = summedError(vert->v, planes);
            if (d - d_next < 1e-20) break;
            d = d_next;
          }

          if (grid) {
            carve::geom::vector<3> v_best = vert->v;
            double d_best = 0.0;
            
            for (size_t axes = 0; axes < 8; ++axes) {
              carve::geom::vector<3> v = vert->v;
              for (size_t N = 0; N < 3; ++N) {
                if (constraint & (1 << N)) continue;
                if (axes & (1<<N)) {
                  v.v[N] = ceil(v.v[N] / grid) * grid;
                } else {
                  v.v[N] = floor(v.v[N] / grid) * grid;
                }
              }
              double d = summedError(v, planes);
              if (axes == 0 || d < d_best) {
                v_best = v;
                d_best = d;
              }
            }

            vert->v = v_best;
          }
        }
      }



      size_t simplify(meshset_t *meshset,
                      double min_colinearity,
                      double min_delta_v,
                      double min_normal_angle,
                      double min_length) {
        size_t modifications = 0;
        size_t n_flip, n_merge;

        initEdgeInfo(meshset);

        std::cerr << "initial merge" << std::endl;
        modifications = collapseEdges(meshset, EdgeMerger(0.0));

        do {
          n_flip = n_merge = 0;
          std::cerr << "flip conservative" << std::endl;
          n_flip = flipEdges(meshset, FlippableConservative());
          std::cerr << "flip" << std::endl;
          n_flip += flipEdges(meshset, Flippable(min_colinearity, min_delta_v, min_normal_angle));
          std::cerr << "merge" << std::endl;
          n_merge = collapseEdges(meshset, EdgeMerger(min_length));
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
