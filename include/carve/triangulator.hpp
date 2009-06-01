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

#include <list>
#include <vector>
#include <algorithm>

#include <carve/carve.hpp>

#include <carve/geom2d.hpp>

namespace carve {
  namespace triangulate {

    /**
     * \class order_h_loops
     * \brief Provides an ordering of hole loops based upon a single projected axis.
     *
     * @tparam project_t A functor which converts vertices to a 2d projection.
     * @tparam hole_t A collection of vertices.
     */
    template<typename project_t, typename hole_t>
    class order_h_loops {
      const project_t &project;
      int axis;
    public:

      /** 
       * 
       * @param _project The projection functor.
       * @param _axis The axis of the 2d projection upon which hole loops are ordered.
       */
      order_h_loops(const project_t &_project, int _axis) : project(_project), axis(_axis) { }

      bool operator()(const std::pair<const hole_t *, typename hole_t::const_iterator> &a,
                      const std::pair<const hole_t *, typename hole_t::const_iterator> &b) const {
        return project(*(a.second)).v[axis] < project(*(b.second)).v[axis];
      }
    };



    /**
     * \class heap_ordering
     * \brief Provides an ordering of vertex indicies in a polygon loop according to proximity to a vertex.
     *
     * @tparam project_t A functor which converts vertices to a 2d projection.
     * @tparam vert_t A vertex type.
     */
    template<typename project_t, typename vert_t>
    class heap_ordering {
      const project_t &project;
      const std::vector<vert_t> &loop;
      const carve::geom2d::P2 p;

    public:
      /** 
       * 
       * @param _project A functor which converts vertices to a 2d projection.
       * @param _loop The polygon loop which indices address.
       * @param _vert The vertex from which distance is measured.
       * 
       */
      heap_ordering(const project_t &_project, const std::vector<vert_t> &_loop, const vert_t &_vert) : project(_project), loop(_loop), p(_project(_vert)) { }

      bool operator()(size_t a, size_t b) const {
        return carve::geom::distance2(p, project(loop[a])) > carve::geom::distance2(p, project(loop[b]));
      }
    };



    /** 
     * \brief Given a polygon loop and a hole loop, and attachment points, insert the hole loop vertices into the polygon loop.
     * 
     * @param[in,out] f_loop The polygon loop to incorporate the hole into.
     * @param f_loop_attach[in] The index of the vertex of the polygon
     *                          loop that the hole is to be attached
     *                          to.
     * @param hole_attach[in] A pair consisting of a pointer to a hole
     *                        container and an iterator into that
     *                        container reflecting the point of
     *                        attachment of the hole.
     */
    template<typename vert_t, typename hole_t>
    void patchHoleIntoPolygon(std::vector<vert_t> &f_loop,
                           unsigned f_loop_attach,
                           const std::pair<const hole_t *, typename hole_t::const_iterator> &hole_attach) {
      // join the vertex curr of the polygon loop to the hole at
      // h_loop_connect
      f_loop.insert(f_loop.begin() + f_loop_attach + 1, hole_attach.first->size() + 2, vert_t());
      typename std::vector<vert_t>::iterator f = f_loop.begin() + f_loop_attach;

      typename hole_t::const_iterator h = hole_attach.second;

      while (h != hole_attach.first->end()) {
        *++f = *h++;
      }

      h = hole_attach.first->begin();
      typename hole_t::const_iterator he = hole_attach.second; ++he;
      while (h != he) {
        *++f = *h++;
      }

      *++f = f_loop[f_loop_attach];
    }



    /** 
     * \brief Merge a set of holes into a polygon.
     *
     * Take a polygon loop and a collection of hole loops, and patch
     * the hole loops into the polygon loop, returning a vector of
     * vertices from the polygon and holes which describe a new
     * polygon boundary with no holes, through the addition of edges
     * joining the polygon loop to the holes.
     * 
     * This may be applied to arbitrary vertex data (generally
     * carve::geom3d::Vertex pointers), but a projection function must
     * be supplied to convert vertices to coordinates in 2-space, in
     * which the work is performed.
     *
     * @tparam project_t A functor which converts vertices to a 2d
     *                   projection.
     * @tparam polygon_container_t A container type that represents an
     *                          ordered set of vertices.
     * @tparam hole_container_t A container type that represents a
     *                          possibly unordered set of ordered sets
     *                          of vertices.
     * @param project The projection functor.
     * @param f_loop The polygon loop into which holes are to be
     *               incorporated.
     * @param h_loops The set of hole loops to be incorporated.
     * 
     * @return A vector of vertices (of type polygon_container_t::value_type).
     */
    template<typename project_t, typename polygon_container_t, typename hole_container_t>
    static std::vector<typename polygon_container_t::value_type>
    incorporateHolesIntoPolygon(const project_t &project,
                                const polygon_container_t &f_loop,
                                const hole_container_t &h_loops) {
      typedef typename polygon_container_t::value_type vert_t;
      typedef typename hole_container_t::value_type hole_t;
      typedef typename polygon_container_t::const_iterator polygon_iter;
      typedef typename hole_container_t::const_iterator hole_iter;
      typedef typename hole_t::const_iterator hole_vert_iter;

      size_t N = f_loop.size();

      // work out how much space to reserve for the patched in holes.
      for (hole_iter i = h_loops.begin(); i != h_loops.end(); ++i) {
        N += 2 + (*i).size();
      }
    
      // this is the vector that we will build the result in.
      std::vector<vert_t> current_f_loop;
      current_f_loop.reserve(N);

      std::vector<size_t> f_loop_heap;
      f_loop_heap.reserve(N);

      for (unsigned i = 0; i < f_loop.size(); ++i) {
        current_f_loop.push_back(f_loop[i]);
      }

      std::vector<std::pair<const hole_t *, hole_vert_iter> > h_loop_min_vertex;

      h_loop_min_vertex.reserve(h_loops.size());

      // find the major axis for the holes - this is the axis that we
      // will sort on for finding vertices on the polygon to join
      // holes up to.
      //
      // it might also be nice to also look for whether it is better
      // to sort ascending or descending.
      // 
      // another trick that could be used is to modify the projection
      // by 90 degree rotations or flipping about an axis. just as
      // long as we keep the carve::geom3d::Vector pointers for the
      // real data in sync, everything should be ok. then we wouldn't
      // need to accomodate axes or sort order in the main loop.

      // find the bounding box of all the holes.
      bool first = true;
      double min_x, min_y, max_x, max_y;
      for (hole_iter i = h_loops.begin(); i != h_loops.end(); ++i) {
        const hole_t &hole(*i);
        for (hole_vert_iter j = hole.begin(); j != hole.end(); ++j) {
          carve::geom2d::P2 curr = project(*j);
          if (first) {
            min_x = max_x = curr.x;
            min_y = max_y = curr.y;
            first = false;
          } else {
            min_x = std::min(min_x, curr.x);
            min_y = std::min(min_y, curr.y);
            max_x = std::max(max_x, curr.x);
            max_y = std::max(max_y, curr.y);
          }
        }
      }

      // choose the axis for which the bbox is largest.
      int axis = (max_x - min_x) > (max_y - min_y) ? 0 : 1;

      // for each hole, find the minimum vertex in the chosen axis.
      for (hole_iter i = h_loops.begin(); i != h_loops.end(); ++i) {
        const hole_t &hole = *i;
        carve::geom2d::P2 best, curr;
        hole_vert_iter best_i, curr_i;
        best_i = curr_i = hole.begin();
        best = project(*best_i);
        for (++curr_i; curr_i != hole.end(); ++curr_i) {
          curr = project(*curr_i);
          if (curr.v[axis] < best.v[axis]) {
            best = curr;
            best_i = curr_i;
          }
        }
        h_loop_min_vertex.push_back(std::make_pair(&hole, best_i));
      }

      // sort the holes by the minimum vertex.
      std::sort(h_loop_min_vertex.begin(), h_loop_min_vertex.end(), order_h_loops<project_t, hole_t>(project, axis));

      // now, for each hole, find a vertex in the current polygon loop that it can be joined to.
      for (unsigned i = 0; i < h_loop_min_vertex.size(); ++i) {
        // the index of the vertex in the hole to connect.
        hole_vert_iter  h_loop_connect = h_loop_min_vertex[i].second;
        carve::geom2d::P2 hole_min = project(*h_loop_connect);

        f_loop_heap.clear();
        // we order polygon loop vertices that may be able to be connected
        // to the hole vertex by their distance to the hole vertex
        heap_ordering<project_t,vert_t> _heap_ordering(project, current_f_loop, *h_loop_connect);

        for (size_t j = 0; j < current_f_loop.size(); ++j) {
          // it is guaranteed that there exists a polygon vertex with
          // coord < the min hole coord chosen, which can be joined to
          // the min hole coord without crossing the polygon
          // boundary. also, because we merge holes in ascending
          // order, it is also true that this join can never cross
          // another hole (and that doesn't need to be tested for).
          if (project(current_f_loop[j]).v[axis] < hole_min.v[axis]) {
            f_loop_heap.push_back(j);
            std::push_heap(f_loop_heap.begin(), f_loop_heap.end(), _heap_ordering);
          }
        }

        // we are going to test each potential (according to the
        // previous test) polygon vertex as a candidate join. we order
        // by closeness to the hole vertex, so that the join we make
        // is as small as possible. to test, we need to check the
        // joining line segment does not cross any other line segment
        // in the current polygon loop (excluding those that have the
        // vertex that we are attempting to join with as an endpoint).
        while (f_loop_heap.size()) {
          std::pop_heap(f_loop_heap.begin(), f_loop_heap.end(), _heap_ordering);
          size_t curr = f_loop_heap.back();
          f_loop_heap.pop_back();
          // test the candidate join from current_f_loop[curr] to hole_min

          carve::geom2d::LineSegment2 test(hole_min, project(current_f_loop[curr]));

          size_t v1, v2;
          for (v1 = current_f_loop.size() - 1, v2 = 0; v2 != current_f_loop.size(); v1 = v2++) {
            // XXX: need to test vertices, not indices, because they may
            // be duplicated.
            if (current_f_loop[v1] == current_f_loop[curr] ||
                current_f_loop[v2] == current_f_loop[curr]) continue;
            carve::geom2d::LineSegment2 test2(project(current_f_loop[v1]), project(current_f_loop[v2]));
            carve::LineIntersectionClass ic = carve::geom2d::lineSegmentIntersection(test, test2).iclass;
            if (ic > 0) {
              // intersection; failed.
              goto intersection;
            }
          }

          patchHoleIntoPolygon(current_f_loop, curr, h_loop_min_vertex[i]);
          goto merged;

        intersection:;
        }
        ASSERT(!!!"didn't manage to link up hole!");

      merged:;
      }

      return current_f_loop;
    }



    struct tri_idx {
      unsigned a, b, c;
      tri_idx(unsigned _a, unsigned _b, unsigned _c) : a(_a), b(_b), c(_c) {
      }
    };

    void triangulate(const std::vector<carve::geom2d::P2> &poly,
                     std::vector<tri_idx> &result);
  }
}
