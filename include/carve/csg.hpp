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

#include <carve/geom3d.hpp>

#include <carve/poly.hpp>

#include <carve/collection_types.hpp>
#include <carve/classification.hpp>
#include <carve/iobj.hpp>
#include <carve/faceloop.hpp>
#include <carve/intersection.hpp>

namespace carve {
  namespace csg {

    class VertexPool {
      const static unsigned blocksize = 1024;
      typedef std::list<std::vector<carve::poly::Vertex<3> > > pool_t;
      pool_t pool;
    public:
      void reset();
      carve::poly::Vertex<3> *get(const carve::geom3d::Vector &v = carve::geom3d::Vector::ZERO());
      bool inPool(const carve::poly::Vertex<3> *v) const;

      VertexPool();
      ~VertexPool();
    };



    /** 
     * \class CSG
     * \brief The class responsible for the computation of CSG operations.
     * 
     */
    class CSG {
    public:
      struct Hook {
        /** 
         * \class Hook
         * \brief Provides API access to intermediate steps in CSG calculation.
         * 
         */
        virtual void intersectionVertex(const carve::poly::Vertex<3> *vertex,
                                        const IObjPairSet &intersections) {
        }
        virtual void processOutputFace(std::vector<carve::poly::Face<3> *> &faces,
                                       const carve::poly::Face<3> *orig_face,
                                       bool flipped) {
        }
        virtual void resultFace(const carve::poly::Face<3> *new_face,
                                const carve::poly::Face<3> *orig_face,
                                bool flipped) {
        }

        virtual ~Hook() {
        }
      };

        /** 
         * \class Hooks
         * \brief Management of API hooks.
         * 
         */
      class Hooks {
      public:
        enum {
          RESULT_FACE_HOOK         = 0,
          PROCESS_OUTPUT_FACE_HOOK = 1,
          INTERSECTION_VERTEX_HOOK = 2,
          HOOK_MAX                 = 3,

          RESULT_FACE_BIT          = 0x0001,
          PROCESS_OUTPUT_FACE_BIT  = 0x0002, 
          INTERSECTION_VERTEX_BIT  = 0x0004
       };

        std::vector<std::list<Hook *> > hooks;

        bool hasHook(unsigned hook_num);

        void intersectionVertex(const carve::poly::Vertex<3> *vertex,
                                const IObjPairSet &intersections);

        void processOutputFace(std::vector<carve::poly::Face<3> *> &faces,
                               const carve::poly::Face<3> *orig_face,
                               bool flipped);

        void resultFace(const carve::poly::Face<3> *new_face,
                        const carve::poly::Face<3> *orig_face,
                        bool flipped);

        void registerHook(Hook *hook, unsigned hook_bits);
        void unregisterHook(Hook *hook);

        void reset();

        Hooks();
        ~Hooks();
      };

        /** 
         * \class Collector
         * \brief Base class for objects responsible for selecting result from which form the result polyhedron.
         * 
         */
      class Collector {
        Collector(const Collector &);
        Collector &operator=(const Collector &);

      protected:

      public:
        virtual void collect(FaceLoopGroup *group, CSG::Hooks &) =0;
        virtual carve::poly::Polyhedron *done(CSG::Hooks &) =0;

        Collector() {}
        virtual ~Collector() {}
      };

    private:
      /// The computed intersection data.
      Intersections intersections;

      /// A map from intersection point to a set of intersections
      /// represented by pairs of intersection objects.
      VertexIntersections vertex_intersections;

      /// A pool from which temporary vertices are allocated. Also
      /// provides testing for pool membership.
      VertexPool vertex_pool;

      void init();

      void dumpIntersections();
      void makeVertexIntersections();

      const carve::poly::Vertex<3> *chooseWeldPoint(const VSet &equivalent);
      const carve::poly::Vertex<3> *weld(const VSet &equivalent);

      void groupIntersections();

      /** 
       * \brief Generate tables of intersecting pairs of faces.
       *
       * @param[out] vmap A mapping from vertex pointer to intersection point.
       * @param[out] emap A mapping from edge pointer to intersection points.
       * @param[out] fmap A mapping from face pointer to intersection points.
       * @param[out] fmap_rev A mapping from intersection points to face pointers.
       */
      void intersectingFacePairs(VVMap &vmap, EVSMap &emap, FVSMap &fmap, VFSMap &fmap_rev);

      void generateVertexEdgeIntersections(const carve::poly::Polyhedron *a, const carve::poly::Polyhedron *b);

      /** 
       * \brief Generate vertex-vertex, vertex-edge and edge-edge intersections between poly \a a and poly \a b.
       * 
       * @param a Polyhedron a.
       * @param b Polyhedron b.
       */
      void generateEdgeEdgeIntersections(const carve::poly::Polyhedron  *a, const carve::poly::Polyhedron  *b);

      /** 
       * \brief Generate vertex-face and edge-face intersections between poly \a a and poly \a b.
       * 
       * @param a Polyhedron a.
       * @param b Polyhedron b.
       */
      void generateEdgeFaceIntersections(const carve::poly::Polyhedron  *a, const carve::poly::Polyhedron  *b);

      void determinePotentiallyInteractingOctreeNodes(const carve::poly::Polyhedron *a, const carve::poly::Polyhedron *b);

      /** 
       * \brief Compute all points of intersection between poly \a a and poly \a b
       * 
       * @param a Polyhedron a.
       * @param b Polyhedron b.
       */
      void generateIntersections(const carve::poly::Polyhedron  *a, const carve::poly::Polyhedron  *b);

      /** 
       * \brief Divide edges in \a edges that are intersected by polyhedron \a poly
       * 
       * @param edges The edges to divide.
       * @param[in] poly The polyhedron to divide against.
       * @param[in] emap A mapping from edge pointer to intersections.
       * @param[out] face_edges A mapping from edge pointer to sets of ordered vertices.
       */
      void divideEdges(const std::vector<carve::poly::Edge<3> > &edges,
                       const carve::poly::Polyhedron  *poly,
                       const EVSMap &emap,
                       EVVMap &face_edges);

      /** 
       * \brief From the intersection points of pairs of intersecting faces, compute intersection edges.
       * 
       * @param[out] face_edges The edges which form the line of intersection on each face.
       * @param[out] eclass Classification information about created edges.
       * @param[in] fmap A mapping from face pointer to intersection vertices.
       * @param[in] fmap_rev The inverse of \a fmap.
       */
      void makeFaceEdges(FV2SMap &face_edges, EdgeClassification &eclass, FVSMap &fmap, VFSMap &fmap_rev);

      friend void classifyEasyFaces(FaceLoopList &face_loops,
                                    VertexClassification &vclass,
                                    const  carve::poly::Polyhedron  *other_poly,
                                    int other_poly_num,
                                    CSG &csg,
                                    CSG::Collector &collector);

      size_t generateFaceLoops(const carve::poly::Polyhedron  *poly,
                               const VVMap &vmap,
                               const FV2SMap &face_split_edges,
                               const EVVMap &divided_edges,
                               FaceLoopList &face_loops_out);


      // intersect_group.cpp

      /** 
       * \brief Build a loop edge mapping from a list of face loops.
       * 
       * @param[in] loops A list of face loops.
       * @param[in] edge_count A hint as to the number of edges in \a loops.
       * @param[out] edge_map The calculated map of edges to loops.
       */
      void makeEdgeMap(const FaceLoopList &loops, size_t edge_count, LoopEdges &edge_map);

      /** 
       * \brief Divide a list of face loops into groups that are connected by at least one edge not present in \a no_cross.
       * 
       * @param[in,out] face_loops The list of loops (will be emptied as a side effect)
       * @param[in] loop_edges A loop edge map used for traversing connected loops.
       * @param[in] no_cross A set of edges not to cross.
       * @param[out] out_loops A list of grouped face loops.
       */
      void groupFaceLoops(FaceLoopList &face_loops,
                          const LoopEdges &loop_edges,
                          const V2Set &no_cross,
                          FLGroupList &out_loops);

      /** 
       * \brief Find the set of edges shared between two edge maps.
       * 
       * @param[in] edge_map_a The first edge map.
       * @param[in] edge_map_b The second edge map.
       * @param[out] shared_edges The resulting set of common edges.
       */
      void findSharedEdges(LoopEdges &edge_map_a, LoopEdges &edge_map_b, V2Set &shared_edges);


      // intersect_classify_edge.cpp

      /** 
       * 
       * 
       * @param shared_edges 
       * @param vclass 
       * @param poly_a 
       * @param a_loops_grouped 
       * @param a_edge_map 
       * @param poly_b 
       * @param b_loops_grouped 
       * @param b_edge_map 
       * @param collector 
       */
      void classifyFaceGroupsEdge(const V2Set &shared_edges,
                                  VertexClassification &vclass,
                                  const carve::poly::Polyhedron  *poly_a,
                                  FLGroupList &a_loops_grouped,
                                  const LoopEdges &a_edge_map,
                                  const carve::poly::Polyhedron  *poly_b,
                                  FLGroupList &b_loops_grouped,
                                  const LoopEdges &b_edge_map,
                                  CSG::Collector &collector);

      // intersect_classify_group.cpp

      /** 
       * 
       * 
       * @param shared_edges 
       * @param vclass 
       * @param poly_a 
       * @param a_loops_grouped 
       * @param a_edge_map 
       * @param poly_b 
       * @param b_loops_grouped 
       * @param b_edge_map 
       * @param collector 
       */
      void classifyFaceGroups(const V2Set &shared_edges,
                              VertexClassification &vclass,
                              const carve::poly::Polyhedron  *poly_a, 
                              FLGroupList &a_loops_grouped,
                              const LoopEdges &a_edge_map,
                              const carve::poly::Polyhedron  *poly_b,
                              FLGroupList &b_loops_grouped,
                              const LoopEdges &b_edge_map,
                              CSG::Collector &collector);

      // intersect_half_classify_group.cpp

      /** 
       * 
       * 
       * @param shared_edges 
       * @param vclass 
       * @param poly_a 
       * @param a_loops_grouped 
       * @param a_edge_map 
       * @param poly_b 
       * @param b_loops_grouped 
       * @param b_edge_map 
       * @param FaceClass 
       * @param b_out 
       */
      void halfClassifyFaceGroups(const V2Set &shared_edges,
                                  VertexClassification &vclass,
                                  const carve::poly::Polyhedron  *poly_a, 
                                  FLGroupList &a_loops_grouped,
                                  const LoopEdges &a_edge_map,
                                  const carve::poly::Polyhedron  *poly_b,
                                  FLGroupList &b_loops_grouped,
                                  const LoopEdges &b_edge_map,
                                  std::list<std::pair<FaceClass, carve::poly::Polyhedron  *> > &b_out);

      // intersect.cpp

      /** 
       * \brief The main calculation method for CSG.
       * 
       * @param[in] a Polyhedron a
       * @param[in] b Polyhedron b
       * @param[out] vclass 
       * @param[out] eclass 
       * @param[out] a_face_loops 
       * @param[out] b_face_loops 
       * @param[out] a_edge_count 
       * @param[out] b_edge_count 
       */
      void calc(const carve::poly::Polyhedron  *a,
                const carve::poly::Polyhedron  *b,
                VertexClassification &vclass,
                EdgeClassification &eclass,
                FaceLoopList &a_face_loops,
                FaceLoopList &b_face_loops,
                size_t &a_edge_count,
                size_t &b_edge_count);

    public:
      /**
       * \enum OP
       * \brief Enumeration of the supported CSG operations.
       */
      enum OP {
        UNION,                  /**< in a or b. */
        INTERSECTION,           /**< in a and b. */
        A_MINUS_B,              /**< in a, but not b. */
        B_MINUS_A,              /**< in b, but not a. */
        SYMMETRIC_DIFFERENCE    /**< in a or b, but not both. */
      };

      /**
       * \enum CLASSIFY_TYPE
       * \brief The type of classification algorithm to use.
       */
      enum CLASSIFY_TYPE {
        CLASSIFY_NORMAL,        /**< Normal (group) classifier. */
        CLASSIFY_EDGE           /**< Edge classifier. */
      };

      CSG::Hooks hooks;         /**< The manager for calculation hooks. */

      CSG();
      ~CSG();

      /** 
       * \brief Compute a CSG operation between two polyhedra, \a a and \a b.
       * 
       * @param a Polyhedron a
       * @param b Polyhedron b
       * @param collector The collector (determines the CSG operation performed)
       * @param shared_edges A pointer to a set that will be populated with shared edges (if not NULL).
       * @param classify_type The type of classifier to use.
       * 
       * @return 
       */
      carve::poly::Polyhedron  *compute(const carve::poly::Polyhedron  *a,
                                        const carve::poly::Polyhedron  *b,
                                        CSG::Collector &collector,
                                        V2Set *shared_edges = NULL,
                                        CLASSIFY_TYPE classify_type = CLASSIFY_NORMAL);

      /** 
       * \brief Compute a CSG operation between two closed polyhedra, \a a and \a b.
       * 
       * @param a Polyhedron a
       * @param b Polyhedron b
       * @param op The CSG operation (A collector is created automatically).
       * @param shared_edges A pointer to a set that will be populated with shared edges (if not NULL).
       * @param classify_type The type of classifier to use.
       * 
       * @return 
       */
      carve::poly::Polyhedron  *compute(const carve::poly::Polyhedron  *a,
                                        const carve::poly::Polyhedron  *b,
                                        OP op,
                                        V2Set *shared_edges = NULL,
                                        CLASSIFY_TYPE classify_type = CLASSIFY_NORMAL);

      void slice(const carve::poly::Polyhedron  *a,
                 const carve::poly::Polyhedron  *b,
                 std::list<carve::poly::Polyhedron  *> &a_sliced,
                 std::list<carve::poly::Polyhedron  *> &b_sliced,
                 V2Set *shared_edges = NULL);

      bool sliceAndClassify(const carve::poly::Polyhedron  *closed,
                            const carve::poly::Polyhedron  *open,
                            std::list<std::pair<FaceClass, carve::poly::Polyhedron  *> > &result,
                            V2Set *shared_edges = NULL);
    };
  }
}
