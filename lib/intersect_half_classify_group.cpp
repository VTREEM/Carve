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
#include <carve/debug_hooks.hpp>

#include <list>
#include <set>
#include <iostream>

#include <algorithm>
#include <assert.h>

#include "intersect_common.hpp"
#include "intersect_classify_common.hpp"
#include "intersect_classify_common_impl.hpp"

namespace carve {
  namespace csg {

    namespace {
      struct GroupPoly : public CSG::Collector {
        const carve::poly::Polyhedron *want_groups_from;
        std::list<std::pair<FaceClass, carve::poly::Polyhedron *> > &out;

        GroupPoly(const carve::poly::Polyhedron *poly,
                  std::list<std::pair<FaceClass, carve::poly::Polyhedron *> > &_out) : CSG::Collector(), want_groups_from(poly), out(_out) {
        }

        virtual ~GroupPoly() {
        }

        virtual void collect(FaceLoopGroup *grp, CSG::Hooks &hooks) {
          if (grp->face_loops.head->orig_face->owner != want_groups_from) return;

          std::list<ClassificationInfo> &cinfo = (grp->classification);
          if (cinfo.size() == 0) {
            std::cerr << "WARNING! group " << grp << " has no classification info!" << std::endl;
            return;
          }
          // XXX: check all the cinfo elements for consistency.
          FaceClass fc = cinfo.front().classification;

          std::vector<carve::poly::Face<3> > faces;
          faces.reserve(grp->face_loops.size());
          for (FaceLoop *loop = grp->face_loops.head; loop != NULL; loop = loop->next) {
            faces.push_back(carve::poly::Face<3>());
            faces.back().init(loop->orig_face, loop->vertices, false);
          }

          out.push_back(std::make_pair(fc, new carve::poly::Polyhedron(faces)));
        }

        virtual carve::poly::Polyhedron *done(CSG::Hooks &hooks) {
          return NULL;
        }
      };

      class FaceMaker {
      public:
 
        bool pointOn(VertexClassification &vclass, FaceLoop *f, size_t index) const {
          return vclass[f->vertices[index]].cls[0] == POINT_ON;
        }

        void explain(FaceLoop *f, size_t index, PointClass pc) const {
#if defined(DEBUG)
          std::cerr << "face loop " << f << " from poly b is easy because vertex " << index << " (" << f->vertices[index]->v << ") is " << ENUM(pc) << std::endl;
#endif
        }
      };

      class HalfClassifyFaceGroups {
      public:
        std::list<std::pair<FaceClass, carve::poly::Polyhedron *> > &b_out;
        CSG::Hooks &hooks;

        HalfClassifyFaceGroups(std::list<std::pair<FaceClass, carve::poly::Polyhedron *> > &c, CSG::Hooks &h) : b_out(c), hooks(h) {
        }

        void classifySimple(FLGroupList &a_loops_grouped,
                            FLGroupList &b_loops_grouped,
                            VertexClassification &vclass,
                            const carve::poly::Polyhedron *poly_a,
                            const carve::poly::Polyhedron *poly_b) const {
          GroupPoly group_poly(poly_b, b_out);
          performClassifySimpleOnFaceGroups(a_loops_grouped, b_loops_grouped, poly_a, poly_b, group_poly, hooks);
#if defined(DEBUG)
          std::cerr << "after removal of simple on groups: " << b_loops_grouped.size() << " b groups" << std::endl;
#endif
        }

        void classifyEasy(FLGroupList &a_loops_grouped,
                          FLGroupList &b_loops_grouped,
                          VertexClassification &vclass,
                          const carve::poly::Polyhedron *poly_a,
                          const carve::poly::Polyhedron *poly_b) const {
          GroupPoly group_poly(poly_b, b_out);
          performClassifyEasyFaceGroups(b_loops_grouped, poly_a, vclass, FaceMaker(), group_poly, hooks);
#if defined(DEBUG)
          std::cerr << "after removal of easy groups: " << b_loops_grouped.size() << " b groups" << std::endl;
#endif
        }

        void classifyHard(FLGroupList &a_loops_grouped,
                          FLGroupList &b_loops_grouped,
                          VertexClassification &vclass,
                          const carve::poly::Polyhedron *poly_a,
                          const carve::poly::Polyhedron *poly_b) const {
          GroupPoly group_poly(poly_b, b_out);
          performClassifyHardFaceGroups(b_loops_grouped, poly_a, FaceMaker(), group_poly, hooks);
#if defined(DEBUG)
          std::cerr << "after removal of hard groups: " << b_loops_grouped.size() << " b groups" << std::endl;
#endif

        }

        void faceLoopWork(FLGroupList &a_loops_grouped,
                          FLGroupList &b_loops_grouped,
                          VertexClassification &vclass,
                          const carve::poly::Polyhedron *poly_a,
                          const carve::poly::Polyhedron *poly_b) const {
          GroupPoly group_poly(poly_b, b_out);
          performFaceLoopWork(poly_a, b_loops_grouped, *this, group_poly, hooks);
        }

        void postRemovalCheck(FLGroupList &a_loops_grouped,
                              FLGroupList &b_loops_grouped) const {
#if defined(DEBUG)
          std::cerr << "after removal of on groups: " << b_loops_grouped.size() << " b groups" << std::endl;
#endif
        }

        bool faceLoopSanityChecker(FaceLoopGroup &i) const {
          return false;
          return i.face_loops.size() != 1;
        }

        void finish(FLGroupList &a_loops_grouped,FLGroupList &b_loops_grouped) const {
#if defined(DEBUG)
          if (a_loops_grouped.size() || b_loops_grouped.size())
            std::cerr << "UNCLASSIFIED! a=" << a_loops_grouped.size() << ", b=" << b_loops_grouped.size() << std::endl;
#endif   
        }
      };
    }

    void CSG::halfClassifyFaceGroups(const V2Set &shared_edges,
                                     VertexClassification &vclass,
                                     const carve::poly::Polyhedron *poly_a,                           
                                     FLGroupList &a_loops_grouped,
                                     const LoopEdges &a_edge_map,
                                     const carve::poly::Polyhedron *poly_b,
                                     FLGroupList &b_loops_grouped,
                                     const LoopEdges &b_edge_map,
                                     std::list<std::pair<FaceClass, carve::poly::Polyhedron *> > &b_out) {
      HalfClassifyFaceGroups classifier(b_out, hooks);
      GroupPoly group_poly(poly_b, b_out);
      performClassifyFaceGroups(a_loops_grouped, b_loops_grouped, vclass, poly_a, poly_b, classifier, group_poly, hooks);
    }

  }
}
