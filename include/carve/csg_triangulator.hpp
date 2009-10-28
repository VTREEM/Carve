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

#include <carve/csg.hpp>
#include <carve/tag.hpp>
#include <carve/triangulator.hpp>

namespace carve {
  namespace csg {

    namespace detail {
      template<bool with_improvement>
      class CarveTriangulator : public csg::CSG::Hook {

      public:
        CarveTriangulator() {
        }

        virtual ~CarveTriangulator() {
        }

        virtual void processOutputFace(std::vector<poly::Face<3> *> &faces,
                                       const poly::Face<3> *orig,
                                       bool flipped) {
          std::vector<poly::Face<3> *> out_faces;

          size_t n_tris = 0;
          for (size_t f = 0; f < faces.size(); ++f) {
            ASSERT(faces[f]->vertices.size() >= 3);
            n_tris += faces[f]->vertices.size() - 2;
          }

          out_faces.reserve(n_tris);

          for (size_t f = 0; f < faces.size(); ++f) {
            poly::Face<3> *face = faces[f];

            if (face->vertices.size() == 3) {
              out_faces.push_back(face);
              continue;
            }

            std::vector<triangulate::tri_idx> result;

            triangulate::triangulate(poly::p2_adapt_project<3>(face->project), face->vertices, result);
            if (with_improvement) {
              triangulate::improve(poly::p2_adapt_project<3>(face->project), face->vertices, result);
            }

            std::vector<const poly::Vertex<3> *> fv;
            fv.resize(3);
            for (size_t i = 0; i < result.size(); ++i) {
              fv[0] = face->vertices[result[i].a];
              fv[1] = face->vertices[result[i].b];
              fv[2] = face->vertices[result[i].c];
              out_faces.push_back(face->create(fv, false));
            }
            delete face;
          }
          std::swap(faces, out_faces);
        }
      };
    }

    typedef detail::CarveTriangulator<false> CarveTriangulator;
    typedef detail::CarveTriangulator<true> CarveTriangulatorWithImprovement;

    class CarveTriangulationImprover : public csg::CSG::Hook {
    public:
      CarveTriangulationImprover() {
      }

      virtual ~CarveTriangulationImprover() {
      }

      virtual void processOutputFace(std::vector<poly::Face<3> *> &faces,
                                     const poly::Face<3> *orig,
                                     bool flipped) {
        if (faces.size() == 1) return;

        // doing improvement as a separate hook is much messier than
        // just incorporating it into the triangulation hook.

        typedef std::map<const poly::Vertex<3> *, size_t> vert_map_t;
        std::vector<poly::Face<3> *> out_faces;
        vert_map_t vert_map;

        out_faces.reserve(faces.size());

        poly::p2_adapt_project<3> projector(faces[0]->project);

        std::vector<triangulate::tri_idx> result;

        for (size_t f = 0; f < faces.size(); ++f) {
          poly::Face<3> *face = faces[f];
          if (face->vertices.size() != 3) {
            out_faces.push_back(face);
          } else {
            triangulate::tri_idx tri;
            for (size_t i = 0; i < 3; ++i) {
              size_t v = 0;
              vert_map_t::iterator j = vert_map.find(face->vertices[i]);
              if (j == vert_map.end()) {
                v = vert_map.size();
                vert_map[face->vertices[i]] = v;
              } else {
                v = (*j).second;
              }
              tri.v[i] = v;
            }
            result.push_back(tri);
            delete face;
          }
        }

        std::vector<const poly::Vertex<3> *> verts;
        verts.resize(vert_map.size());
        for (vert_map_t::iterator i = vert_map.begin(); i != vert_map.end(); ++i) {
          verts[(*i).second] = (*i).first;
        }
 
        triangulate::improve(projector, verts, result);

        std::vector<const poly::Vertex<3> *> fv;
        fv.resize(3);
        for (size_t i = 0; i < result.size(); ++i) {
          fv[0] = verts[result[i].a];
          fv[1] = verts[result[i].b];
          fv[2] = verts[result[i].c];
          out_faces.push_back(orig->create(fv, false));
        }

        std::swap(faces, out_faces);
      }
    };

    class CarveTriangulationQuadMerger : public csg::CSG::Hook {
      typedef std::map<V2, F2> edge_map_t;

    public:
      CarveTriangulationQuadMerger() {
      }

      virtual ~CarveTriangulationQuadMerger() {
      }

      double scoreQuad(edge_map_t::iterator i, edge_map_t &edge_map) {
        if (!(*i).second.first || !(*i).second.second) return -1;
      }

      poly::Face<3> *mergeQuad(edge_map_t::iterator i, edge_map_t &edge_map) {
      }

      void recordEdge(const poly::Vertex<3> *v1,
                      const poly::Vertex<3> *v2,
                      const poly::Face<3> *f,
                      edge_map_t &edge_map) {
        if (v1 < v2) {
          edge_map[V2(v1, v2)].first = f;
        } else {
          edge_map[V2(v2, v1)].second = f;
        }
      }

      virtual void processOutputFace(std::vector<poly::Face<3> *> &faces,
                                     const poly::Face<3> *orig,
                                     bool flipped) {
        if (faces.size() == 1) return;

        std::vector<poly::Face<3> *> out_faces;
        edge_map_t edge_map;

        out_faces.reserve(faces.size());

        poly::p2_adapt_project<3> projector(faces[0]->project);

        for (size_t f = 0; f < faces.size(); ++f) {
          poly::Face<3> *face = faces[f];
          if (face->vertices.size() != 3) {
            out_faces.push_back(face);
          } else {
            recordEdge(face->vertices[0], face->vertices[1], face, edge_map);
            recordEdge(face->vertices[1], face->vertices[2], face, edge_map);
            recordEdge(face->vertices[2], face->vertices[0], face, edge_map);
          }
        }

        for (edge_map_t::iterator i = edge_map.begin(); i != edge_map.end();) {
          if ((*i).second.first && (*i).second.second) {
            ++i;
          } else {
            edge_map.erase(i++);
          }
        }

        while (edge_map.size()) {
          edge_map_t::iterator i = edge_map.begin();
          edge_map_t::iterator best = i;
          double best_score = scoreQuad(i, edge_map);
          for (++i; i != edge_map.end(); ++i) {
            double score = scoreQuad(i, edge_map);
            if (score > best_score) best = i;
          }
          if (best_score < 0) break;
          out_faces.push_back(mergeQuad(best, edge_map));
        }

        if (edge_map.size()) {
          tagable::tag_begin();
          for (edge_map_t::iterator i = edge_map.begin(); i != edge_map.end(); ++i) {
            poly::Face<3> *a = const_cast<poly::Face<3> *>((*i).second.first);
            poly::Face<3> *b = const_cast<poly::Face<3> *>((*i).second.first);
            if (a && a->tag_once()) out_faces.push_back(a);
            if (b && b->tag_once()) out_faces.push_back(b);
          }
        }

        std::swap(faces, out_faces);
      }
    };
  }
}
