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

#include <carve/geom3d.hpp>

#include <carve/vertex_decl.hpp>
#include <carve/edge_decl.hpp>
#include <carve/face_decl.hpp>
#include <carve/octree_decl.hpp>
#include <carve/collection_types.hpp>

#include <assert.h>
#include <list>

namespace carve {
  namespace poly {

    struct Polyhedron {
    private:
      Polyhedron();
      Polyhedron(const Polyhedron &);
      Polyhedron &operator=(const Polyhedron &);

      bool initSpatialIndex();

    public:
      std::vector<bool> manifold_is_closed;
      std::vector<bool> manifold_is_negative;

      struct Connectivity {
        std::vector<std::vector<const Edge *> > vertex_to_edge;
        std::vector<std::vector<const Face *> > vertex_to_face;
        std::vector<std::vector<const Face *> > edge_to_face;
      } connectivity;

      std::vector<Vertex> poly_vertices;
      std::vector<Edge> edges;
      std::vector<Face> faces;

      carve::geom3d::AABB aabb;
      carve::csg::Octree octree;

      void invertAll();
      void invert(const std::vector<bool> &selected_manifolds);

      void invert(int m_id);
      void invert();

      Polyhedron *makeCopy(const std::vector<bool> &selected_manifolds) const;

      Polyhedron *makeCopy(int m_id) const;
      Polyhedron *makeCopy() const;

      const Face *connectedFace(const Face *, const Edge *) const;

      bool buildEdges();
      bool sortVertices();
      bool markManifolds();
      bool calcManifoldEmbedding();

      bool init();
      void faceRecalc();

      static Polyhedron *make(std::vector<Face> &_faces, bool _recalc = false);
      static Polyhedron *make(std::list<Face> &_faces, bool _recalc = false);

      static void collectFaceVertices(std::vector<Face> &faces,
                                      std::vector<Vertex> &vertices,
                                      carve::csg::VVMap &vmap);

      static void collectFaceVertices(std::vector<Face> &faces,
                                      std::vector<Vertex> &vertices);

      Polyhedron(std::vector<Face> &_faces, std::vector<Vertex> &_vertices, bool _recalc = false);
      Polyhedron(std::vector<Face> &_faces, bool _recalc = false);
      Polyhedron(std::list<Face> &_faces, bool _recalc = false);

      void commonFaceInit(bool _recalc);

      Polyhedron(const std::vector<carve::geom3d::Vector> &vertices, int n_faces, const std::vector<int> &face_indices);
      ~Polyhedron();

      void testVertexAgainstClosedManifolds(const carve::geom3d::Vector &v,
                                            std::map<int, PointClass> &result,
                                            bool ignore_orentation) const;

      PointClass containsVertex(const carve::geom3d::Vector &v,
                                const Face **hit_face = NULL,
                                bool even_odd = false,
                                int manifold_id = -1) const;

      void findEdgesNear(const carve::geom3d::LineSegment &l, std::vector<const Edge*> &edges) const;
      void findEdgesNear(const carve::geom3d::Vector &v, std::vector<const Edge*> &edges) const;
      void findEdgesNear(const Face &face, std::vector<const Edge*> &edges) const;
      void findEdgesNear(const Edge &edge, std::vector<const Edge*> &edges) const;

      void findFacesNear(const carve::geom3d::LineSegment &l, std::vector<const Face*> &faces) const;
      void findFacesNear(const Edge &edge, std::vector<const Face*> &faces) const;

      template<typename T>
      int vertexToEdges(const Vertex *v, T result) const;

      template<typename T>
      int edgeManifolds(const Edge *e, T result) const;

      template<typename T>
      int vertexManifolds(const Vertex *v, T result) const;

      template<typename T>
      int _faceNeighbourhood(const Face *f, int depth, T *result) const;

      template<typename T>
      int faceNeighbourhood(const Face *f, int depth, T result) const;

      template<typename T>
      int faceNeighbourhood(const Edge *e, int m_id, int depth, T result) const;

      template<typename T>
      int faceNeighbourhood(const Vertex *v, int m_id, int depth, T result) const;

      template<typename T>
      int edgeToFaces(const Edge *e, T result) const;

      template<typename T>
      int vertexToFaces(const Vertex *v, T result) const;

      void transform(const carve::math::Matrix &xform);

      template<typename T>
      void transform(const T &xform);

      void print(std::ostream &) const;

      ptrdiff_t vertexToIndex_fast(const Vertex *v) const;
      ptrdiff_t vertexToIndex(const Vertex *v) const;

      ptrdiff_t edgeToIndex_fast(const Edge *e) const;
      ptrdiff_t edgeToIndex(const Edge *e) const;

      size_t manifoldCount() const;

      bool hasOpenManifolds() const;

      void canonicalize();
    };

    std::ostream &operator<<(std::ostream &, const Polyhedron &);

  }
}
