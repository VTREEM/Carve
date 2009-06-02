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

#include <carve/geom.hpp>
#include <carve/aabb.hpp>
#include <carve/vertex_decl.hpp>
#include <carve/edge_decl.hpp>
#include <carve/face_decl.hpp>

namespace carve {

  namespace space {

    const static double SLACK_FACTOR = 1.0009765625;
    const static unsigned MAX_SPLIT_DEPTH = 32;

    struct nodedata_FaceEdge {
      std::vector<const carve::poly::Face<3> *> faces;
      std::vector<const carve::poly::Edge<3> *> edges;

      static bool overlap_test(const carve::geom::aabb<3> &aabb, const carve::poly::Face<3> *face) {
      }

      static bool overlap_test(const carve::geom::aabb<3> &aabb, const carve::poly::Edge<3> *edge) {
      }

      void add(const carve::poly::Face<3> *face) {
        faces.push_back(face);
      }

      void add(const carve::poly::Edge<3> *edge) {
        edges.push_back(edge);
      }

      template<typename iter_t>
      void _fetch(iter_t &iter, const carve::poly::Edge<3> *) {
        std::copy(edges.begin(), edges.end(), iter);
      }

      template<typename iter_t>
      void _fetch(iter_t &iter, const carve::poly::Face<3> *) {
        std::copy(faces.begin(), faces.end(), iter);
      }

      template<typename iter_t>
      void fetch(iter_t &iter) {
        return _fetch(iter, std::iterator_traits<iter_t>::value_type);
      }
    };

    template<unsigned n_dim, typename nodedata_t>
    class SpatialSubdivTree {

      typedef carve::geom::aabb<n_dim> aabb_t;
      typedef carve::geom::vector<n_dim> vector_t;

    public:
      class Node {
      public:
        Node *parent;
        Node *children;

        vector_t min;
        vector_t max;

        aabb_t aabb;

        nodedata_t data;

      private:
        Node(const Node &node); // undefined.
        Node &operator=(const Node &node); // undefined.

        Node() {
        }

        inline aabb_t makeAABB() const {
          vector_t centre = 0.5 * (min + max);
          vector_t size = SLACK_FACTOR * 0.5 * (max - min);
          return aabb_t(centre, size);
        }

        void setup(Node *_parent, const vector_t &_min, const vector_t &_max) {
          parent = _parent;
          min = _min;
          max = _max;
          aabb = makeAABB();
        }

        void alloc_children() {
          vector_t mid = 0.5 * (min + max);
          children = new Node[1 << n_dim];
          for (size_t i = 0; i < (1 << n_dim); ++i) {
            vector_t new_min, new_max;
            for (size_t c = 0; c < n_dim; ++c) {
              if (i & (1 << c)) {
                new_min.v[c] = min.v[c];
                new_max.v[c] = mid.v[c];
              } else {
                new_min.v[c] = mid.v[c];
                new_max.v[c] = max.v[c];
              }
            }
            children[i].setup(this, new_min, new_max);
          }
        }

        void dealloc_children() {
          delete [] children;
        }

      public:

        inline bool isLeaf() const { return children == NULL; }

        Node(Node *_parent, const vector_t &_min, const vector_t &_max) : parent(_parent), min(_min), max(_max), children(NULL) {
          aabb = makeAABB();
        }

        ~Node() {
          dealloc_children();
        }

        bool split() {
          if (isLeaf()) {
            alloc_children();
            // redistribute geometry to children.
          }
          return isLeaf();
        }

        // XXX: how do we manage node contents in n-dim?
        // std::vector<const carve::poly::Face<3> *> faces;
        // std::vector<const carve::poly::Edge<3> *> edges;
        // std::vector<const carve::poly::Vertex<3> *> vertices;

        // bool mightContain(const carve::poly::Face<3> &face);
        // bool mightContain(const carve::poly::Edge<3> &edge);
        // bool mightContain(const carve::poly::Vertex<3> &p);

        // bool hasGeometry();

        // template <class T>
        // void putInside(const T &input, Node *child, T &output);

      };



      Node *root;

      SpatialSubdivTree(const vector_t &_min, const vector_t &_max) : root(new Node(NULL, _min, _max)) {
      }

      ~SpatialSubdivTree() {
        delete root;
      }

      struct no_filter {
        template<typename obj_t>
        bool operator()(const obj_t &obj) const {
          return true;
        }
      };

      struct tag_filter {
        template<typename obj_t>
        bool operator()(const obj_t &obj) const {
          return obj.tag_once();
        }
      };

      template<typename obj_t, typename iter_t, typename filter_t>
      void _findObjectsNear(Node *node, const obj_t &object, iter_t &output, filter_t filter) {
        if (!node->isLeaf()) {
          for (size_t i = 0; i < (1 << n_dim); ++i) {
            if (nodedata_t::overlap_test(node->children[i].aabb, object)) {
              _findObjectsNear(node->children + i, object, output, filter);
            }
          }
          return;
        }
        node->data.fetch(output);
      }

      // in order to be used as an input, aabb_t::intersect(const obj_t &) must exist.
      template<typename obj_t, typename iter_t, typename filter_t>
      void findObjectsNear(const obj_t &object, iter_t output, filter_t filter) {
        if (!nodedata_t::overlap_test(root->aabb, object)) return;
        _findObjectsNear(root, object, output, filter);
      }
    };

  }
}
