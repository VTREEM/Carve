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

namespace carve {
  namespace csg {
    template<typename filter_t>
    void Octree::doFindEdges(const carve::poly::Face &f,
                             Node *node,
                             std::vector<const carve::poly::Edge *> &out,
                             unsigned depth,
                             filter_t filter) const {
      if (node == NULL) {
        return;
      }

      if (node->aabb.intersects(f.aabb) && node->aabb.intersectsPlane(f.plane_eqn)) {
        if (node->hasChildren()) {
          for (int i = 0; i < 8; ++i) {
            doFindEdges(f, node->children[i], out, depth + 1, filter);
          }
        } else {
          if (depth < MAX_SPLIT_DEPTH && node->edges.size() > EDGE_SPLIT_THRESHOLD) {
            if (!node->split()) {
              for (int i = 0; i < 8; ++i) {
                doFindEdges(f, node->children[i], out, depth + 1, filter);
              }
              return;
            }
          }
          for (std::vector<const carve::poly::Edge*>::const_iterator it = node->edges.begin(), e = node->edges.end(); it != e; ++it) {
            if ((*it)->tag_once()) {
              if (filter(*it)) {
                out.push_back(*it);
              }
            }
          }
        }
      }
    }

    template<typename filter_t>
    void Octree::findEdgesNear(const carve::poly::Face &f, std::vector<const carve::poly::Edge *> &out, filter_t filter) const {
      tagable::tag_begin();
      doFindEdges(f, root, out, 0, filter);
    }

    inline void Octree::findEdgesNear(const carve::poly::Face &f, std::vector<const carve::poly::Edge *> &out) const {
      return findEdgesNear(f, out, no_filter());
    }

    template <typename FUNC>
    void Octree::doIterate(int level, Node *node, const FUNC &f) const{
      f(level, node);
      if (node->hasChildren()) {
        for (int i = 0; i < 8; ++i) {
          doIterate(level + 1, node->children[i], f);
        }
      }
    }

    template <typename FUNC>
    void Octree::iterateNodes(const FUNC &f) const {
      doIterate(0, root, f);
    }

  }
}
