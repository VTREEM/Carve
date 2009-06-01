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

namespace carve {
  namespace csg {

    typedef std::pair<const carve::poly::Vertex<3> *, const carve::poly::Vertex<3> *> V2;

    static inline V2 unoriented_edge(const carve::poly::Vertex<3> *a, const carve::poly::Vertex<3> *b) {
      return V2(std::min(a, b), std::max(a, b));
    }

    static inline V2 flip(const V2 &v) {
      return V2(v.second, v.first);
    }

    typedef std::unordered_map<const carve::poly::Edge<3> *, std::vector<const carve::poly::Vertex<3> *>, carve::poly::hash_edge_ptr> EVVMap;

    typedef std::unordered_set<const carve::poly::Vertex<3> *, carve::poly::hash_vertex_ptr> VSet;
    typedef std::set<const carve::poly::Vertex<3> *> VSetSmall;

    typedef std::unordered_set<V2, carve::poly::hash_vertex_ptr> V2Set;
    typedef std::set<V2> V2SetSmall;

    typedef std::unordered_set<const carve::poly::Face<3> *, carve::poly::hash_face_ptr> FSet;
    typedef std::set<const carve::poly::Face<3> *> FSetSmall;

    typedef std::unordered_map<const carve::poly::Vertex<3> *, const carve::poly::Vertex<3> *, carve::poly::hash_vertex_ptr> VVMap;
    typedef std::unordered_map<const carve::poly::Vertex<3> *, VSetSmall, carve::poly::hash_vertex_ptr> VVSMap;
    typedef std::unordered_map<const carve::poly::Vertex<3> *, std::vector<const carve::poly::Edge<3> *>, carve::poly::hash_vertex_ptr> VEMap;
    typedef std::unordered_map<const carve::poly::Vertex<3> *, std::vector<const carve::poly::Face<3> *>, carve::poly::hash_vertex_ptr> VFMap;
    typedef std::unordered_map<const carve::poly::Vertex<3> *, FSetSmall, carve::poly::hash_vertex_ptr> VFSMap;

    typedef std::unordered_map<const carve::poly::Edge<3> *, VSetSmall, carve::poly::hash_edge_ptr> EVSMap;
    typedef std::unordered_map<const carve::poly::Edge<3> *, std::vector<const carve::poly::Face<3> *>, carve::poly::hash_edge_ptr> EFMap;

    typedef std::unordered_map<const carve::poly::Face<3> *, V2SetSmall, carve::poly::hash_face_ptr> FV2SMap;

    typedef std::unordered_map<const carve::poly::Face<3> *, VSetSmall, carve::poly::hash_face_ptr> FVSMap;

    std::ostream &operator<<(std::ostream &o, const FSet &s);
  }
}
