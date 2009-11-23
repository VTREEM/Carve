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

#include <carve/polyhedron_base.hpp>

namespace carve {
  namespace csg {

    typedef std::pair<const carve::poly::Geometry<3>::vertex_t *, const carve::poly::Geometry<3>::vertex_t *> V2;
    typedef std::pair<const carve::poly::Geometry<3>::face_t *, const carve::poly::Geometry<3>::face_t *> F2;

    static inline V2 ordered_edge(const carve::poly::Geometry<3>::vertex_t *a, const carve::poly::Geometry<3>::vertex_t *b) {
      return V2(std::min(a, b), std::max(a, b));
    }

    static inline V2 flip(const V2 &v) {
      return V2(v.second, v.first);
    }

    typedef std::unordered_map<const carve::poly::Geometry<3>::edge_t *, std::vector<const carve::poly::Geometry<3>::vertex_t *>, carve::poly::hash_edge_ptr> EVVMap;

    typedef std::unordered_set<const carve::poly::Geometry<3>::vertex_t *, carve::poly::hash_vertex_ptr> VSet;
    typedef std::set<const carve::poly::Geometry<3>::vertex_t *> VSetSmall;

    typedef std::unordered_set<V2, carve::poly::hash_vertex_ptr> V2Set;
    typedef std::set<V2> V2SetSmall;

    typedef std::unordered_set<const carve::poly::Geometry<3>::face_t *, carve::poly::hash_face_ptr> FSet;
    typedef std::set<const carve::poly::Geometry<3>::face_t *> FSetSmall;

    typedef std::unordered_map<const carve::poly::Geometry<3>::vertex_t *, const carve::poly::Geometry<3>::vertex_t *, carve::poly::hash_vertex_ptr> VVMap;
    typedef std::unordered_map<const carve::poly::Geometry<3>::vertex_t *, VSetSmall, carve::poly::hash_vertex_ptr> VVSMap;
    typedef std::unordered_map<const carve::poly::Geometry<3>::vertex_t *, std::vector<const carve::poly::Geometry<3>::edge_t *>, carve::poly::hash_vertex_ptr> VEMap;
    typedef std::unordered_map<const carve::poly::Geometry<3>::vertex_t *, std::vector<const carve::poly::Geometry<3>::face_t *>, carve::poly::hash_vertex_ptr> VFMap;
    typedef std::unordered_map<const carve::poly::Geometry<3>::vertex_t *, FSetSmall, carve::poly::hash_vertex_ptr> VFSMap;

    typedef std::unordered_map<const carve::poly::Geometry<3>::edge_t *, VSetSmall, carve::poly::hash_edge_ptr> EVSMap;
    typedef std::unordered_map<const carve::poly::Geometry<3>::edge_t *, std::vector<const carve::poly::Geometry<3>::face_t *>, carve::poly::hash_edge_ptr> EFMap;

    typedef std::unordered_map<const carve::poly::Geometry<3>::face_t *, V2SetSmall, carve::poly::hash_face_ptr> FV2SMap;

    typedef std::unordered_map<const carve::poly::Geometry<3>::face_t *, VSetSmall, carve::poly::hash_face_ptr> FVSMap;

    std::ostream &operator<<(std::ostream &o, const FSet &s);
  }
}
