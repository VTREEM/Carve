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

namespace std {
  template<unsigned ndim>
  inline void swap(carve::poly::Face<ndim> &a, carve::poly::Face<ndim> &b) {
    std::swap(a.vertices, b.vertices);
    std::swap(a.edges, b.edges);
    std::swap(a.owner, b.owner);
    std::swap(a.aabb, b.aabb);
    std::swap(a.plane_eqn, b.plane_eqn);
    std::swap(a.manifold_id, b.manifold_id);
    std::swap(a.project, b.project);
    std::swap(a.unproject, b.unproject);
  }
}

namespace carve {
  namespace poly {
  }
}
