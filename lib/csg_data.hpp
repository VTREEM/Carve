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

#include <carve/csg.hpp>

#include "internal_collection_types.hpp"

struct carve::csg::detail::Data {
//        * @param[out] vmap A mapping from vertex pointer to intersection point.
//        * @param[out] emap A mapping from edge pointer to intersection points.
//        * @param[out] fmap A mapping from face pointer to intersection points.
//        * @param[out] fmap_rev A mapping from intersection points to face pointers.
  // map from intersected vertex to intersection point.
  carve::csg::VVMap vmap;

  // map from intersected edge to intersection points.
  carve::detail::EVSMap emap;

  // map from intersected face to intersection points.
  carve::detail::FVSMap fmap;

  // map from intersection point to intersected faces.
  carve::detail::VFSMap fmap_rev;

  // created by divideEdges().
  // holds, for each edge, a 
  carve::detail::EVVMap divided_edges;

  // created by faceSplitEdges.
  carve::detail::FV2SMap face_split_edges;
};
