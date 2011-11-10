// Begin License:
// Copyright (C) 2006-2011 Tobias Sargeant (tobias.sargeant@gmail.com).
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

#include <gtest/gtest.h>

#if defined(HAVE_CONFIG_H)
#  include <carve_config.h>
#endif

#include <carve/carve.hpp>
#include <carve/csg.hpp>
#include <carve/input.hpp>

#include <map>

static carve::mesh::MeshSet<3> *makeCube(const carve::math::Matrix &transform) {
  carve::input::PolyhedronData data;

  data.addVertex(transform * carve::geom::VECTOR(+1.0, +1.0, +1.0));
  data.addVertex(transform * carve::geom::VECTOR(-1.0, +1.0, +1.0));
  data.addVertex(transform * carve::geom::VECTOR(-1.0, -1.0, +1.0));
  data.addVertex(transform * carve::geom::VECTOR(+1.0, -1.0, +1.0));
  data.addVertex(transform * carve::geom::VECTOR(+1.0, +1.0, -1.0));
  data.addVertex(transform * carve::geom::VECTOR(-1.0, +1.0, -1.0));
  data.addVertex(transform * carve::geom::VECTOR(-1.0, -1.0, -1.0));
  data.addVertex(transform * carve::geom::VECTOR(+1.0, -1.0, -1.0));

  data.addFace(0, 1, 2, 3);
  data.addFace(7, 6, 5, 4);
  data.addFace(0, 4, 5, 1);
  data.addFace(1, 5, 6, 2);
  data.addFace(2, 6, 7, 3);
  data.addFace(3, 7, 4, 0);

  return new carve::mesh::MeshSet<3>(data.points, data.getFaceCount(), data.faceIndices);
}

struct ResultFaceHook : public carve::csg::CSG::Hook {
  std::map<const carve::mesh::MeshSet<3> *, int> &counter;

  ResultFaceHook(std::map<const carve::mesh::MeshSet<3> *, int> &_counter) : counter(_counter) {
  }

  virtual void resultFace(const carve::mesh::MeshSet<3>::face_t *output_face,
      const carve::mesh::MeshSet<3>::face_t *source_face,
      bool flipped) {
    const carve::mesh::MeshSet<3> *source_poly = static_cast<const carve::mesh::MeshSet<3> *>(source_face->mesh->meshset);
    counter[source_poly]++;
  }
};

TEST(HookTest, ResultFace) {
  std::map<const carve::mesh::MeshSet<3> *, int> counter;
  carve::csg::CSG csg;

  carve::mesh::MeshSet<3> *a = makeCube(carve::math::Matrix::SCALE(+5, +5, .5));
  carve::mesh::MeshSet<3> *b = makeCube(carve::math::Matrix::ROT(.5, +1, +1, +1));

  csg.hooks.registerHook(new ResultFaceHook(counter), carve::csg::CSG::Hooks::RESULT_FACE_BIT);
  csg.compute(a, b, carve::csg::CSG::UNION, NULL, carve::csg::CSG::CLASSIFY_EDGE);

  ASSERT_EQ(counter.size(), 2);
  ASSERT_TRUE(counter.find(a) != counter.end());
  ASSERT_TRUE(counter.find(b) != counter.end());

  ASSERT_EQ(counter[a], 6);
  ASSERT_EQ(counter[b], 10);
}
