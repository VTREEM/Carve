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
#include <carve/mesh.hpp>
#include <carve/mesh_impl.hpp>

#include "write_ply.hpp"

#include <vector>

void dumpMeshes(carve::mesh::MeshSet<3> *meshes) {
  std::cout << "*** meshes->meshes.size()=" << meshes->meshes.size() << std::endl;
  for (size_t i = 0; i < meshes->meshes.size(); ++i) {
    std::cout << "====mesh" << i << "====" << std::endl;
    carve::poly::Polyhedron *poly = carve::polyhedronFromMesh(meshes, i);
    writePLY(std::cout, poly, true);
  }
}

void make(const carve::geom::vector<3> *vec, size_t vec_count,
          const size_t *f_idx,
          std::vector<carve::mesh::Vertex<3> > &vertices,
          std::vector<carve::mesh::Face<3> *> &faces) {
  std::copy(vec, vec + vec_count, std::back_inserter(vertices));

  faces.clear();
  for (size_t i = 0; f_idx[i]; i += f_idx[i] + 1) {
    std::vector<carve::mesh::Vertex<3> *> v;
    v.reserve(f_idx[i]);
    for (size_t j = 0; j < f_idx[i]; ++j) {
      v.push_back(&vertices[f_idx[i+j+1]]);
    }
    faces.push_back(new carve::mesh::Face<3>(v.begin(), v.end()));
  }
}

void obj1(std::vector<carve::mesh::Vertex<3> > &vertices, std::vector<carve::mesh::Face<3> *> &faces) {
  const carve::geom::vector<3> vec[] = {
    carve::geom::VECTOR(-1,3,-1),
    carve::geom::VECTOR(3,-1,-1),
    carve::geom::VECTOR(1,3,1),
    carve::geom::VECTOR(1,1,-1),
    carve::geom::VECTOR(1,-1,-1),
    carve::geom::VECTOR(-1,-3,-1),
    carve::geom::VECTOR(3,-1,1),
    carve::geom::VECTOR(-1,1,1),
    carve::geom::VECTOR(1,3,-1),
    carve::geom::VECTOR(-1,1,-1),
    carve::geom::VECTOR(-3,1,-1),
    carve::geom::VECTOR(-3,1,1),
    carve::geom::VECTOR(1,-3,-1),
    carve::geom::VECTOR(-1,-3,1),
    carve::geom::VECTOR(-1,-1,-1),
    carve::geom::VECTOR(-3,-1,1),
    carve::geom::VECTOR(1,-1,1),
    carve::geom::VECTOR(1,-3,1),
    carve::geom::VECTOR(-3,-1,-1),
    carve::geom::VECTOR(1,1,1),
    carve::geom::VECTOR(-1,-1,1),
    carve::geom::VECTOR(3,1,-1),
    carve::geom::VECTOR(-1,3,1),
    carve::geom::VECTOR(3,1,1)
  };
  const size_t f_idx[] = {
    4,23,19,16,6,
    4,6,1,21,23,
    4,23,21,3,19,
    4,19,3,4,16,
    4,16,4,1,6,
    4,1,4,3,21,
    4,7,9,3,19,
    4,19,3,8,2,
    4,2,22,7,19,
    4,7,9,10,11,
    4,20,14,9,7,
    4,22,0,9,7,
    4,3,9,0,8,
    4,2,8,0,22,
    4,7,11,15,20,
    4,14,18,10,9,
    4,11,10,18,15,
    4,15,18,14,20,
    4,16,20,13,17,
    4,17,12,4,16,
    4,16,4,14,20,
    4,20,14,5,13,
    4,13,5,12,17,
    4,12,5,14,4,
    0
  };
  make(vec, sizeof(vec) / sizeof(vec[0]), f_idx, vertices, faces);
}

void obj2(std::vector<carve::mesh::Vertex<3> > &vertices, std::vector<carve::mesh::Face<3> *> &faces) {
  const carve::geom::vector<3> vec[] = {
    carve::geom::VECTOR(3,0.963138923901496135648869767465,1),
    carve::geom::VECTOR(1.41365128765474601379992236616,3,1),
    carve::geom::VECTOR(0.963138923901497356894196855137,-3,0),
    carve::geom::VECTOR(-3,1.41365128765474645788913221622,0),
    carve::geom::VECTOR(-1.68294196961579256566210460733,-2.16120922347255950413114078401,-1),
    carve::geom::VECTOR(-3,3,0),
    carve::geom::VECTOR(-1.68294196961579256566210460733,-2.16120922347255950413114078401,1),
    carve::geom::VECTOR(-2.16120922347255906004193093395,1.68294196961579300975131445739,1),
    carve::geom::VECTOR(2.16120922347255906004193093395,-1.68294196961579300975131445739,1),
    carve::geom::VECTOR(3,-1.41365128765474645788913221622,1),
    carve::geom::VECTOR(-0.963138923901496912804987005075,3,1),
    carve::geom::VECTOR(0,0,0),
    carve::geom::VECTOR(3,-1.41365128765474690197834206629,0),
    carve::geom::VECTOR(3.36588393923158601950262891478,1.08060461173627953002096546697,1),
    carve::geom::VECTOR(-3,-3,0),
    carve::geom::VECTOR(0.602337357879512258485021902743,5.52709316270414507954455984873,1),
    carve::geom::VECTOR(-1.08060461173628019615478024207,3.36588393923158557541341906472,1),
    carve::geom::VECTOR(-1.08060461173628019615478024207,3.36588393923158557541341906472,-1),
    carve::geom::VECTOR(-3,1.41365128765474645788913221622,1),
    carve::geom::VECTOR(-2.16120922347255906004193093395,1.68294196961579300975131445739,0),
    carve::geom::VECTOR(3,0.963138923901496468715777155012,0),
    carve::geom::VECTOR(5.52709316270414507954455984873,-0.602337357879513479730348990415,1),
    carve::geom::VECTOR(3,3,6),
    carve::geom::VECTOR(1.08060461173628019615478024207,-3.36588393923158557541341906472,-1),
    carve::geom::VECTOR(0.602337357879512258485021902743,5.52709316270414507954455984873,-1),
    carve::geom::VECTOR(-2.16120922347255906004193093395,1.68294196961579300975131445739,-1),
    carve::geom::VECTOR(-3,-0.963138923901496135648869767465,1),
    carve::geom::VECTOR(2.16120922347255906004193093395,-1.68294196961579300975131445739,0),
    carve::geom::VECTOR(-0.963138923901497356894196855137,3,0),
    carve::geom::VECTOR(-3.36588393923158601950262891478,-1.08060461173627953002096546697,-1),
    carve::geom::VECTOR(3,-3,0),
    carve::geom::VECTOR(-5.52709316270414507954455984873,0.602337357879513479730348990415,-1),
    carve::geom::VECTOR(-3,-3,6),
    carve::geom::VECTOR(1.68294196961579256566210460733,2.16120922347255950413114078401,-1),
    carve::geom::VECTOR(3.36588393923158601950262891478,1.08060461173627953002096546697,-1),
    carve::geom::VECTOR(0.963138923901496912804987005075,-3,1),
    carve::geom::VECTOR(1.68294196961579256566210460733,2.16120922347255950413114078401,0),
    carve::geom::VECTOR(1.41365128765474601379992236616,3,0),
    carve::geom::VECTOR(-0.602337357879512258485021902743,-5.52709316270414507954455984873,1),
    carve::geom::VECTOR(-3,3,6),
    carve::geom::VECTOR(-5.52709316270414507954455984873,0.602337357879513479730348990415,1),
    carve::geom::VECTOR(3,-3,6),
    carve::geom::VECTOR(0,0,-1),
    carve::geom::VECTOR(-0.602337357879512258485021902743,-5.52709316270414507954455984873,-1),
    carve::geom::VECTOR(0,0,1),
    carve::geom::VECTOR(-1.41365128765474579175531744113,-3,1),
    carve::geom::VECTOR(-1.68294196961579256566210460733,-2.16120922347255950413114078401,0),
    carve::geom::VECTOR(-1.41365128765474601379992236616,-3,0),
    carve::geom::VECTOR(-3,-0.963138923901496468715777155012,0),
    carve::geom::VECTOR(1.68294196961579256566210460733,2.16120922347255950413114078401,1),
    carve::geom::VECTOR(1.08060461173628019615478024207,-3.36588393923158557541341906472,1),
    carve::geom::VECTOR(3,3,0),
    carve::geom::VECTOR(5.52709316270414507954455984873,-0.602337357879513479730348990415,-1),
    carve::geom::VECTOR(2.16120922347255906004193093395,-1.68294196961579300975131445739,-1),
    carve::geom::VECTOR(-3.36588393923158601950262891478,-1.08060461173627953002096546697,1)
  };
  const size_t f_idx[] = {
    4,44,26,18,7,
    4,19,11,44,7,
    4,18,3,19,7,
    4,48,26,44,11,
    4,27,11,44,8,
    4,44,49,36,11,
    4,46,11,44,6,
    4,44,11,28,10,
    4,44,11,20,0,
    4,35,44,11,2,
    4,9,12,27,8,
    4,44,0,9,8,
    4,37,36,49,1,
    4,49,44,10,1,
    4,45,47,46,6,
    4,44,35,45,6,
    4,18,40,54,26,
    6,3,19,25,31,40,18,
    4,40,31,29,54,
    6,29,42,11,48,26,54,
    4,11,42,25,19,
    4,29,31,25,42,
    6,16,17,42,11,28,10,
    6,13,34,42,11,20,0,
    6,23,42,11,2,35,50,
    4,11,42,53,27,
    4,4,46,11,42,
    4,33,36,11,42,
    4,15,16,10,1,
    4,15,24,17,16,
    4,17,24,33,42,
    4,0,9,21,13,
    4,21,52,34,13,
    4,53,42,34,52,
    4,38,43,23,50,
    4,4,42,23,43,
    4,35,45,38,50,
    6,9,12,27,53,52,21,
    6,4,43,38,45,47,46,
    6,15,1,37,36,33,24,
    4,22,39,32,41,
    8,22,41,30,12,9,0,20,51,
    8,5,39,22,51,37,1,10,28,
    8,5,3,18,26,48,14,32,39,
    8,14,47,45,35,2,30,41,32,
    5,30,2,11,27,12,
    5,11,36,37,51,20,
    5,5,28,11,19,3,
    5,14,48,11,46,47,
    4,48,11,19,3,
    4,18,26,48,3,
    4,11,36,37,28,
    4,10,28,37,1,
    4,11,46,47,2,
    4,47,45,35,2,
    4,20,11,27,12,
    4,20,12,9,0,
    0
  };
  make(vec, sizeof(vec) / sizeof(vec[0]), f_idx, vertices, faces);
}

TEST(MeshTest, ComplexMesh) {
  std::vector<carve::mesh::Vertex<3> > vertices;
  std::vector<carve::mesh::Face<3> *> faces;
  obj1(vertices, faces);
  std::vector<carve::mesh::Mesh<3> *> meshes;
  carve::mesh::Mesh<3>::create(faces.begin(), faces.end(), meshes, carve::mesh::MeshOptions());
  ASSERT_EQ(meshes.size(), 1U);

  carve::mesh::MeshSet<3> *mesh = new carve::mesh::MeshSet<3>(vertices, meshes);
  carve::mesh::MeshSet<3> *mesh2 = mesh->clone();
  dumpMeshes(mesh);

  ASSERT_EQ(mesh->meshes.size(), 1U);
  ASSERT_EQ(mesh2->meshes.size(), 1U);

  delete mesh;
  delete mesh2;
}

TEST(MeshTest, ComplexMesh2) {
  std::vector<carve::mesh::Vertex<3> > vertices;
  std::vector<carve::mesh::Face<3> *> faces;
  obj2(vertices, faces);
  std::vector<carve::mesh::Mesh<3> *> meshes;
  carve::mesh::Mesh<3>::create(faces.begin(), faces.end(), meshes, carve::mesh::MeshOptions());
  ASSERT_EQ(meshes.size(), 5U);

  carve::mesh::MeshSet<3> *mesh = new carve::mesh::MeshSet<3>(vertices, meshes);
  carve::mesh::MeshSet<3> *mesh2 = mesh->clone();
  dumpMeshes(mesh);

  ASSERT_EQ(mesh->meshes.size(), 5U);
  ASSERT_EQ(mesh2->meshes.size(), 5U);

  delete mesh;
  delete mesh2;
}

TEST(MeshTest, MeshConstruction2) {
  std::vector<carve::mesh::Vertex<3> > vertices;
  vertices.reserve(9);
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(-1.0, -1.0, -1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(-1.0, +1.0, -1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(+1.0, +1.0, -1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(+1.0, -1.0, -1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(-1.0, -1.0, +1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(-1.0, +1.0, +1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(+1.0, +1.0, +1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(+1.0, -1.0, +1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR( 0.0,  0.0,  0.0)));

  std::vector<carve::mesh::Face<3> *> quadfaces;
  quadfaces.reserve(6);
  quadfaces.push_back(new carve::mesh::Face<3>(&vertices[0], &vertices[1], &vertices[2], &vertices[3]));
  quadfaces.push_back(new carve::mesh::Face<3>(&vertices[0], &vertices[4], &vertices[5], &vertices[1]));
  quadfaces.push_back(new carve::mesh::Face<3>(&vertices[1], &vertices[5], &vertices[6], &vertices[2]));
  quadfaces.push_back(new carve::mesh::Face<3>(&vertices[2], &vertices[6], &vertices[7], &vertices[3]));
  quadfaces.push_back(new carve::mesh::Face<3>(&vertices[3], &vertices[7], &vertices[4], &vertices[0]));
  quadfaces.push_back(new carve::mesh::Face<3>(&vertices[7], &vertices[6], &vertices[5], &vertices[4]));

  ASSERT_EQ(quadfaces[0]->plane.N, carve::geom::VECTOR( 0.0,  0.0, -1.0));
  ASSERT_EQ(quadfaces[1]->plane.N, carve::geom::VECTOR(-1.0,  0.0,  0.0));
  ASSERT_EQ(quadfaces[2]->plane.N, carve::geom::VECTOR( 0.0, +1.0,  0.0));
  ASSERT_EQ(quadfaces[3]->plane.N, carve::geom::VECTOR(+1.0,  0.0,  0.0));
  ASSERT_EQ(quadfaces[4]->plane.N, carve::geom::VECTOR( 0.0, -1.0,  0.0));
  ASSERT_EQ(quadfaces[5]->plane.N, carve::geom::VECTOR( 0.0,  0.0, +1.0));

  std::vector<carve::mesh::Mesh<3> *> quadmeshes;
  carve::mesh::Mesh<3>::create(quadfaces.begin(), quadfaces.end(), quadmeshes, carve::mesh::MeshOptions());
  ASSERT_EQ(quadmeshes.size(), 1U);
  carve::mesh::MeshSet<3> *mesh = new carve::mesh::MeshSet<3>(vertices, quadmeshes);
  dumpMeshes(mesh);
  delete mesh;
}

TEST(MeshTest, MeshConstruction1) {
  std::vector<carve::mesh::Vertex<3> > vertices;
  vertices.reserve(9);
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(-1.0, -1.0, -1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(-1.0, +1.0, -1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(+1.0, +1.0, -1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(+1.0, -1.0, -1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(-1.0, -1.0, +1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(-1.0, +1.0, +1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(+1.0, +1.0, +1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR(+1.0, -1.0, +1.0)));
  vertices.push_back(carve::mesh::Vertex<3>(carve::geom::VECTOR( 0.0,  0.0,  0.0)));

  std::vector<carve::mesh::Face<3> *> trifaces;
  trifaces.reserve(12);

  trifaces.push_back(new carve::mesh::Face<3>(&vertices[0], &vertices[1], &vertices[2]));
  trifaces.push_back(new carve::mesh::Face<3>(&vertices[3], &vertices[0], &vertices[2]));
  trifaces.push_back(new carve::mesh::Face<3>(&vertices[0], &vertices[4], &vertices[5]));
  trifaces.push_back(new carve::mesh::Face<3>(&vertices[1], &vertices[0], &vertices[5]));
  trifaces.push_back(new carve::mesh::Face<3>(&vertices[1], &vertices[5], &vertices[6]));
  trifaces.push_back(new carve::mesh::Face<3>(&vertices[2], &vertices[1], &vertices[6]));
  trifaces.push_back(new carve::mesh::Face<3>(&vertices[2], &vertices[6], &vertices[7]));
  trifaces.push_back(new carve::mesh::Face<3>(&vertices[3], &vertices[2], &vertices[7]));
  trifaces.push_back(new carve::mesh::Face<3>(&vertices[3], &vertices[7], &vertices[4]));
  trifaces.push_back(new carve::mesh::Face<3>(&vertices[0], &vertices[3], &vertices[4]));
  trifaces.push_back(new carve::mesh::Face<3>(&vertices[7], &vertices[6], &vertices[5]));
  trifaces.push_back(new carve::mesh::Face<3>(&vertices[4], &vertices[7], &vertices[5]));

  std::vector<carve::mesh::Mesh<3> *> trimeshes;
  carve::mesh::Mesh<3>::create(trifaces.begin(), trifaces.end(), trimeshes, carve::mesh::MeshOptions());
  ASSERT_EQ(trimeshes.size(), 1U);
  carve::mesh::MeshSet<3> *mesh = new carve::mesh::MeshSet<3>(vertices, trimeshes);
  dumpMeshes(mesh);
  delete mesh;
}
