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

#include <carve/csg.hpp>
#include <carve/poly.hpp>
#include <carve/geom.hpp>

#include <vector>
#include <iostream>

int main()
{
  //create a tetrahedron
  std::vector<carve::mesh::MeshSet<3>::vertex_t> tet_verts;
  std::vector<carve::mesh::MeshSet<3>::face_t *> tet_faces;
  std::vector<carve::mesh::MeshSet<3>::vertex_t *> corners;

  tet_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(0.0, 0.0, 0.0)));
  tet_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(1.0, 0.0, 0.0)));
  tet_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(0.0, 1.0, 0.0)));
  tet_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(0.0, 0.0, 1.0)));
 
  corners.push_back(&tet_verts[0]);
  corners.push_back(&tet_verts[2]);
  corners.push_back(&tet_verts[1]);
  tet_faces.push_back(new carve::mesh::MeshSet<3>::face_t(corners.begin(), corners.end()));

  corners.clear();
  corners.push_back(&tet_verts[0]);
  corners.push_back(&tet_verts[1]);
  corners.push_back(&tet_verts[3]);
  tet_faces.push_back(new carve::mesh::MeshSet<3>::face_t(corners.begin(), corners.end()));

  corners.clear();
  corners.push_back(&tet_verts[0]);
  corners.push_back(&tet_verts[3]);
  corners.push_back(&tet_verts[2]);
  tet_faces.push_back(new carve::mesh::MeshSet<3>::face_t(corners.begin(), corners.end()));

  corners.clear();
  corners.push_back(&tet_verts[1]);
  corners.push_back(&tet_verts[2]);
  corners.push_back(&tet_verts[3]);
  tet_faces.push_back(new carve::mesh::MeshSet<3>::face_t(corners.begin(), corners.end()));

  carve::mesh::MeshSet<3> tetrahedron(tet_faces);

  //create a triangle
  std::vector<carve::mesh::MeshSet<3>::vertex_t> tri_verts;
  std::vector<carve::mesh::MeshSet<3>::face_t *> tri_faces;

  //Vertices
  //crashes if last coordinate set to 1e-8, but ok for 1e-7
  tri_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(-0.3, 0.0, 1e-8)));
  tri_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(1.0, 0.0, 1.1e-8)));
  tri_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(-0.3, 1.0, 1.1e-8)));

  //Face
  corners.clear();
  corners.push_back(&tri_verts[0]);
  corners.push_back(&tri_verts[2]);
  corners.push_back(&tri_verts[1]);
  tri_faces.push_back(new carve::mesh::MeshSet<3>::face_t(corners.begin(), corners.end()));

//  corners.clear();
//  corners.push_back(&tri_verts[0]);
//  corners.push_back(&tri_verts[1]);
//  corners.push_back(&tri_verts[2]);
//  tri_faces.push_back(new carve::mesh::MeshSet<3>::face_t(corners));

  carve::mesh::MeshSet<3> triangle(tri_faces);

  //cut triangle with tetrahedron.
  carve::mesh::MeshSet<3> *is_poly = carve::csg::CSG().compute(&tetrahedron,
                                                               &triangle,
                                                               carve::csg::CSG::INTERSECTION);

  // std::cout << "Tetrahedron is ... \n" << tetrahedron;
  // std::cout << "Triangle is ... \n" << triangle;
  // std::cout << "Intersection is ... \n" << *is_poly;

  return 0;
}
