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


#if defined(HAVE_CONFIG_H)
#  include <carve_config.h>
#endif

#include <carve/csg.hpp>
#include <carve/tree.hpp>
#include <carve/csg_triangulator.hpp>

#include "read_ply.hpp"
#include "write_ply.hpp"

#include "opts.hpp"

#include <fstream>
#include <algorithm>
#include <string>
#include <utility>
#include <set>
#include <iostream>
#include <iomanip>

#include <time.h>
typedef std::vector<std::string>::iterator TOK;



struct Options : public opt::Parser {
  bool improve;
  bool ascii;
  bool obj;
  bool vtk;
  bool canonicalize;

  std::string file;
  
  virtual void optval(const std::string &o, const std::string &v) {
    if (o == "--canonicalize" || o == "-c") { canonicalize = true; return; }
    if (o == "--binary"       || o == "-b") { ascii = false; return; }
    if (o == "--obj"          || o == "-O") { obj = true; return; }
    if (o == "--vtk"          || o == "-V") { vtk = true; return; }
    if (o == "--ascii"        || o == "-a") { ascii = true; return; }
    if (o == "--help"         || o == "-h") { help(std::cout); exit(0); }
    if (o == "--improve"      || o == "-i") { improve = true; return; }
  }

  virtual std::string usageStr() {
    return std::string ("Usage: ") + progname + std::string(" [options] expression");
  };

  virtual void arg(const std::string &a) {
    file = a;
  }

  virtual void help(std::ostream &out) {
    this->opt::Parser::help(out);
  }

  Options() {
    improve = false;
    ascii = true;
    obj = false;
    vtk = false;
    canonicalize = false;
    file = "";

    option("canonicalize", 'c', false, "Canonicalize before output (for comparing output).");
    option("binary",       'b', false, "Produce binary output.");
    option("ascii",        'a', false, "ASCII output (default).");
    option("obj",          'O', false, "Output in .obj format.");
    option("vtk",          'V', false, "Output in .vtk format.");
    option("improve",      'i', false, "Improve triangulation by minimising internal edge lengths.");
    option("help",         'h', false, "This help message.");
  }
};



static Options options;




static bool endswith(const std::string &a, const std::string &b) {
  if (a.size() < b.size()) return false;

  for (unsigned i = a.size(), j = b.size(); j; ) {
    if (tolower(a[--i]) != tolower(b[--j])) return false;
  }
  return true;
}

carve::mesh::MeshSet<3> *readModel(const std::string &file) {
  carve::mesh::MeshSet<3> *poly;

  if (file == "") {
    if (options.obj) {
      poly = readOBJasMesh(std::cin);
    } else if (options.vtk) {
      poly = readVTKasMesh(std::cin);
    } else {
      poly = readPLYasMesh(std::cin);
    }
  } else if (endswith(file, ".ply")) {
    poly = readPLYasMesh(file);
  } else if (endswith(file, ".vtk")) {
    poly = readVTKasMesh(file);
  } else if (endswith(file, ".obj")) {
    poly = readOBJasMesh(file);
  }

  if (poly == NULL) return NULL;

  std::cerr << "loaded polyhedron "
            << poly << " has " << poly->meshes.size()
            << " manifolds (" << std::count_if(poly->meshes.begin(),
                                               poly->meshes.end(),
                                               carve::mesh::Mesh<3>::IsClosed()) << " closed)" << std::endl; 

  std::cerr << "closed:    ";
  for (size_t i = 0; i < poly->meshes.size(); ++i) {
    std::cerr << (poly->meshes[i]->isClosed() ? '+' : '-');
  }
  std::cerr << std::endl;

  std::cerr << "negative:  ";
  for (size_t i = 0; i < poly->meshes.size(); ++i) {
    std::cerr << (poly->meshes[i]->isNegative() ? '+' : '-');
  }
  std::cerr << std::endl;

  return poly;
}

int main(int argc, char **argv) {
  options.parse(argc, argv);

  carve::mesh::MeshSet<3> *poly = readModel(options.file);
  if (!poly) exit(1);

  std::vector<carve::mesh::MeshSet<3>::face_t *> out_faces;

  size_t N = 0;
  for (carve::mesh::MeshSet<3>::face_iter i = poly->faceBegin(); i != poly->faceEnd(); ++i) {
    carve::mesh::MeshSet<3>::face_t *f = *i;
    N += f->nVertices() - 2;
  }
  out_faces.reserve(N);

  for (carve::mesh::MeshSet<3>::face_iter i = poly->faceBegin(); i != poly->faceEnd(); ++i) {
    carve::mesh::MeshSet<3>::face_t *f = *i;
    std::vector<carve::triangulate::tri_idx> result;

    std::vector<carve::mesh::MeshSet<3>::vertex_t *> vloop;
    f->getVertices(vloop);

    carve::triangulate::triangulate(
      carve::mesh::MeshSet<3>::face_t::projection_mapping(f->project),
      vloop,
      result);
    if (options.improve) {
      carve::triangulate::improve(
        carve::mesh::MeshSet<3>::face_t::projection_mapping(f->project),
        vloop,
        carve::mesh::vertex_distance(),
        result);
    }

    for (size_t j = 0; j < result.size(); ++j) {
      out_faces.push_back(
        new carve::mesh::MeshSet<3>::face_t(
          vloop[result[j].a],
          vloop[result[j].b],
          vloop[result[j].c]));
    }
  }

  carve::mesh::MeshSet<3> *result = new carve::mesh::MeshSet<3>(out_faces);

  if (options.canonicalize) result->canonicalize();

  if (options.obj) {
    writeOBJ(std::cout, result);
  } else if (options.vtk) {
    writeVTK(std::cout, result);
  } else {
    writePLY(std::cout, result, options.ascii);
  }

  delete result;
  delete poly;
}
