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

carve::poly::Polyhedron *readModel(const std::string &file) {
  carve::poly::Polyhedron *poly;

  if (endswith(file, ".ply")) {
    poly = readPLY(file);
  } else if (endswith(file, ".vtk")) {
    poly = readVTK(file);
  } else if (endswith(file, ".obj")) {
    poly = readOBJ(file);
  }
  if (poly == NULL) return NULL;

  std::cerr << "loaded polyhedron " << poly << " has "
    << poly->vertices.size() << " vertices "
    << poly->faces.size() << " faces "
    << poly->manifold_is_closed.size() << " manifolds (" << std::count(poly->manifold_is_closed.begin(), poly->manifold_is_closed.end(), true) << " closed)" << std::endl;

  return poly;
}

int main(int argc, char **argv) {
  options.parse(argc, argv);

  if (options.file == "") {
    exit(1);
  }

  carve::poly::Polyhedron *poly = readModel(options.file);

  std::vector<carve::poly::Vertex<3> > out_vertices = poly->vertices;
  std::vector<carve::poly::Face<3> > out_faces;

  size_t N = 0;
  for (size_t i = 0; i < poly->faces.size(); ++i) {
    carve::poly::Face<3> &f = poly->faces[i];
    N += f.vertices.size() - 2;
  }
  out_faces.reserve(N);

  for (size_t i = 0; i < poly->faces.size(); ++i) {
    carve::poly::Face<3> &f = poly->faces[i];
    if (f.vertices.size() == 3) {
      out_faces.push_back(carve::poly::Face<3>(
            &out_vertices[poly->vertexToIndex_fast(f.vertices[0])],
            &out_vertices[poly->vertexToIndex_fast(f.vertices[1])],
            &out_vertices[poly->vertexToIndex_fast(f.vertices[2])]
            ));
    } else {
      std::vector<carve::geom2d::P2> projected;
      projected.reserve(f.vertices.size());
      for (size_t j = 0; j < f.vertices.size(); ++j) {
        projected.push_back(f.project(f.vertices[j]->v));
      }

      std::vector<carve::triangulate::tri_idx> result;
      carve::triangulate::triangulate(projected, result);
      std::cerr << "triangulation, :" << projected.size() << " points -> " << result.size() << " triangles" << std::endl;
      for (size_t j = 0; j < projected.size(); ++j) {
        std::cerr << "    " << j << ": " << projected[j] << std::endl;
      }

      std::map<std::pair<size_t, size_t>, int> pairs;
      for (size_t j = 0; j < result.size(); ++j) {
        std::cerr << "   " << result[j].a << " " << result[j].b << " " << result[j].c << std::endl;
        if (result[j].a < result[j].b) pairs[std::make_pair(result[j].a, result[j].b)]++; else pairs[std::make_pair(result[j].b, result[j].a)]--;
        if (result[j].b < result[j].c) pairs[std::make_pair(result[j].b, result[j].c)]++; else pairs[std::make_pair(result[j].c, result[j].b)]--;
        if (result[j].c < result[j].a) pairs[std::make_pair(result[j].c, result[j].a)]++; else pairs[std::make_pair(result[j].a, result[j].c)]--;
        out_faces.push_back(carve::poly::Face<3>(
              &out_vertices[poly->vertexToIndex_fast(f.vertices[result[j].a])],
              &out_vertices[poly->vertexToIndex_fast(f.vertices[result[j].b])],
              &out_vertices[poly->vertexToIndex_fast(f.vertices[result[j].c])]
              ));
      }
      for (std::map<std::pair<size_t, size_t>, int>::const_iterator
             j = pairs.begin();
           j != pairs.end();
           ++j) {
        if ((*j).second) {
          std::cerr << "   | " << (*j).first.first << "-" << (*j).first.second << " " << (*j).second << std::endl;
        }
      }
    }
  }

  carve::poly::Polyhedron *result = new carve::poly::Polyhedron(out_faces, out_vertices);

  if (options.canonicalize) result->canonicalize();

  if (options.obj) {
    writeOBJ(std::cout, result);
  } else if (options.vtk) {
    writeVTK(std::cout, result);
  } else {
    writePLY(std::cout, result, options.ascii);
  }
}
