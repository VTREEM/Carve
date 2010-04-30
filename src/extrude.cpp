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

#include <carve/geom2d.hpp>
#include <carve/poly.hpp>
#include <carve/input.hpp>
#include <carve/triangulator.hpp>
#include <carve/matrix.hpp>

#include <fstream>
#include <sstream>

#include "write_ply.hpp"

#include <iostream>

void consume(std::istream &in, char ch) {
  while (in.good()) {
    int c;
    c = in.peek();
    if (in.eof()) return;
    if (std::isspace(c)) { in.ignore(); continue; }
    if (c == ch) { in.ignore(); }
    return;
  }
}

void readVals(std::istream &in, std::vector<double> &vals) {
  while (in.good()) {
    int c;

    while (std::isspace(in.peek())) in.ignore();

    c = in.peek();
    if (in.eof()) break;

    if (c != '-' && c != '+' && c != '.' && !std::isdigit(c)) {
      break;
    }

    double v;
    in >> v;
    vals.push_back(v);

    while (std::isspace(in.peek())) in.ignore();
    c = in.peek();

    if (c == ',') {
      in.ignore();
      continue;
    } else if (c == '-' || c == '+' || c == '.' || std::isdigit(c)) {
      continue;
    } else {
      std::cerr << "remain c=" << (char)c << std::endl;
      break;
    }
  }
}

void add(std::vector<carve::geom2d::P2> &points, carve::geom2d::P2 p) {
  if (!points.size() || points.back() != p) {
    points.push_back(p);
  }
}

void parsePath(std::istream &in, std::vector<std::vector<carve::geom2d::P2> > &paths) {
  double cx = 0.0;
  double cy = 0.0;

  std::vector<double> vals;
  std::vector<carve::geom2d::P2> points;
  while (in.good()) {
    char c;
    double x, y;
    in >> c;
    if (in.eof()) break;
    std::cerr << "c=" << c << std::endl;
    std::cerr << "points.size() == " << points.size() << std::endl;
    if (std::isspace(c)) continue;

    vals.clear();
    switch (c) {
      case 'M': {
        readVals(in, vals);
        points.clear();
        for (unsigned i = 0; i < vals.size(); i += 2) {
          add(points, carve::geom::VECTOR(vals[i], vals[i+1]));
        }
        cx = points.back().x;
        cy = points.back().y;
        break;
      }
      case 'm': {
        readVals(in, vals);
        points.clear();
        cx = 0.0;
        cy = 0.0;
        for (unsigned i = 0; i < vals.size(); i += 2) {
          add(points, carve::geom::VECTOR(cx+vals[i], cy+vals[i+1]));
          cx = points.back().x;
          cy = points.back().y;
        }
        break;
      }
      case 'L': {
        readVals(in, vals);
        for (unsigned i = 0; i < vals.size(); i += 2) {
          add(points, carve::geom::VECTOR(vals[i], vals[i+1]));
        }
        cx = points.back().x;
        cy = points.back().y;
        break;
      }
      case 'l': {
        readVals(in, vals);
        for (unsigned i = 0; i < vals.size(); i += 2) {
          add(points, carve::geom::VECTOR(cx+vals[i], cy+vals[i+1]));
          cx = points.back().x;
          cy = points.back().y;
        }
        break;
      }
      case 'H': {
        readVals(in, vals);
        for (unsigned i = 0; i < vals.size(); ++i) {
          add(points, carve::geom::VECTOR(vals[i], cy));
        }
        cx = points.back().x;
        break;
      }
      case 'h': {
        readVals(in, vals);
        for (unsigned i = 0; i < vals.size(); ++i) {
          add(points, carve::geom::VECTOR(cx+vals[i], cy));
          cx = points.back().x;
        }
        break;
      }
      case 'V': {
        readVals(in, vals);
        for (unsigned i = 0; i < vals.size(); ++i) {
          add(points, carve::geom::VECTOR(cx, vals[i]));
        }
        cy = points.back().y;
        break;
      }
      case 'v': {
        readVals(in, vals);
        for (unsigned i = 0; i < vals.size(); ++i) {
          add(points, carve::geom::VECTOR(cx, cy+vals[i]));
          cy = points.back().y;
        }
        break;
      }
      case 'Z':
      case 'z': {

        std::cerr << "path coords: " << std::endl;
        for (size_t i = 0; i < points.size(); ++i) {
          std::cerr << " " << i << ": " << points[i].x << "," << points[i].y;
        }
        if (points.back() == points.front()) points.pop_back();
        std::cerr << std::endl;
        paths.push_back(points);                     
        cx = 0.0;
        cy = 0.0;
        break;
      }
    }
  }
}

void parsePath(std::string &s, std::vector<std::vector<carve::geom2d::P2> > &paths) {
  std::istringstream in(s);
  return parsePath(in, paths);
}

carve::poly::Polyhedron *extrude(const std::vector<carve::geom2d::P2> &path, carve::geom3d::Vector dir) {
  const unsigned N = path.size();
  std::cerr << "N=" << N << std::endl;

  carve::input::PolyhedronData data;
  std::vector<unsigned> fwd, rev;
  fwd.reserve(N);
  rev.reserve(N);

  std::map<carve::geom3d::Vector, size_t> vert_idx;
  std::set<std::pair<size_t, size_t> > edges;
      
  for (size_t i = 0; i < path.size(); ++i) {
    carve::geom3d::Vector v;
    std::map<carve::geom3d::Vector, size_t>::iterator j;

    v = carve::geom::VECTOR(path[i].x, -path[i].y, 0.0);
    j = vert_idx.find(v);
    if (j == vert_idx.end()) {
      data.addVertex(v);
      fwd.push_back(vert_idx[v] = data.getVertexCount()-1);
    } else {
      fwd.push_back((*j).second);
    }

    v = carve::geom::VECTOR(path[i].x, -path[i].y, 0.0) + dir;
    j = vert_idx.find(v);
    if (j == vert_idx.end()) {
      data.addVertex(v);
      rev.push_back(vert_idx[v] = data.getVertexCount()-1);
    } else {
      rev.push_back((*j).second);
    }
  }

  data.addFace(fwd.begin(), fwd.end());
  data.addFace(rev.rbegin(), rev.rend());

  for (size_t i = 0; i < path.size()-1; ++i) {
    edges.insert(std::make_pair(fwd[i+1], fwd[i]));
  }
  edges.insert(std::make_pair(fwd[0], fwd[N-1]));

  for (size_t i = 0; i < path.size()-1; ++i) {
    if (edges.find(std::make_pair(fwd[i], fwd[i+1])) == edges.end()) {
      data.addFace(fwd[i+1], fwd[i], rev[i], rev[i+1]);
    }
  }
  if (edges.find(std::make_pair(fwd[N-1], fwd[0])) == edges.end()) {
    data.addFace(fwd[0], fwd[N-1], rev[N-1], rev[0]);
  }

  return new carve::poly::Polyhedron(data.points, data.getFaceCount(), data.faceIndices);
}

int main(int argc, char **argv) {
  typedef std::vector<carve::geom2d::P2> loop_t;

  std::ifstream in(argv[1]);
  unsigned file_num = 0;
  while (in.good()) {
    std::string s;
    std::getline(in, s);
    if (in.eof()) break;
    std::vector<std::vector<carve::geom2d::P2> > paths;
    parsePath(s, paths);

    std::cerr << "paths.size()=" << paths.size() << std::endl;

    if (paths.size() > 1) {
      std::vector<std::pair<size_t, size_t> > result;
      std::vector<carve::geom2d::P2> merged;

      result = carve::triangulate::incorporateHolesIntoPolygon(paths);

      merged.reserve(result.size());
      for (size_t i = 0; i < result.size(); ++i) {
        merged.push_back(paths[result[i].first][result[i].second]);
      }
      paths.clear();
      paths.push_back(merged);
    }

    std::vector<carve::geom2d::P2> &path = paths[0];
    std::cerr << "path coords: " << std::endl;
    for (size_t i = 0; i < path.size(); ++i) {
      std::cerr << " " << path[i].x << "," << path[i].y;
    }
    std::cerr << std::endl;

    std::cerr << "path.size()=" << path.size() << std::endl;
    carve::poly::Polyhedron *p = extrude(path, carve::geom::VECTOR(0,0,100));
    p->transform(carve::math::Matrix::TRANS(-p->aabb.pos));
    std::ostringstream outf;
    outf << "file_" << file_num++ << ".ply";
    std::ofstream out(outf.str().c_str());
    writePLY(out, p, false);
    delete p;
  }
  return 0;
}
