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
#include <carve/timing.hpp>

#include "intersect_common.hpp"



void carve::csg::CSG::makeEdgeMap(const carve::csg::FaceLoopList &loops,
                                  size_t edge_count,
                                  carve::csg::LoopEdges &edge_map) {
#if defined(UNORDERED_COLLECTIONS_SUPPORT_RESIZE)
  edge_map.resize(edge_count);
#endif

  for (carve::csg::FaceLoop *i = loops.head; i; i = i->next) {
    edge_map.addFaceLoop(i);
    i->group = NULL;
  }
}

#include <carve/polyline.hpp>

void carve::csg::CSG::findSharedEdges(carve::csg::LoopEdges &edge_map_a,
                                      carve::csg::LoopEdges &edge_map_b,
                                      carve::csg::V2Set &shared_edges) {
  for (carve::csg::LoopEdges::const_iterator
         i = edge_map_a.begin(), e = edge_map_a.end();
       i != e;
       ++i) {
    carve::csg::LoopEdges::const_iterator j = edge_map_b.find((*i).first);
    if (j != edge_map_b.end()) {
      shared_edges.insert((*i).first);
    }
  }

#if defined(DEBUG)
  carve::csg::VVSMap edge_graph;

  for (carve::csg::V2Set::const_iterator i = shared_edges.begin(); i != shared_edges.end(); ++i) {
    edge_graph[(*i).first].insert((*i).second);
    edge_graph[(*i).second].insert((*i).first);
  }

  std::cerr << "*** testing consistency of edge graph" << std::endl;
  for (carve::csg::VVSMap::const_iterator i = edge_graph.begin(); i != edge_graph.end(); ++i) {
    if ((*i).second.size() > 2) {
      std::cerr << "branch at: " << (*i).first << std::endl;
    }
    if ((*i).second.size() == 1) {
      std::cerr << "endpoint at: " << (*i).first << std::endl;
      std::cerr << "coordinate: " << (*i).first->v << std::endl;
    }
  }

  {
    carve::line::PolylineSet intersection_graph;
    intersection_graph.vertices.resize(edge_graph.size());
    std::map<const carve::poly::Vertex<3> *, size_t> vmap;

    size_t j = 0;
    for (carve::csg::VVSMap::const_iterator i = edge_graph.begin(); i != edge_graph.end(); ++i) {
      intersection_graph.vertices[j].v = (*i).first->v;
      vmap[(*i).first] = j++;
    }

    while (edge_graph.size()) {
      carve::csg::VVSMap::iterator prior_i = edge_graph.begin();
      const carve::poly::Vertex<3> *prior = (*prior_i).first;
      std::vector<size_t> connected;
      connected.push_back(vmap[prior]);
      while (prior_i != edge_graph.end() && (*prior_i).second.size()) {
        const carve::poly::Vertex<3> *next = *(*prior_i).second.begin();
        carve::csg::VVSMap::iterator next_i = edge_graph.find(next);
        assert(next_i != edge_graph.end());
        connected.push_back(vmap[next]);
        (*prior_i).second.erase(next);
        (*next_i).second.erase(prior);
        if (!(*prior_i).second.size()) { edge_graph.erase(prior_i); prior_i = edge_graph.end(); }
        if (!(*next_i).second.size()) { edge_graph.erase(next_i); next_i = edge_graph.end(); }
        prior_i = next_i;
        prior = next;
      }
      bool closed = connected.front() == connected.back();
      for (size_t k = 0; k < connected.size(); ++k) {
        std::cerr << " " << connected[k];
      }
      std::cerr << std::endl;
      intersection_graph.addPolyline(closed, connected.begin(), connected.end());
    }

#if defined(DEBUG_WRITE_PLY_DATA)
    void writePLY(std::string &out_file, const carve::line::PolylineSet *lines, bool ascii);
    std::string out("/tmp/intersection.ply");
    writePLY(out, &intersection_graph, true);
#endif
  }

  std::cerr << "*** edge graph consistency test done" << std::endl;
#endif
}



#if defined(DEBUG)
static carve::poly::Polyhedron *groupToPolyhedron(const carve::csg::FaceLoopGroup &grp) {
  const carve::csg::FaceLoopList &fl = grp.face_loops;
  std::vector<carve::poly::Face<3> > faces;
  faces.reserve(fl.size());
  for (carve::csg::FaceLoop *f = fl.head; f; f = f->next) {
    faces.push_back(carve::poly::Face<3>());
    faces.back().init(f->orig_face, f->vertices, false);
  }
  carve::poly::Polyhedron *poly = new carve::poly::Polyhedron(faces);

  return poly;
}
#endif



void carve::csg::CSG::groupFaceLoops(carve::csg::FaceLoopList &face_loops,
                                     const carve::csg::LoopEdges &loop_edges,
                                     const carve::csg::V2Set &no_cross,
                                     carve::csg::FLGroupList &out_loops) {
  static carve::TimingName GROUP_FACE_LOOPS("groupFaceLoops()");

  carve::TimingBlock block(GROUP_FACE_LOOPS);

  int tag_num = 0;
  while (face_loops.size()) {
    out_loops.push_back(FaceLoopGroup());
    carve::csg::FaceLoopGroup &group = (out_loops.back());
    carve::csg::FaceLoopList &curr = (group.face_loops);
    carve::csg::V2Set &perim = (group.perimeter);

    carve::csg::FaceLoop *expand = face_loops.head;

    expand->group = &group;
    face_loops.remove(expand);
    curr.append(expand);

    while (expand) {
      std::vector<const carve::poly::Vertex<3> *> &loop = (expand->vertices);
      const carve::poly::Vertex<3> *v1, *v2;

      v1 = loop.back();
      for (size_t i = 0; i < loop.size(); ++i) {
        v2 = loop[i];

        carve::csg::V2Set::const_iterator nc = no_cross.find(std::make_pair(v1, v2));
        if (nc == no_cross.end()) {
          carve::csg::LoopEdges::const_iterator j;

          j = loop_edges.find(std::make_pair(v1, v2));
          if (j != loop_edges.end()) {
            for (std::list<carve::csg::FaceLoop *>::const_iterator
                   k = (*j).second.begin(), ke = (*j).second.end();
                 k != ke; ++k) {
              if ((*k)->group != NULL) continue;
              face_loops.remove((*k));
              curr.append((*k));
              (*k)->group = &group;
            }
          }

          j = loop_edges.find(std::make_pair(v2, v1));
          if (j != loop_edges.end()) {
            for (std::list<carve::csg::FaceLoop *>::const_iterator
                   k = (*j).second.begin(), ke = (*j).second.end();
                 k != ke; ++k) {
              if ((*k)->group != NULL) continue;
              face_loops.remove((*k));
              curr.append((*k));
              (*k)->group = &group;
            }
          }
        } else {
          perim.insert(std::make_pair(v1, v2));
        }
        v1 = v2;
      }
      expand = expand->next;
    }
    tag_num++;

#if defined(DEBUG)
    {
      carve::poly::Polyhedron *poly = groupToPolyhedron(group);
      char buf[128];
      sprintf(buf, "/tmp/group-%p.ply", &curr);
      std::string out(buf);
      void writePLY(std::string &out_file, const carve::poly::Polyhedron *poly, bool ascii);
      writePLY(out, poly, false);
      delete poly;
    }
#endif
  }
}
