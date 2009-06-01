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
#include <carve/pointset.hpp>

#include <list>
#include <set>
#include <iostream>

#include <algorithm>
#include <assert.h>

#include "intersect_debug.hpp"
#include "intersect_common.hpp"
#include "intersect_classify_common.hpp"

#include "csg_collector.hpp"

#include <carve/timing.hpp>

typedef carve::poly::Polyhedron poly_t;

carve::csg::VertexPool::VertexPool() {
}

carve::csg::VertexPool::~VertexPool() {
}

void carve::csg::VertexPool::reset() {
  pool.clear();
}

poly_t::vertex_t *carve::csg::VertexPool::get(const carve::geom3d::Vector &v) {
  if (!pool.size() || pool.back().size() == blocksize) {
    pool.push_back(std::vector<poly_t::vertex_t>());
    pool.back().reserve(blocksize);
  }
  pool.back().push_back(poly_t::vertex_t(v));
  return &pool.back().back();
}

bool carve::csg::VertexPool::inPool(const poly_t::vertex_t *v) const {
  for (pool_t::const_iterator i = pool.begin(); i != pool.end(); ++i) {
    if (v >= &(i->front()) && v <= &(i->back())) return true;
  }
  return false;
}



namespace {
  static inline void HSV2RGB(float H, float S, float V, float &r, float &g, float &b) {
    H = 6.0f * H;
    if (S < 5.0e-6) {
      r = g = b = V; return;
    } else {
      int i = (int)H;
      float f = H - i;
      float p1 = V * (1.0f - S);
      float p2 = V * (1.0f - S * f);
      float p3 = V * (1.0f - S * (1.0f - f));
      switch (i) {
      case 0:  r = V;  g = p3; b = p1; return;
      case 1:  r = p2; g = V;  b = p1; return;
      case 2:  r = p1; g = V;  b = p3; return;
      case 3:  r = p1; g = p2; b = V;  return;
      case 4:  r = p3; g = p1; b = V;  return;
      case 5:  r = V;  g = p1; b = p2; return;
      }
    }
    r = g = b = 0.0;
  }



  /** 
   * \brief Sort a range [\a beg, \a end) of vertices in order of increasing dot product of vertex - \a base on \dir.
   * 
   * @tparam[in] T a forward iterator type.
   * @param[in] dir The direction in which to sort vertices.
   * @param[in] base 
   * @param[in] beg The start of the vertex range to sort.
   * @param[in] end The end of the vertex range to sort.
   * @param[out] out The sorted vertex result.
   * @param[in] size_hint A hint regarding the size of the output
   *            vector (to avoid needing to be able to calculate \a
   *            end - \a beg).
   */
  template<typename T>
  void orderVertices(const carve::geom3d::Vector &dir, const carve::geom3d::Vector &base,
                     T beg, const T end, std::vector<const poly_t::vertex_t *> &out,
                     size_t size_hint = 1) {
    typedef std::vector<std::pair<double, const poly_t::vertex_t *> > DVVector;
    std::vector<std::pair<double, const poly_t::vertex_t *> > ordered_vertices;
    ordered_vertices.reserve(size_hint);

    for (; beg != end; ++beg) {
      const poly_t::vertex_t *v = (*beg);
      ordered_vertices.push_back(std::make_pair(carve::geom::dot(v->v - base, dir), v));
    }

    std::sort(ordered_vertices.begin(), ordered_vertices.end());

    out.clear();
    out.reserve(ordered_vertices.size());
    for (DVVector::const_iterator
           i = ordered_vertices.begin(), e = ordered_vertices.end();
         i != e;
         ++i) {
      out.push_back((*i).second);
    }
  }



  /** 
   * 
   * 
   * @param dir 
   * @param base 
   * @param beg 
   * @param end 
   */
  template<typename T>
  void selectOrderingProjection(carve::geom3d::Vector &dir, carve::geom3d::Vector &base,
                                T beg, const T end) {
    double dx, dy, dz;
    const poly_t::vertex_t *min_x, *min_y, *min_z, *max_x, *max_y, *max_z;
    if (beg == end) return;
    min_x = max_x = min_y = max_y = min_z = max_z = *beg++;
    for (; beg != end; ++beg) {
      if (min_x->v.x > (*beg)->v.x) min_x = *beg;
      if (min_y->v.y > (*beg)->v.y) min_y = *beg;
      if (min_z->v.z > (*beg)->v.z) min_z = *beg;
      if (max_x->v.x < (*beg)->v.x) max_x = *beg;
      if (max_y->v.y < (*beg)->v.y) max_y = *beg;
      if (max_z->v.z < (*beg)->v.z) max_z = *beg;
    }

    dx = max_x->v.x - min_x->v.x;
    dy = max_y->v.y - min_y->v.y;
    dz = max_z->v.z - min_z->v.z;

    if (dx > dy) {
      if (dx > dz) {
        dir = max_x->v - min_x->v; base = min_x->v;
      } else {
        dir = max_z->v - min_z->v; base = min_z->v;
      }
    } else {
      if (dy > dz) {
        dir = max_y->v - min_y->v; base = min_y->v;
      } else {
        dir = max_z->v - min_z->v; base = min_z->v;
      }
    }
  }

}



#if defined(DEBUG)
std::ostream &carve::csg::operator<<(std::ostream &o, const carve::csg::IObj &a) {
  switch (a.obtype) {
  case carve::csg::IObj::OBTYPE_NONE:   o << "NONE{}"; break;
  case carve::csg::IObj::OBTYPE_VERTEX: o << "VERT{" << a.vertex << "}"; break;
  case carve::csg::IObj::OBTYPE_EDGE:   o << "EDGE{" << a.edge << "}"; break;
  case carve::csg::IObj::OBTYPE_FACE:   o << "FACE{" << a.face << "}"; break;
  }
  return o;
}



struct dump_data {
  const poly_t::vertex_t *i_pt;
  carve::csg::IObj i_src;
  carve::csg::IObj i_tgt;
  dump_data(const poly_t::vertex_t *_i_pt,
            carve::csg::IObj _i_src,
            carve::csg::IObj _i_tgt) : i_pt(_i_pt), i_src(_i_src), i_tgt(_i_tgt) {
  }
};



struct dump_sort {
  bool operator()(const dump_data &a, const dump_data &b) const {
    if (a.i_pt->v.x < b.i_pt->v.x) return true;
    if (a.i_pt->v.x > b.i_pt->v.x) return false;
    if (a.i_pt->v.y < b.i_pt->v.y) return true;
    if (a.i_pt->v.y > b.i_pt->v.y) return false;
    if (a.i_pt->v.z < b.i_pt->v.z) return true;
    if (a.i_pt->v.z > b.i_pt->v.z) return false;
    return false;
  }
};
#endif



void carve::csg::CSG::dumpIntersections() {
#if defined(DEBUG)
  std::vector<dump_data> temp;

  for (Intersections::const_iterator
         i = intersections.begin(),
         ie = intersections.end();
       i != ie;
       ++i) {
    const IObj &i_src = ((*i).first);

    for (Intersections::mapped_type::const_iterator
           j = (*i).second.begin(),
           je = (*i).second.end();
         j != je;
         ++j) {
      const IObj &i_tgt = ((*j).first);
      const poly_t::vertex_t *i_pt = ((*j).second);
      temp.push_back(dump_data(i_pt, i_src, i_tgt));
    }
  }

  std::sort(temp.begin(), temp.end(), dump_sort());

  for (size_t i = 0; i < temp.size(); ++i) {
    const IObj &i_src = temp[i].i_src;
    const IObj &i_tgt = temp[i].i_tgt;
    std::cerr << "INTERSECTION: " << temp[i].i_pt << " (" << temp[i].i_pt->v << ") is ";
    std::cerr << i_src << ".." << i_tgt << std::endl;
  }

  std::vector<carve::geom3d::Vector> vertices;
  for (Intersections::const_iterator
         i = intersections.begin(),
         ie = intersections.end();
       i != ie;
       ++i) {
    const IObj &i_src = ((*i).first);

    for (Intersections::mapped_type::const_iterator
           j = (*i).second.begin(),
           je = (*i).second.end();
         j != je;
         ++j) {
      const poly_t::vertex_t *i_pt = ((*j).second);
      vertices.push_back(i_pt->v);
    }
  }
  carve::point::PointSet points(vertices);

#if defined(DEBUG_WRITE_PLY_DATA)
  void writePLY(std::string &out_file, const carve::point::PointSet *points, bool ascii);
  std::string out("/tmp/intersection-points.ply");
  writePLY(out, &points, true);
#endif
#endif
}



bool carve::csg::CSG::Hooks::hasHook(unsigned hook_num) {
  return hooks[hook_num].size() > 0;
}

void carve::csg::CSG::Hooks::intersectionVertex(const poly_t::vertex_t *vertex,
                                                const IObjPairSet &intersections) {
  for (std::list<Hook *>::iterator j = hooks[INTERSECTION_VERTEX_HOOK].begin();
       j != hooks[INTERSECTION_VERTEX_HOOK].end();
       ++j) {
    (*j)->intersectionVertex(vertex, intersections);
  }
}

void carve::csg::CSG::Hooks::processOutputFace(std::vector<poly_t::face_t *> &faces,
                                               const poly_t::face_t *orig_face,
                                               bool flipped) {
  for (std::list<Hook *>::iterator j = hooks[PROCESS_OUTPUT_FACE_HOOK].begin();
       j != hooks[PROCESS_OUTPUT_FACE_HOOK].end();
       ++j) {
    (*j)->processOutputFace(faces, orig_face, flipped);
  }
}

void carve::csg::CSG::Hooks::resultFace(const poly_t::face_t *new_face,
                                        const poly_t::face_t *orig_face,
                                        bool flipped) {
  for (std::list<Hook *>::iterator j = hooks[RESULT_FACE_HOOK].begin();
       j != hooks[RESULT_FACE_HOOK].end();
       ++j) {
    (*j)->resultFace(new_face, orig_face, flipped);
  }
}

void carve::csg::CSG::Hooks::registerHook(Hook *hook, unsigned hook_bits) {
  for (unsigned i = 0; i < HOOK_MAX; ++i) {
    if (hook_bits & (1U << i)) {
      hooks[i].push_back(hook);
    }
  }
}

void carve::csg::CSG::Hooks::unregisterHook(Hook *hook) {
  for (unsigned i = 0; i < HOOK_MAX; ++i) {
    hooks[i].erase(std::remove(hooks[i].begin(), hooks[i].end(), hook), hooks[i].end());
  }
}

void carve::csg::CSG::Hooks::reset() {
  for (unsigned i = 0; i < HOOK_MAX; ++i) {
    for (std::list<Hook *>::iterator j = hooks[i].begin(); j != hooks[i].end(); ++j) {
      delete (*j);
    }
    hooks[i].clear();
  }
}

carve::csg::CSG::Hooks::Hooks() : hooks() {
  hooks.resize(HOOK_MAX);
}
 
carve::csg::CSG::Hooks::~Hooks() {
  reset();
}



void carve::csg::CSG::makeVertexIntersections() {
  static carve::TimingName FUNC_NAME("CSG::makeVertexIntersections()");
  carve::TimingBlock block(FUNC_NAME);
  vertex_intersections.clear();
  for (Intersections::const_iterator
         i = intersections.begin(),
         ie = intersections.end();
       i != ie;
       ++i) {
    const IObj &i_src = ((*i).first);

    for (Intersections::mapped_type::const_iterator
           j = (*i).second.begin(),
           je = (*i).second.end();
         j != je;
         ++j) {
      const IObj &i_tgt = ((*j).first);
      const poly_t::vertex_t *i_pt = ((*j).second);

      vertex_intersections[i_pt].insert(std::make_pair(i_src, i_tgt));
    }
  }
}



const poly_t::vertex_t *carve::csg::CSG::chooseWeldPoint(const carve::csg::VSet &equivalent) {
  // XXX: choose a better weld point.
  if (!equivalent.size()) return NULL;

  for (carve::csg::VSet::const_iterator
         i = equivalent.begin(), e = equivalent.end();
       i != e;
       ++i) {
    if (!vertex_pool.inPool((*i))) return (*i);
  }
  return *equivalent.begin();
}



const poly_t::vertex_t *carve::csg::CSG::weld(const carve::csg::VSet &equivalent) {
  const poly_t::vertex_t *weld_point = chooseWeldPoint(equivalent);

#if defined(DEBUG)
  std::cerr << "weld: " << equivalent.size() << " vertices ( ";
  for (carve::csg::VSet::const_iterator
         i = equivalent.begin(), e = equivalent.end();
       i != e;
       ++i) {
    const poly_t::vertex_t *v = (*i);
    std::cerr << " " << v;
  }
  std::cerr << ") to " << weld_point << std::endl;
#endif

  if (!weld_point) return NULL;

  carve::csg::VertexIntersections::mapped_type &weld_tgt = (vertex_intersections[weld_point]);

  for (carve::csg::VSet::const_iterator
         i = equivalent.begin(), e = equivalent.end();
       i != e;
       ++i) {
    const poly_t::vertex_t *v = (*i);

    if (v != weld_point) {
      carve::csg::VertexIntersections::iterator j = vertex_intersections.find(v);

      if (j != vertex_intersections.end()) {
        weld_tgt.insert((*j).second.begin(), (*j).second.end());
        vertex_intersections.erase(j);
      }
    }
  }
  return weld_point;
}



void carve::csg::CSG::groupIntersections() {
  static carve::TimingName GROUP_INTERSECTONS("groupIntersections()");

  carve::TimingBlock block(GROUP_INTERSECTONS);
  
  std::vector<const poly_t::vertex_t *> vertices;
  VVSMap graph;
#if 1 || defined(DEBUG)
  std::cerr << "groupIntersections()" << ": vertex_intersections.size()==" << vertex_intersections.size() << std::endl;
#endif

  vertices.reserve(vertex_intersections.size());
  for (carve::csg::VertexIntersections::const_iterator
         i = vertex_intersections.begin(),
         e = vertex_intersections.end();
       i != e;
       ++i) 
    {
      vertices.push_back((*i).first);
    }
  carve::geom3d::AABB aabb;
  aabb.fit(vertices.begin(), vertices.end(), carve::poly::vec_adapt_vertex_ptr());
  Octree vertex_intersections_octree;
  vertex_intersections_octree.setBounds(aabb);

  vertex_intersections_octree.addVertices(vertices);
      
  std::vector<const poly_t::vertex_t *> out;
  for (size_t i = 0, l = vertices.size(); i != l; ++i) {
    // let's find all the vertices near this one. 
    out.clear();
    vertex_intersections_octree.findVerticesNearAllowDupes(vertices[i]->v, out);

    for (size_t j = 0; j < out.size(); ++j) {
      if (vertices[i] != out[j] && carve::geom::equal(vertices[i]->v, out[j]->v)) {
#if defined(DEBUG)
        std::cerr << "EQ: " << vertices[i] << "," << out[j] << " " << vertices[i]->v << "," << out[j]->v << std::endl;
#endif
        graph[vertices[i]].insert(out[j]);
        graph[out[j]].insert(vertices[i]);
      }
    }
  }

  VSet visited, open;
  while (graph.size()) {
    visited.clear();
    open.clear();
    VVSMap::iterator i = graph.begin();
    open.insert((*i).first);
    while (open.size()) {
      VSet::iterator t = open.begin();
      const poly_t::vertex_t *o = (*t);
      open.erase(t);
      i = graph.find(o);
      ASSERT(i != graph.end());
      visited.insert(o);
      for (VVSMap::mapped_type::const_iterator
             j = (*i).second.begin(),
             je = (*i).second.end();
           j != je;
           ++j) {
        if (visited.count((*j)) == 0) {
          open.insert((*j));
        }
      }
      graph.erase(i);
    }
    weld(visited);
  }
}



void carve::csg::CSG::intersectingFacePairs(VVMap &vmap, EVSMap &emap, FVSMap &fmap, VFSMap &fmap_rev) {
  static carve::TimingName FUNC_NAME("CSG::intersectingFacePairs()");
  carve::TimingBlock block(FUNC_NAME);

  // iterate over all intersection points.
  for (carve::csg::VertexIntersections::const_iterator
         i = vertex_intersections.begin(),
         ie = vertex_intersections.end();
       i != ie;
       ++i) {
    const poly_t::vertex_t *i_pt = ((*i).first);
    carve::csg::VFSMap::mapped_type &face_set = (fmap_rev[i_pt]);

    // for all pairs of intersecting objects at this point
    for (carve::csg::VertexIntersections::mapped_type::const_iterator
           j = (*i).second.begin(),
           je = (*i).second.end();
         j != je;
         ++j) {
      const carve::csg::IObj &i_src = ((*j).first);
      const carve::csg::IObj &i_tgt = ((*j).second);

      // work out the faces involved. this updates fmap_rev.
      intersections.facesForObject(i_src, face_set);
      intersections.facesForObject(i_tgt, face_set);

      // record the intersection with respect to any involved vertex.
      if (i_src.obtype == IObj::OBTYPE_VERTEX) vmap[i_src.vertex] = i_pt;
      if (i_tgt.obtype == IObj::OBTYPE_VERTEX) vmap[i_tgt.vertex] = i_pt;

      // record the intersection with respect to any involved edge.
      if (i_src.obtype == IObj::OBTYPE_EDGE) emap[i_src.edge].insert(i_pt);
      if (i_tgt.obtype == IObj::OBTYPE_EDGE) emap[i_tgt.edge].insert(i_pt);
    }

    // record the intersection with respect to each face.
    for (carve::csg::VFSMap::mapped_type::const_iterator k = face_set.begin(), ke = face_set.end(); k != ke; ++k) {
      const poly_t::face_t *f = (*k);
      fmap[f].insert(i_pt);
    }
  }
}



void carve::csg::CSG::generateVertexEdgeIntersections(const poly_t *a, const poly_t *b) {
  static carve::TimingName FUNC_NAME("CSG::generateVertexEdgeIntersections()");
  carve::TimingBlock block(FUNC_NAME);

  std::vector<const poly_t::edge_t *> edges_in_b;
  for (size_t va_i = 0, va_l = a->vertices.size(); va_i != va_l; ++va_i) {
    const poly_t::vertex_t *v = &(a->vertices[va_i]);
    if (a->connectivity.vertex_to_face[a->vertexToIndex_fast(v)].size() == 0) {
      continue;
    }
    b->findEdgesNear(v->v, edges_in_b);
    // std::cerr << "testing vertex: " << v << " " << v->v << std::endl;

    for (size_t eb_i = 0, eb_l = edges_in_b.size(); eb_i != eb_l; ++eb_i) {
      const poly_t::edge_t *edge_b = edges_in_b[eb_i];
      const poly_t::vertex_t *ev1 = edge_b->v1, *ev2 = edge_b->v2;
      // std::cerr << "  aganist edge: " << edge_b << " [" << ev1 << "," << ev2 << "] " << ev1->v << " " << ev2->v << std::endl;
      if (intersections.intersects(v, edge_b)) {
        // std::cerr << "    already intersected" << std::endl;
        continue;
      }


      if (std::min(ev1->v.x, ev2->v.x) - carve::EPSILON > v->v.x ||
          std::max(ev1->v.x, ev2->v.x) + carve::EPSILON < v->v.x ||
          std::min(ev1->v.y, ev2->v.y) - carve::EPSILON > v->v.y ||
          std::max(ev1->v.y, ev2->v.y) + carve::EPSILON < v->v.y ||
          std::min(ev1->v.z, ev2->v.z) - carve::EPSILON > v->v.z ||
          std::max(ev1->v.z, ev2->v.z) + carve::EPSILON < v->v.z) {
        continue;
      }

      if (distance2(ev1->v, v->v) < carve::EPSILON2) {
        // vertex-vertex intersection
        intersections.record(IObj(ev1), IObj(v), v);
        // std::cerr << "INTERSECT(VV) " << v << "-" << ev1 << std::endl;
      } else if (distance2(ev2->v, v->v) < carve::EPSILON2) {
        // vertex-vertex intersection
        intersections.record(IObj(ev2), IObj(v), v);
        // std::cerr << "INTERSECT(VV) " << v << "-" << ev2 << std::endl;
      } else {
        double a = cross(ev2->v - ev1->v, v->v - ev1->v).length2();
        double b = (ev2->v - ev1->v).length2();
        if (a < b * carve::EPSILON2) {
          // vertex-edge intersection
          intersections.record(IObj(edge_b), IObj(v), v);
          // std::cerr << "INTERSECT(VE) " << v << "-" << edge_b << std::endl;
        }
      }

    }
  }
}



void carve::csg::CSG::generateEdgeEdgeIntersections(const poly_t *a, const poly_t *b) {
  static carve::TimingName FUNC_NAME("CSG::generateEdgeEdgeIntersections()");
  carve::TimingBlock block(FUNC_NAME);

  std::vector<const poly_t::edge_t *> edges_in_b;
  for (size_t ea_i = 0, ea_l = a->edges.size(); ea_i != ea_l; ++ea_i) {
    const poly_t::edge_t *edge_a = &a->edges[ea_i];
    const poly_t::vertex_t *v1 = edge_a->v1, *v2 = edge_a->v2;

    b->findEdgesNear(*edge_a, edges_in_b);

    for (size_t eb_i = 0, eb_l = edges_in_b.size(); eb_i != eb_l; ++eb_i) {
      const poly_t::edge_t *edge_b = edges_in_b[eb_i];
      const poly_t::vertex_t *v3 = edge_b->v1, *v4 = edge_b->v2;

      if (intersections.intersects(edge_a, edge_b)) {
        continue;
      }

      if (std::max(v3->v.x, v4->v.x) + carve::EPSILON < std::min(v1->v.x, v2->v.x) - carve::EPSILON ||
          std::max(v1->v.x, v2->v.x) + carve::EPSILON < std::min(v3->v.x, v4->v.x) - carve::EPSILON) continue;
      if (std::max(v3->v.y, v4->v.y) + carve::EPSILON < std::min(v1->v.y, v2->v.y) - carve::EPSILON ||
          std::max(v1->v.y, v2->v.y) + carve::EPSILON < std::min(v3->v.y, v4->v.y) - carve::EPSILON) continue;
      if (std::max(v3->v.z, v4->v.z) + carve::EPSILON < std::min(v1->v.z, v2->v.z) - carve::EPSILON ||
          std::max(v1->v.z, v2->v.z) + carve::EPSILON < std::min(v3->v.z, v4->v.z) - carve::EPSILON) continue;

      carve::geom3d::Vector p1, p2;
      double mu1, mu2;

      switch (carve::geom3d::rayRayIntersection(carve::geom3d::Ray(v2->v - v1->v, v1->v),
                                                carve::geom3d::Ray(v4->v - v3->v, v3->v),
                                                p1, p2, mu1, mu2)) {
      case RR_INTERSECTION: {
        // edges intersect

        carve::geom3d::Vector p1, p2;
        double mu1, mu2;

        // std::cerr << "edge intersect: " << v1 << " " << v2 << "   " << v3 << " " << v4 << std::endl;
        if (!carve::geom3d::rayRayIntersection(carve::geom3d::Ray(v2->v - v1->v, v1->v),
                                               carve::geom3d::Ray(v4->v - v3->v, v3->v),
                                               p1, p2, mu1, mu2) ||
            !carve::geom::equal(p1, p2)) {
          continue;
        }

        if (mu1 >= 0.0 && mu1 <= 1.0 && mu2 >= 0.0 && mu2 <= 1.0) {
          IObj o1, o2;
          const poly_t::vertex_t *p;

          o1 = IObj(edge_a);
          o2 = IObj(edge_b);
          p = vertex_pool.get((p1 + p2) / 2.0);
          intersections.record(o1, o2, p);
        }
      }
      case RR_PARALLEL: {
        // edges parallel. any intersection of this type should have
        // been handled by generateVertexEdgeIntersections().
        break;
      }
      case RR_DEGENERATE: {
        throw carve::exception("degenerate edge");
        break;
      }
      case RR_NO_INTERSECTION: {
        break;
      }
      }
    }
  }
}



void carve::csg::CSG::generateEdgeFaceIntersections(const poly_t *a, const poly_t *b) {
  static carve::TimingName FUNC_NAME("CSG::generateEdgeFaceIntersections()");
  carve::TimingBlock block(FUNC_NAME);

  carve::csg::FSet if_e;

  std::vector<const poly_t::edge_t *> edges_in_b;

  for (size_t fa_i = 0, fa_l = a->faces.size(); fa_i != fa_l; ++fa_i) {
    const poly_t::face_t &face_a = a->faces[fa_i];

    b->findEdgesNear(face_a, edges_in_b);

    // vertex-face intersections.
    for (size_t eb_i = 0, eb_l = edges_in_b.size(); eb_i != eb_l; ++eb_i) {
      const poly_t::edge_t *edge_b = (edges_in_b[eb_i]);
      if (edge_b == NULL) continue;
      double d1 = carve::geom::distance(face_a.plane_eqn, edge_b->v1->v);
      double d2 = carve::geom::distance(face_a.plane_eqn, edge_b->v2->v);

      // shortcircuit: does the edge cross the face?
      if (std::max(d1, d2) < -carve::EPSILON || std::min(d1, d2) > carve::EPSILON) { edges_in_b[eb_i] = NULL; continue; }

      if (fabs(d1) < carve::EPSILON &&
          !intersections.intersects(edge_b->v1, &face_a) &&
          face_a.containsPoint(edge_b->v1->v)) {
        intersections.record(edge_b->v1, &face_a, edge_b->v1);
        for (size_t eb_j = eb_i + 1; eb_j != eb_l; ++eb_j) {
          if (edges_in_b[eb_j] &&
              (edges_in_b[eb_j]->v1 == edge_b->v1 ||
               edges_in_b[eb_j]->v2 == edge_b->v1)) edges_in_b[eb_j] = NULL;
        }
        edges_in_b[eb_i] = NULL;
      }

      if (fabs(d2) < carve::EPSILON &&
          !intersections.intersects(edge_b->v2, &face_a) &&
          face_a.containsPoint(edge_b->v2->v)) {
        intersections.record(edge_b->v2, &face_a, edge_b->v2);
        for (size_t eb_j = eb_i + 1; eb_j != eb_l; ++eb_j) {
          if (edges_in_b[eb_j] &&
              (edges_in_b[eb_j]->v1 == edge_b->v2 ||
               edges_in_b[eb_j]->v2 == edge_b->v2)) edges_in_b[eb_j] = NULL;
        }
        edges_in_b[eb_i] = NULL;
      }
    }

    // edge-face intersections.
    for (size_t eb_i = 0, eb_l = edges_in_b.size(); eb_i != eb_l; ++eb_i) {
      const poly_t::edge_t *edge_b = (edges_in_b[eb_i]);
      if (edge_b == NULL) continue;

      if (intersections.intersects(edge_b, &face_a)) continue;

      carve::geom3d::Vector p;
      if (face_a.simpleLineSegmentIntersection(carve::geom3d::LineSegment(edge_b->v1->v, edge_b->v2->v), p)) {
        intersections.record(edge_b, &face_a, vertex_pool.get(p));
      }
    }
  }
}

static void dump_octree_stats(carve::csg::Octree::Node *node, size_t depth, std::string indent = "") {
  if (node->is_leaf) {
    std::cerr << indent << node << "." << depth << " " << node->faces.size() << " faces " << node->edges.size() << " edges " << node->vertices.size() << " vertices" << std::endl;
  } else {
    std::cerr << indent << node << "." << depth << std::endl;
    for (size_t i = 0; i < 8; ++i) {
      dump_octree_stats(node->children[i], depth+1, indent + "  ");
    }
  }
}

struct hash_node {
  size_t operator()(const carve::csg::Octree::Node * const &node) const { return size_t(node); }
};

void carve::csg::CSG::determinePotentiallyInteractingOctreeNodes(const poly_t *a, const poly_t *b) {
}

void carve::csg::CSG::generateIntersections(const poly_t *a, const poly_t *b) {
  generateVertexEdgeIntersections(a, b);
  generateVertexEdgeIntersections(b, a);

#if defined(DEBUG)
  std::cerr << "generateEdgeEdgeIntersections" << std::endl;
#endif
  generateEdgeEdgeIntersections(a, b);

#if defined(DEBUG)
  std::cerr << std::endl << "after edge-edge" << std::endl;
  dumpIntersections();
  std::cerr << std::endl;
#endif

#if defined(DEBUG)
  std::cerr << "generateEdgeFaceIntersections" << std::endl;
#endif
  generateEdgeFaceIntersections(a, b);
  generateEdgeFaceIntersections(b, a);

#if defined(DEBUG)
  std::cerr << std::endl << "after edge-face" << std::endl;
  dumpIntersections();
  std::cerr << std::endl;
#endif
 
#if defined(DEBUG)
  std::cerr << "makeVertexIntersections" << std::endl;
#endif
  makeVertexIntersections();

#if defined(DEBUG)
  std::cerr << "  intersections.size() " << intersections.size() << std::endl;
  map_histogram(intersections);
  std::cerr << "  vertex_intersections.size() " << vertex_intersections.size() << std::endl;
  map_histogram(vertex_intersections);
#endif

#if defined(DEBUG) && defined(DEBUG_DRAW_INTERSECTIONS)
  HOOK(drawIntersections(vertex_intersections););
#endif

#if defined(DEBUG)
  // std::cerr << "groupIntersections" << std::endl;
#endif
  //groupIntersections();

#if defined(DEBUG)
  std::cerr << "  intersections.size() " << intersections.size() << std::endl;
  std::cerr << "  vertex_intersections.size() " << vertex_intersections.size() << std::endl;
#endif

  // notify about intersections.
  if (hooks.hasHook(Hooks::INTERSECTION_VERTEX_HOOK)) {
    for (VertexIntersections::const_iterator i = vertex_intersections.begin();
         i != vertex_intersections.end();
         ++i) {
      hooks.intersectionVertex((*i).first, (*i).second);
    }
  }

  // from here on, only vertex_intersections is used for intersection
  // information.

  // intersections still contains the vertex_to_face map. maybe that
  // should be moved out into another class.
  static_cast<Intersections::super>(intersections).clear();
}



carve::csg::CSG::CSG() {
}



carve::csg::CSG::~CSG() {
}



void carve::csg::CSG::divideEdges(const std::vector<poly_t::edge_t > &edges,
                                  const poly_t *poly,
                                  const carve::csg::EVSMap &emap,
                                  carve::csg::EVVMap &face_edges) {
  static carve::TimingName FUNC_NAME("CSG::divideEdges()");
  carve::TimingBlock block(FUNC_NAME);

  for (std::vector<poly_t::edge_t >::const_iterator
         i = edges.begin(), e = edges.end();
       i != e;
       ++i) {
    const poly_t::edge_t *edge = (&(*i));
    carve::csg::EVSMap::const_iterator ei = emap.find(edge);
    if (ei != emap.end()) {
      const carve::csg::EVSMap::mapped_type &vertices = ((*ei).second);
      std::vector<const poly_t::vertex_t *> &verts = (face_edges[edge]);
      orderVertices(edge->v2->v - edge->v1->v, edge->v1->v,
                    vertices.begin(), vertices.end(),
                    verts, vertices.size());
    }
  }
}



void carve::csg::CSG::makeFaceEdges(carve::csg::FV2SMap &face_edges,
                                    carve::csg::EdgeClassification &eclass,
                                    carve::csg::FVSMap &fmap,
                                    carve::csg::VFSMap &fmap_rev) {
  carve::csg::FSet face_b_set;
  for (carve::csg::FVSMap::const_iterator
         i = fmap.begin(), ie = fmap.end();
       i != ie;
       ++i) {
    const poly_t::face_t *face_a = (*i).first;
    const carve::csg::FVSMap::mapped_type &face_a_intersections = ((*i).second);

    face_b_set.clear();

    // work out the set of faces from the opposing polyhedron that intersect face_a.
    for (carve::csg::FVSMap::mapped_type::const_iterator
           j = face_a_intersections.begin(), je = face_a_intersections.end();
         j != je;
         ++j) {
      for (carve::csg::VFSMap::mapped_type::const_iterator
             k = fmap_rev[*j].begin(), ke = fmap_rev[*j].end();
           k != ke;
           ++k) {
        const poly_t::face_t *face_b = (*k);
        if (face_a != face_b && face_b->owner != face_a->owner) {
          face_b_set.insert(face_b);
        }
      }
    }

    // run through each intersecting face.
    for (carve::csg::FSet::const_iterator
           j = face_b_set.begin(), je = face_b_set.end();
         j != je;
         ++j) {
      const poly_t::face_t *face_b = (*j);
      const carve::csg::FVSMap::mapped_type &face_b_intersections = (fmap[face_b]);
      carve::csg::VSet vertices;

      // record the points of intersection between face_a and face_b
      for (carve::csg::FVSMap::mapped_type::const_iterator
             k = face_a_intersections.begin(), ke = face_a_intersections.end();
           k != ke;
           ++k) {
        const poly_t::vertex_t *v = *k;
        if (face_b_intersections.find(v) != face_b_intersections.end()) {
          vertices.insert(v);
        }
      }

#if defined(DEBUG)
      std::cerr << "face pair: "
                << face_a << ":" << face_b
                << " N(verts) " << vertices.size() << std::endl;
      for (VSet::const_iterator i = vertices.begin(), e = vertices.end(); i != e; ++i) {
        std::cerr << (*i) << " " << (*i)->v << " ("
                  << carve::geom::distance(face_a->plane_eqn, (*i)->v) << ","
                  << carve::geom::distance(face_b->plane_eqn, (*i)->v) << ")"
                  << std::endl;
        //ASSERT(carve::geom3d::distance(face_a->plane_eqn, *(*i)) < EPSILON);
        //ASSERT(carve::geom3d::distance(face_b->plane_eqn, *(*i)) < EPSILON);
      }
#endif

      // if there are two points of intersection, then the added edge is simple to determine.
      if (vertices.size() == 2) {
        VSet::const_iterator t = vertices.begin();
        const poly_t::vertex_t *v1 = *t++;
        const poly_t::vertex_t *v2 = *t++;
        carve::geom3d::Vector c = (v1->v + v2->v) / 2;

        // determine whether the midpoint of the implied edge is contained in face_a and face_b

#if defined(DEBUG)
        std::cerr << "face_a->vertices.size() = " << face_a->vertices.size() << " face_a->containsPoint(c) = " << face_a->containsPoint(c) << std::endl;
        std::cerr << "face_b->vertices.size() = " << face_b->vertices.size() << " face_b->containsPoint(c) = " << face_b->containsPoint(c) << std::endl;
#endif

        if (face_a->containsPoint(c) && face_b->containsPoint(c)) {
#if defined(DEBUG)
          std::cerr << "adding edge: " << v1 << "-" << v2 << std::endl;
#if defined(DEBUG_DRAW_FACE_EDGES)
          HOOK(drawEdge(v1, v2, 1, 1, 1, 1, 1, 1, 1, 1, 2.0););
#endif
#endif
          // record the edge, with class information.
          if (v1 > v2) std::swap(v1, v2);
          eclass[unoriented_edge(v1, v2)] = EC2(EDGE_ON, EDGE_ON);
          face_edges[face_a].insert(std::make_pair(v1, v2));
          face_edges[face_b].insert(std::make_pair(v1, v2));
        }
        continue;
      }

      // otherwise, it's more complex.
      carve::geom3d::Vector base, dir;
      std::vector<const poly_t::vertex_t *> ordered;

      // skip coplanar edges. this simplifies the resulting
      // mesh. eventually all coplanar face regions of two polyhedra
      // must reach a point where they are no longer coplanar (or the
      // polyhedra are identical).
      if (!facesAreCoplanar(face_a, face_b)) {
        // order the intersection vertices (they must lie along a
        // vector, as the faces aren't coplanar).
        selectOrderingProjection(dir, base, vertices.begin(), vertices.end());
        orderVertices(dir, base, vertices.begin(), vertices.end(), ordered, vertices.size());

        // for each possible edge in the ordering, test the midpoint,
        // and record if it's contained in face_a and face_b.
        for (int k = 0, ke = (int)ordered.size() - 1; k < ke; ++k) {
          const poly_t::vertex_t *v1 = ordered[k];
          const poly_t::vertex_t *v2 = ordered[k + 1];
          carve::geom3d::Vector c = (v1->v + v2->v) / 2;

#if defined(DEBUG)
          std::cerr << "testing edge: " << v1 << "-" << v2 << " at " << c << std::endl;
          std::cerr << "a: " << face_a->containsPoint(c) << " b: " << face_b->containsPoint(c) << std::endl;
          std::cerr << "face_a->containsPoint(c): " << face_a->containsPoint(c) << std::endl;
          std::cerr << "face_b->containsPoint(c): " << face_b->containsPoint(c) << std::endl;
#endif

          if (face_a->containsPoint(c) && face_b->containsPoint(c)) {
#if defined(DEBUG)
            std::cerr << "adding edge: " << v1 << "-" << v2 << std::endl;
#if defined(DEBUG_DRAW_FACE_EDGES)
            HOOK(drawEdge(v1, v2, .5, .5, .5, 1, .5, .5, .5, 1, 2.0););
#endif
#endif
            // record the edge, with class information.
            if (v1 > v2) std::swap(v1, v2);
            eclass[unoriented_edge(v1, v2)] = EC2(EDGE_ON, EDGE_ON);
            face_edges[face_a].insert(std::make_pair(v1, v2));
            face_edges[face_b].insert(std::make_pair(v1, v2));
          }
        }
      }
    }
  }
}



/** 
 * 
 * 
 * @param fll 
 */
static void checkFaceLoopIntegrity(carve::csg::FaceLoopList &fll) {
  static carve::TimingName FUNC_NAME("CSG::checkFaceLoopIntegrity()");
  carve::TimingBlock block(FUNC_NAME);

  std::unordered_map<carve::csg::V2, int, carve::poly::hash_vertex_ptr> counts;
  for (carve::csg::FaceLoop *fl = fll.head; fl; fl = fl->next) {
    std::vector<const poly_t::vertex_t *> &loop = (fl->vertices);
    const poly_t::vertex_t *v1, *v2;
    v1 = loop[loop.size() - 1];
    for (unsigned i = 0; i < loop.size(); ++i) {
      v2 = loop[i];
      if (v1 < v2) {
        counts[std::make_pair(v1, v2)]++;
      } else {
        counts[std::make_pair(v2, v1)]--;
      }
      v1 = v2;
    }
  }
  for (std::unordered_map<carve::csg::V2, int, carve::poly::hash_vertex_ptr>::const_iterator
         x = counts.begin(), xe = counts.end(); x != xe; ++x) {
    if ((*x).second) {
      std::cerr << "FACE LOOP ERROR: " << (*x).first.first << "-" << (*x).first.second << " : " << (*x).second << std::endl;
    }
  }
}



/** 
 * 
 * 
 * @param a 
 * @param b 
 * @param vclass 
 * @param eclass 
 * @param a_face_loops 
 * @param b_face_loops 
 * @param a_edge_count 
 * @param b_edge_count 
 * @param hooks 
 */
void carve::csg::CSG::calc(const poly_t *a,
                           const poly_t *b,
                           carve::csg::VertexClassification &vclass,
                           carve::csg::EdgeClassification &eclass,
                           carve::csg::FaceLoopList &a_face_loops,
                           carve::csg::FaceLoopList &b_face_loops,
                           size_t &a_edge_count,
                           size_t &b_edge_count) {
  carve::csg::VVMap vmap;
  carve::csg::EVSMap emap;
  carve::csg::FV2SMap face_split_edges;
  carve::csg::EVVMap divided_edges;
  carve::csg::FVSMap fmap;
  carve::csg::VFSMap fmap_rev;

#if defined(DEBUG)
  std::cerr << "init" << std::endl;
#endif
  init();

  generateIntersections(a, b);

#if defined(DEBUG)
  std::cerr << "intersectingFacePairs" << std::endl;
#endif
  intersectingFacePairs(vmap, emap, fmap, fmap_rev);

#if defined(DEBUG)
  std::cerr << "emap:" << std::endl;
  map_histogram(emap);
  std::cerr << "fmap:" << std::endl;
  map_histogram(fmap);
  std::cerr << "fmap_rev:" << std::endl;
  map_histogram(fmap_rev);
#endif

  // std::cerr << "removeCoplanarFaces" << std::endl;
  // fp_intersections.removeCoplanarFaces();

#if defined(DEBUG) && defined(DEBUG_DRAW_OCTREE)
  HOOK(drawOctree(a->octree););
  HOOK(drawOctree(b->octree););
#endif

#if defined(DEBUG)
  std::cerr << "divideEdges" << std::endl;
#endif
  divideEdges(a->edges, b, emap, divided_edges);
  divideEdges(b->edges, a, emap, divided_edges);

#if defined(DEBUG)
  std::cerr << "makeFaceEdges" << std::endl;
#endif
  makeFaceEdges(face_split_edges, eclass, fmap, fmap_rev);

#if defined(DEBUG)
  std::cerr << "generateFaceLoops" << std::endl;
#endif
  a_edge_count = generateFaceLoops(a, vmap, face_split_edges, divided_edges, a_face_loops);
  b_edge_count = generateFaceLoops(b, vmap, face_split_edges, divided_edges, b_face_loops);

#if defined(DEBUG)
  std::cerr << "generated " << a_edge_count << " edges for poly a" << std::endl;
  std::cerr << "generated " << b_edge_count << " edges for poly b" << std::endl;
#endif

  checkFaceLoopIntegrity(a_face_loops);
  checkFaceLoopIntegrity(b_face_loops);

#if defined(DEBUG)
  std::cerr << "classify" << std::endl;
#endif
  // initialize some classification information.
  for (std::vector<poly_t::vertex_t>::const_iterator
         i = a->vertices.begin(), e = a->vertices.end(); i != e; ++i) {
    vclass[map_vertex(vmap, &(*i))].cls[0] = POINT_ON;
  }
  for (std::vector<poly_t::vertex_t>::const_iterator
         i = b->vertices.begin(), e = b->vertices.end(); i != e; ++i) {
    vclass[map_vertex(vmap, &(*i))].cls[1] = POINT_ON;
  }
  for (VertexIntersections::const_iterator
         i = vertex_intersections.begin(), e = vertex_intersections.end(); i != e; ++i) {
    vclass[(*i).first] = PC2(POINT_ON, POINT_ON);
  }

#if defined(DEBUG)
  std::cerr << divided_edges.size() << " edges are split" << std::endl;
  std::cerr << face_split_edges.size() << " faces are split" << std::endl;

  std::cerr << "poly a: " << a_face_loops.size() << " face loops" << std::endl;
  std::cerr << "poly b: " << b_face_loops.size() << " face loops" << std::endl;
#endif

  // std::cerr << "OCTREE A:" << std::endl;
  // dump_octree_stats(a->octree.root, 0);
  // std::cerr << "OCTREE B:" << std::endl;
  // dump_octree_stats(b->octree.root, 0);
}



/** 
 * 
 * 
 * @param shared_edges 
 * @param result_list 
 * @param shared_edge_ptr 
 */
void returnSharedEdges(carve::csg::V2Set &shared_edges, 
                       std::list<poly_t *> &result_list,
                       carve::csg::V2Set *shared_edge_ptr) {
  // need to convert shared edges to point into result
  typedef std::map<carve::geom3d::Vector, poly_t::vertex_t *> remap_type;
  remap_type remap;
  for (std::list<poly_t *>::iterator list_it =
         result_list.begin(); list_it != result_list.end(); list_it++) {
    poly_t *result = *list_it;
    if (result) {
      for (std::vector<poly_t::vertex_t>::iterator it =
             result->vertices.begin(); it != result->vertices.end(); it++) {
        remap.insert(std::make_pair((*it).v, &(*it)));
      }
    }
  }
  for (carve::csg::V2Set::iterator it = shared_edges.begin(); 
       it != shared_edges.end(); it++) {
    remap_type::iterator first_it = remap.find(((*it).first)->v);
    remap_type::iterator second_it = remap.find(((*it).second)->v);
    assert(first_it != remap.end() && second_it != remap.end());
    shared_edge_ptr->insert(std::make_pair(first_it->second, second_it->second));
  }
}



/** 
 * 
 * 
 * @param a 
 * @param b 
 * @param collector 
 * @param hooks 
 * @param shared_edges_ptr 
 * @param classify_type 
 * 
 * @return 
 */
poly_t *carve::csg::CSG::compute(const poly_t *a,
                                                  const poly_t *b,
                                                  carve::csg::CSG::Collector &collector,
                                                  carve::csg::V2Set *shared_edges_ptr,
                                                  CLASSIFY_TYPE classify_type) {
  static carve::TimingName FUNC_NAME("CSG::compute");
  carve::TimingBlock block(FUNC_NAME);

  VertexClassification vclass;
  EdgeClassification eclass;

  FLGroupList a_loops_grouped;
  FLGroupList b_loops_grouped;

  FaceLoopList a_face_loops;
  FaceLoopList b_face_loops;

  size_t a_edge_count;
  size_t b_edge_count;

  {
    static carve::TimingName FUNC_NAME("CSG::compute - calc()");
    carve::TimingBlock block(FUNC_NAME);
    calc(a, b, vclass, eclass,a_face_loops, b_face_loops, a_edge_count, b_edge_count);
  }

  LoopEdges a_edge_map;
  LoopEdges b_edge_map;

  {
    static carve::TimingName FUNC_NAME("CSG::compute - makeEdgeMap()");
    carve::TimingBlock block(FUNC_NAME);
    makeEdgeMap(a_face_loops, a_edge_count, a_edge_map);
    makeEdgeMap(b_face_loops, b_edge_count, b_edge_map);

  }
  
  {
    static carve::TimingName FUNC_NAME("CSG::compute - sortFaceLoopLists()");
    carve::TimingBlock block(FUNC_NAME);
    a_edge_map.sortFaceLoopLists();
    b_edge_map.sortFaceLoopLists();
  }

  V2Set shared_edges;
  
  {
    static carve::TimingName FUNC_NAME("CSG::compute - findSharedEdges()");
    carve::TimingBlock block(FUNC_NAME);
    findSharedEdges(a_edge_map, b_edge_map, shared_edges);
  }

  {
    static carve::TimingName FUNC_NAME("CSG::compute - groupFaceLoops()");
    carve::TimingBlock block(FUNC_NAME);
    groupFaceLoops(a_face_loops, a_edge_map, shared_edges, a_loops_grouped);
    groupFaceLoops(b_face_loops, b_edge_map, shared_edges, b_loops_grouped);
#if defined(DEBUG)
    std::cerr << "*** a_loops_grouped.size(): " << a_loops_grouped.size() << std::endl;
    std::cerr << "*** b_loops_grouped.size(): " << b_loops_grouped.size() << std::endl;
#endif
  }

#if defined(DEBUG) && defined(DEBUG_DRAW_GROUPS)
  {
    float n = 1.0 / (a_loops_grouped.size() + b_loops_grouped.size() + 1);
    float H = 0.0, S = 1.0, V = 1.0;
    float r, g, b;
    for (FLGroupList::const_iterator i = a_loops_grouped.begin(); i != a_loops_grouped.end(); ++i) {
      HSV2RGB(H, S, V, r, g, b); H += n;
      drawFaceLoopList((*i).face_loops, r, g, b, 1.0, r * .5, g * .5, b * .5, 1.0, true);
    }
    for (FLGroupList::const_iterator i = b_loops_grouped.begin(); i != b_loops_grouped.end(); ++i) {
      HSV2RGB(H, S, V, r, g, b); H += n;
      drawFaceLoopList((*i).face_loops, r, g, b, 1.0, r * .5, g * .5, b * .5, 1.0, true);
    }

    for (FLGroupList::const_iterator i = a_loops_grouped.begin(); i != a_loops_grouped.end(); ++i) {
      drawFaceLoopListWireframe((*i).face_loops);
    }
    for (FLGroupList::const_iterator i = b_loops_grouped.begin(); i != b_loops_grouped.end(); ++i) {
      drawFaceLoopListWireframe((*i).face_loops);
    }
  }

#endif

  switch (classify_type) {
  case CLASSIFY_EDGE:
    classifyFaceGroupsEdge(shared_edges,
                           vclass,
                           a,
                           a_loops_grouped,
                           a_edge_map,
                           b,
                           b_loops_grouped,
                           b_edge_map,
                           collector);
    break;
  case CLASSIFY_NORMAL:
    classifyFaceGroups(shared_edges,
                       vclass,
                       a,
                       a_loops_grouped,
                       a_edge_map,
                       b,
                       b_loops_grouped,
                       b_edge_map,
                       collector);
    break;
  }

  poly_t *result = collector.done(hooks);
  if (result != NULL && shared_edges_ptr != NULL) {
    std::list<poly_t *> result_list;
    result_list.push_back(result);
    returnSharedEdges(shared_edges, result_list, shared_edges_ptr);
  }
  return result;
}



/** 
 * 
 * 
 * @param a 
 * @param b 
 * @param op 
 * @param hooks 
 * @param shared_edges 
 * @param classify_type 
 * 
 * @return 
 */
poly_t *carve::csg::CSG::compute(const poly_t *a,
                                                  const poly_t *b,
                                                  carve::csg::CSG::OP op,
                                                  carve::csg::V2Set *shared_edges,
                                                  CLASSIFY_TYPE classify_type) {
  Collector *coll = makeCollector(op, a, b);
  if (!coll) return NULL;

  poly_t *result = compute(a, b, *coll, shared_edges, classify_type);
     
  delete coll;

  return result;
}



/** 
 * 
 * 
 * @param closed 
 * @param open 
 * @param FaceClass 
 * @param result 
 * @param hooks 
 * @param shared_edges_ptr 
 * 
 * @return 
 */
bool carve::csg::CSG::sliceAndClassify(const poly_t *closed,
                                       const poly_t *open,
                                       std::list<std::pair<FaceClass, poly_t *> > &result,
                                       carve::csg::V2Set *shared_edges_ptr) {
  if (closed->hasOpenManifolds()) return false;
  carve::csg::VertexClassification vclass;
  carve::csg::EdgeClassification eclass;

  carve::csg::FLGroupList a_loops_grouped;
  carve::csg::FLGroupList b_loops_grouped;

  carve::csg::FaceLoopList a_face_loops;
  carve::csg::FaceLoopList b_face_loops;

  size_t a_edge_count;
  size_t b_edge_count;

  calc(closed, open, vclass, eclass,a_face_loops, b_face_loops, a_edge_count, b_edge_count);

  carve::csg::LoopEdges a_edge_map;
  carve::csg::LoopEdges b_edge_map;

  makeEdgeMap(a_face_loops, a_edge_count, a_edge_map);
  makeEdgeMap(b_face_loops, b_edge_count, b_edge_map);
  
  carve::csg::V2Set shared_edges;
  
  findSharedEdges(a_edge_map, b_edge_map, shared_edges);
  
  groupFaceLoops(a_face_loops, a_edge_map, shared_edges, a_loops_grouped);
  groupFaceLoops(b_face_loops, b_edge_map, shared_edges, b_loops_grouped);

  halfClassifyFaceGroups(shared_edges,
                         vclass,
                         closed,
                         a_loops_grouped,
                         a_edge_map,
                         open,
                         b_loops_grouped,
                         b_edge_map,
                         result);

  if (shared_edges_ptr != NULL) {
    std::list<poly_t *> result_list;
    for (std::list<std::pair<FaceClass, poly_t *> >::iterator it = result.begin(); it != result.end(); it++) {
      result_list.push_back(it->second);
    }
    returnSharedEdges(shared_edges, result_list, shared_edges_ptr);
  }
  return true;
}



/** 
 * 
 * 
 * @param a 
 * @param b 
 * @param a_sliced 
 * @param b_sliced 
 * @param hooks 
 * @param shared_edges_ptr 
 */
void carve::csg::CSG::slice(const poly_t *a,
                            const poly_t *b,
                            std::list<poly_t *> &a_sliced,
                            std::list<poly_t *> &b_sliced,
                            carve::csg::V2Set *shared_edges_ptr) {
  carve::csg::VertexClassification vclass;
  carve::csg::EdgeClassification eclass;

  carve::csg::FLGroupList a_loops_grouped;
  carve::csg::FLGroupList b_loops_grouped;

  carve::csg::FaceLoopList a_face_loops;
  carve::csg::FaceLoopList b_face_loops;

  size_t a_edge_count;
  size_t b_edge_count;

  calc(a, b, vclass, eclass,a_face_loops, b_face_loops, a_edge_count, b_edge_count);

  carve::csg::LoopEdges a_edge_map;
  carve::csg::LoopEdges b_edge_map;
      
  makeEdgeMap(a_face_loops, a_edge_count, a_edge_map);
  makeEdgeMap(b_face_loops, b_edge_count, b_edge_map);
  
  carve::csg::V2Set shared_edges;
  
  findSharedEdges(a_edge_map, b_edge_map, shared_edges);
  
  groupFaceLoops(a_face_loops, a_edge_map, shared_edges, a_loops_grouped);
  groupFaceLoops(b_face_loops, b_edge_map, shared_edges, b_loops_grouped);

  for (carve::csg::FLGroupList::iterator
         i = a_loops_grouped.begin(), e = a_loops_grouped.end();
       i != e; ++i) {
    carve::csg::FaceLoopList &fl = ((*i).face_loops);
    std::vector<poly_t::face_t > faces;
    faces.reserve(fl.size());
    for (FaceLoop *f = fl.head; f; f = f->next) {
      faces.push_back(poly_t::face_t());
      faces.back().init(f->orig_face, f->vertices, false);
    }
    a_sliced.push_back(new poly_t(faces));
  }

  for (carve::csg::FLGroupList::iterator
         i = b_loops_grouped.begin(), e = b_loops_grouped.end();
       i != e; ++i) {
    carve::csg::FaceLoopList &fl = ((*i).face_loops);
    std::vector<poly_t::face_t > faces;
    faces.reserve(fl.size());
    for (FaceLoop *f = fl.head; f; f = f->next) {
      faces.push_back(poly_t::face_t());
      faces.back().init(f->orig_face, f->vertices, false);
    }
    b_sliced.push_back(new poly_t(faces));
  }
  if (shared_edges_ptr != NULL) {
    std::list<poly_t *> result_list;
    result_list.insert(result_list.end(), a_sliced.begin(), a_sliced.end());
    result_list.insert(result_list.end(), b_sliced.begin(), b_sliced.end());
    returnSharedEdges(shared_edges, result_list, shared_edges_ptr);
  }
}



/** 
 * 
 * 
 */
void carve::csg::CSG::init() {
  intersections.clear();
  vertex_intersections.clear();
  vertex_pool.reset();
}
