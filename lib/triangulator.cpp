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
#include <carve/triangulator.hpp>

#include <algorithm>

namespace {

#if defined(DEBUG)
      void dumpPoly(const std::vector<carve::geom2d::P2> &points,
                    const std::vector<carve::triangulate::tri_idx> &result);
#endif 
  struct vertex_info;

  double ear_angle(const carve::geom2d::P2 &prev, const carve::geom2d::P2 &curr, const carve::geom2d::P2 &next) {
    double a = carve::geom2d::atan2(prev - curr) - carve::geom2d::atan2(next - curr);
    if (a < 0) a += M_PI * 2;
    return a;
  }

  bool isLeft(const vertex_info *a,
               const vertex_info *b,
               const vertex_info *c);

  // 52% of execution time.
  bool pointInTriangle(const vertex_info *a,
                       const vertex_info *b,
                       const vertex_info *c,
                       const vertex_info *d) {
    return !isLeft(a, c, d) && !isLeft(b, a, d) && !isLeft(c, b, d);
  }

  struct vertex_info {
    vertex_info *prev;
    vertex_info *next;
    const carve::geom2d::P2 &p;
    size_t idx;
    double score;
    bool convex;
    bool failed;

    vertex_info(const carve::geom2d::P2 &_p, size_t _idx) :
        prev(NULL), next(NULL),
        p(_p), idx(_idx),
        score(0.0), convex(false) {
    }

    static double triScore(const vertex_info *p, const vertex_info *v, const vertex_info *n) {
#if 0
      bool convex = isLeft(p, v, n);
      if (!convex) return -1e-5;

      double a1 = carve::geom2d::atan2(p->p - v->p) - carve::geom2d::atan2(n->p - v->p);
      double a2 = carve::geom2d::atan2(v->p - n->p) - carve::geom2d::atan2(p->p - n->p);
      if (a1 < 0) a1 += M_PI * 2;
      if (a2 < 0) a2 += M_PI * 2;

      std::cerr << "a1: " << a1 << " a2: " << a2 << " a3: " << (M_PI - a1 - a2) << std::endl;
      std::cerr << "...score: " << std::min(a1, std::min(a2, M_PI - a1 - a2)) / (M_PI / 3) << std::endl;
      return std::min(a1, std::min(a2, M_PI - a1 - a2)) / (M_PI / 3);
#endif

#if 1
      // range: 0 - 1
      double a, b, c;
      bool convex = isLeft(p, v, n);
      if (!convex) return -1e-5;
      a = (n->p - v->p).length();
      b = (p->p - n->p).length();
      c = (v->p - p->p).length();
      if (a < 1e-10 || b < 1e-10 || c < 1e-10) return 0.0;
      return std::max(std::min((a+b)/c, std::min((a+c)/b, (b+c)/a)) - 1.0, 0.0);
#endif
    }

    double calcScore() const {
#if 0
      double this_tri = triScore(prev, this, next);
      return this_tri;
#endif

#if 1
      double this_tri = triScore(prev, this, next);
      double next_tri = triScore(prev, next, next->next);
      double prev_tri = triScore(prev->prev, prev, next);
      double next_delta = next_tri - next->score;
      double prev_delta = prev_tri - prev->score;
#if defined(DEBUG)
      std::cerr << "calcScore: " << this_tri << ", " << prev_delta << ", " << next_delta << std::endl;
#endif
      return this_tri + std::max(next_tri, prev_tri) * .5;
#endif

#if 0
      double score = triScore(prev, this, next);

      // penalise ears that will require producing a sliver triangle.
      double a1, a2;
      a1 = carve::geom2d::atan2(prev->p - next->p);
      a2 = carve::geom2d::atan2(next->next->p - next->p);
      if (fabs(a1 - a2) < 1e-5) score -= .5;

      a1 = carve::geom2d::atan2(next->p - prev->p);
      a2 = carve::geom2d::atan2(prev->prev->p - prev->p);
      if (fabs(a1 - a2) < 1e-5) score -= .5;

      return score;
#endif
    }

    void recompute() {
      score = calcScore();
      convex = isLeft(prev, this, next);
      failed = false;
    }

    bool isCandidate() const {
      return convex && !failed;
    }

    void remove() {
      next->prev = prev;
      prev->next = next;
    }
  };

  struct vertex_info_ordering {
    bool operator()(const vertex_info *a, const vertex_info *b) const {
      return a->score < b->score;
    }
  };

  struct vertex_info_l2norm_inc_ordering {
    const vertex_info *v;
    vertex_info_l2norm_inc_ordering(const vertex_info *_v) : v(_v) {
    }
    bool operator()(const vertex_info *a, const vertex_info *b) const {
      return carve::geom::distance2(v->p, a->p) > carve::geom::distance2(v->p, b->p);
    }
  };

  bool isLeft(const vertex_info *a,
              const vertex_info *b,
              const vertex_info *c) {
    if (a->idx < b->idx && b->idx < c->idx) {
      return carve::geom2d::orient2d(a->p, b->p, c->p) > 0.0;
    } else if (a->idx < c->idx && c->idx < b->idx) {
      return carve::geom2d::orient2d(a->p, c->p, b->p) < 0.0;
    } else if (b->idx < a->idx && a->idx < c->idx) {
      return carve::geom2d::orient2d(b->p, a->p, c->p) < 0.0;
    } else if (b->idx < c->idx && c->idx < a->idx) {
      return carve::geom2d::orient2d(b->p, c->p, a->p) > 0.0;
    } else if (c->idx < a->idx && a->idx < b->idx) {
      return carve::geom2d::orient2d(c->p, a->p, b->p) > 0.0;
    } else {
      return carve::geom2d::orient2d(c->p, b->p, a->p) < 0.0;
    }
  }

  // =====================================================================
  class EarQueue {
    std::vector<vertex_info *> queue;

    void checkheap() {
#ifdef __GNUC__
      ASSERT(std::__is_heap(queue.begin(), queue.end(), vertex_info_ordering()));
#endif
    }

  public:
    EarQueue() {
    }

    size_t size() const {
      return queue.size();
    }

    void push(vertex_info *v) {
      checkheap();
      queue.push_back(v);
      std::push_heap(queue.begin(), queue.end(), vertex_info_ordering());
    }

    vertex_info *pop() {
      checkheap();
      std::pop_heap(queue.begin(), queue.end(), vertex_info_ordering());
      vertex_info *v = queue.back();
      queue.pop_back();
      return v;
    }

    void remove(vertex_info *v) {
      checkheap();
      ASSERT(std::find(queue.begin(), queue.end(), v) != queue.end());
      double score = v->score;
      if (v != queue[0]) {
        v->score = queue[0]->score + 1;
        std::make_heap(queue.begin(), queue.end(), vertex_info_ordering());
      }
      ASSERT(v == queue[0]);
      std::pop_heap(queue.begin(), queue.end(), vertex_info_ordering());
      ASSERT(queue.back() == v);
      queue.pop_back();
      v->score = score;
    }

    void changeScore(vertex_info *v, double score) {
      checkheap();
      ASSERT(std::find(queue.begin(), queue.end(), v) != queue.end());
      if (v->score != score) {
        v->score = score;
        std::make_heap(queue.begin(), queue.end(), vertex_info_ordering());
      }
    }
  };




  bool inCone(const vertex_info *a,
              const vertex_info *b,
              const vertex_info *c,
              const carve::geom2d::P2 &p) {
    return
      carve::geom2d::orient2d(a->p, b->p, p) > 0.0 &&
      carve::geom2d::orient2d(b->p, c->p, p) > 0.0;
  }

  // 39% of execution time
  void updateVertex(vertex_info *v, EarQueue &vq) {
    double spre = v->score;
    bool qpre = v->isCandidate();
    v->recompute();
    bool qpost = v->isCandidate();
    double spost = v->score;

    v->score = spre;

    if (qpre) {
      if (qpost) {
        if (v->score != spre) {
          vq.changeScore(v, spost);
        }
      } else {
        vq.remove(v);
      }
    } else {
      if (qpost) {
        vq.push(v);
      }
    }
  }

  // 60% of execution time
  bool isClipable(vertex_info *v) {
    for (vertex_info *v_test = v->next->next; v_test != v->prev; v_test = v_test->next) {
      if (!v_test->convex &&
          v_test->p != v->prev->p &&
          v_test->p != v->next->p &&
          pointInTriangle(v->prev, v, v->next, v_test)) {
        return false;
      }
    }
    return true;
  }

  size_t removeDegeneracies(vertex_info *&begin, std::vector<carve::triangulate::tri_idx> &result) {
    vertex_info *v = begin;
    vertex_info *n;
    size_t count = 0;
    do {
      if (v->p == v->next->p) {
        result.push_back(carve::triangulate::tri_idx(v->idx, v->next->idx, v->next->next->idx));
        n = v->next;
        if (n == begin) begin = n->next;
        n->remove();
        count++;
        delete n;
      } else if (v->p == v->next->next->p) {
        if (v->next->p == v->next->next->p ||
            (!v->convex && !v->next->next->convex)) {
          result.push_back(carve::triangulate::tri_idx(v->idx, v->next->idx, v->next->next->idx));
          // leave v where it is, so that the next degeneracy gets clipped.
          v->next->remove();
          count++;
        }
      } else {
        v = v->next;
      }
    } while (v != begin);
    return count;
  }

  int windingNumber(vertex_info *begin, const carve::geom2d::P2 &point) {
    int wn = 0;

    vertex_info *v = begin;
    do {
      if (v->p.y <= point.y) {
        if (v->next->p.y > point.y && carve::geom2d::orient2d(v->p, v->next->p, point) > 0.0) {
          ++wn;
        }
      } else {
        if (v->next->p.y <= point.y && carve::geom2d::orient2d(v->p, v->next->p, point) < 0.0) {
          --wn;
        }
      }
      v = v->next;
    } while (v != begin);

    return wn;
  }

  bool findDiagonal(vertex_info *begin, vertex_info *&v1, vertex_info *&v2) {
    vertex_info *t;
    std::vector<vertex_info *> heap;

    v1 = begin;
    do {
      heap.clear();

      for (v2 = v1->next->next; v2 != v1->prev; v2 = v2->next) {
        if (!inCone(v1->prev, v1, v1->next, v2->p) ||
            !inCone(v2->prev, v2, v2->next, v1->p)) continue;

        heap.push_back(v2);
        std::push_heap(heap.begin(), heap.end(), vertex_info_l2norm_inc_ordering(v1));
      }

      while (heap.size()) {
        std::pop_heap(heap.begin(), heap.end(), vertex_info_l2norm_inc_ordering(v1));
        v2 = heap.back(); heap.pop_back();

#if defined(DEBUG)
        std::cerr << "testing: " << v1 << " - " << v2 << std::endl;
        std::cerr << "  length = " << (v2->p - v1->p).length() << std::endl;
        std::cerr << "  pos: " << v1->p << " - " << v2->p << std::endl;
#endif
        // test whether v1-v2 is a valid diagonal.
        double v_min_x = std::min(v1->p.x, v2->p.x);
        double v_max_x = std::max(v1->p.x, v2->p.x);

        bool intersected = false;

        for (t = v1->next; !intersected && t != v1->prev; t = t->next) {
          vertex_info *u = t->next;
          if (t == v2 || u == v2) continue;

          double l1 = carve::geom2d::orient2d(v1->p, v2->p, t->p);
          double l2 = carve::geom2d::orient2d(v1->p, v2->p, u->p);

          if ((l1 > 0.0 && l2 > 0.0) || (l1 < 0.0 && l2 < 0.0)) {
            // both on the same side; no intersection
            continue;
          }

          double dx13 = v1->p.x - t->p.x;
          double dy13 = v1->p.y - t->p.y;
          double dx43 = u->p.x - t->p.x;
          double dy43 = u->p.y - t->p.y;
          double dx21 = v2->p.x - v1->p.x;
          double dy21 = v2->p.y - v1->p.y;
          double ua_n = dx43 * dy13 - dy43 * dx13;
          double ub_n = dx21 * dy13 - dy21 * dx13;
          double u_d  = dy43 * dx21 - dx43 * dy21;

          if (carve::math::ZERO(u_d)) {
            // parallel
            if (carve::math::ZERO(ua_n)) {
              // colinear
              if (std::max(t->p.x, u->p.x) >= v_min_x && std::min(t->p.x, u->p.x) <= v_max_x) {
                // colinear and intersecting
                intersected = true;
              }
            }
          } else {
            // not parallel
            double ua = ua_n / u_d;
            double ub = ub_n / u_d;

            if (0.0 <= ua && ua <= 1.0 && 0.0 <= ub && ub <= 1.0) {
              intersected = true;
            }
          }
#if defined(DEBUG)
          if (intersected) {
            std::cerr << "  failed on edge: " << t << " - " << u << std::endl;
            std::cerr << "    pos: " << t->p << " - " << u->p << std::endl;
          }
#endif
        }

        if (!intersected) {
          // test whether midpoint winding == 1

          carve::geom2d::P2 mid = (v1->p + v2->p) / 2;
          if (windingNumber(begin, mid) == 1) {
            // this diagonal is ok
            return true;
          }
        }
      }

      // couldn't find a diagonal from v1 that was ok.
      v1 = v1->next;
    } while (v1 != begin);
    return false;
  }

  bool doTriangulate(vertex_info *begin, std::vector<carve::triangulate::tri_idx> &result);

  bool splitAndResume(vertex_info *begin, std::vector<carve::triangulate::tri_idx> &result) {
    vertex_info *v1, *v2;

    if (!findDiagonal(begin, v1, v2)) return false;

    vertex_info *v1_copy = new vertex_info(*v1);
    vertex_info *v2_copy = new vertex_info(*v2);

    v1->next = v2;
    v2->prev = v1;

    v1_copy->next->prev = v1_copy;
    v2_copy->prev->next = v2_copy;

    v1_copy->prev = v2_copy;
    v2_copy->next = v1_copy;

    bool r1 = doTriangulate(v1, result);
    bool r2 =  doTriangulate(v1_copy, result);
    return r1 && r2;
  }

  bool doTriangulate(vertex_info *begin, std::vector<carve::triangulate::tri_idx> &result) {
    EarQueue vq;

    vertex_info *v = begin;
    size_t remain = 0;
    do {
      if (v->isCandidate()) vq.push(v);
      v = v->next;
      remain++;
    } while (v != begin);

    std::cerr << "doTriangulate; remain=" << remain << std::endl;

    while (vq.size()) {
      vertex_info *v = vq.pop();
      if (!isClipable(v)) {
        v->failed = true;
        continue;
      }

    continue_clipping:
      vertex_info *n = v->next;
      vertex_info *p = v->prev;

      result.push_back(carve::triangulate::tri_idx(v->prev->idx, v->idx, v->next->idx));

#if defined(DEBUG)
      {
        std::vector<carve::geom2d::P2> temp;
        temp.push_back(v->prev->p);
        temp.push_back(v->p);
        temp.push_back(v->next->p);
        std::cerr << "clip " << v << " " << v->idx << " area = " << carve::geom2d::signedArea(temp) << std::endl;
      }
#endif
      v->remove();
      remain--;
      if (v == begin) begin = v->next;
      delete v;

      updateVertex(n, vq);
      updateVertex(p, vq);

      bool swapped = false;
      std::cerr << "  continuation options: next; ";
      if (n->isCandidate()) std::cerr << n->score << " "; else std::cerr << "---";
      std::cerr << " prev; ";
      if (p->isCandidate()) std::cerr << p->score << " "; else std::cerr << "---";
      std::cerr << std::endl;
      if (n->score < p->score) { std::swap(n, p); swapped = true; }

      if (n->score > 0.25 && n->isCandidate() && isClipable(n)) {
        std::cerr << "continue, " << (swapped ? "prev" : "next") << std::endl;
        vq.remove(n);
        v = n;
        goto continue_clipping;
      }

      if (p->score > 0.25 && p->isCandidate() && isClipable(p)) {
        std::cerr << "continue, " << (swapped ? "next" : "prev") << std::endl;
        vq.remove(p);
        v = p;
        goto continue_clipping;
      }

      std::cerr << "looking for new start point" << std::endl;
#if defined(DEBUG)
    {
      std::cerr << "remain = " << remain << std::endl;
      std::vector<carve::triangulate::tri_idx> dummy;
      std::vector<carve::geom2d::P2> dummy_p;
      vertex_info *v = begin;
      do {
        dummy_p.push_back(v->p);
        v = v->next;
      } while (v != begin);
      dumpPoly(dummy_p, dummy);
    }
#endif
    }

    std::cerr << "doTriangulate complete; remain=" << remain << std::endl;

    if (remain < 3) {
      return true;
    }


    if (remain > 3) {
      remain -= removeDegeneracies(begin, result);
    }
    if (remain == 3) {
      result.push_back(carve::triangulate::tri_idx(begin->idx, begin->next->idx, begin->next->next->idx));
      return true;
    } else if (remain > 3) {
      // must split the remainder and recurse.
      if (splitAndResume(begin, result)) return true;
    }

    return false;
  }

#if defined(DEBUG)
  void dumpPoly(const std::vector<carve::geom2d::P2> &points,
                const std::vector<carve::triangulate::tri_idx> &result) {
    double minx = points[0].x, maxx = points[0].x;
    double miny = points[0].y, maxy = points[0].y;

    for (size_t i = 1; i < points.size(); ++i) {
      minx = std::min(points[i].x, minx); maxx = std::max(points[i].x, maxx);
      miny = std::min(points[i].y, miny); maxy = std::max(points[i].y, maxy);
    }

    double width = maxx - minx + 10;
    double height = maxy - miny + 10;
  
    std::cerr << "\
<?xml version=\"1.0\"?>\n\
<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\
<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"" << width << 
"\" height=\"" << height << "\">\n \
";

    std::cerr << "<polygon fill=\"rgb(0,0,0)\" stroke=\"blue\" stroke-width=\"1\" points=\"";
    for (size_t i = 0; i < points.size(); ++i) {
      if (i) std::cerr << ' ';
      double x, y;
      x = points[i].x - minx + 5;
      y = points[i].y - miny + 5;
      std::cerr << x << ',' << y;
    }
    std::cerr << "\" />" << std::endl;

    for (size_t i = 0; i < result.size(); ++i) {
      std::cerr << "<polygon fill=\"rgb(255,255,255)\" stroke=\"black\" stroke-width=\"1\" points=\"";
      double x, y;
      x = points[result[i].a].x - minx + 5;
      y = points[result[i].a].y - miny + 5;
      std::cerr << x << ',' << y << ' ';
      x = points[result[i].b].x - minx + 5;
      y = points[result[i].b].y - miny + 5;
      std::cerr << x << ',' << y << ' ';
      x = points[result[i].c].x - minx + 5;
      y = points[result[i].c].y - miny + 5;
      std::cerr << x << ',' << y;
      std::cerr << "\" />" << std::endl;
    }

    std::cerr << "</svg>" << std::endl;
  }
#endif
}

void carve::triangulate::triangulate(const std::vector<carve::geom2d::P2> &poly,
                                     std::vector<carve::triangulate::tri_idx> &result) {
  std::vector<vertex_info *> vinfo;
  const size_t N = poly.size();

#if defined(DEBUG)
  dumpPoly(poly, result);
  std::cerr << "TRIANGULATION BEGINS" << std::endl;
#endif
  result.clear();
  if (N < 3) {
    return;
  }

  result.reserve(poly.size() - 2);

  if (N == 3) {
    result.push_back(tri_idx(0, 1, 2));
    return;
  }

  vinfo.resize(N);

  vinfo[0] = new vertex_info(poly[0], 0);
  for (size_t i = 1; i < N-1; ++i) {
    vinfo[i] = new vertex_info(poly[i], i);
    vinfo[i]->prev = vinfo[i-1];
    vinfo[i-1]->next = vinfo[i];
  }
  vinfo[N-1] = new vertex_info(poly[N-1], N-1);
  vinfo[N-1]->prev = vinfo[N-2];
  vinfo[N-1]->next = vinfo[0];
  vinfo[0]->prev = vinfo[N-1];
  vinfo[N-2]->next = vinfo[N-1];

  for (size_t i = 0; i < N; ++i) {
    vinfo[i]->recompute();
  }

  vertex_info *begin = vinfo[0];

  //srandom(11);
  //for (int i = 0; i < random() % vinfo.size(); ++i) begin = begin->next;

  removeDegeneracies(begin, result);
  doTriangulate(begin, result);

#if defined(DEBUG)
  dumpPoly(poly, result);
  std::cerr << "TRIANGULATION ENDS" << std::endl;
#endif
}
