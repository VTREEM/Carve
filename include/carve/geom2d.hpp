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

#include <carve/carve.hpp>

#include <carve/math.hpp>
#include <carve/math_constants.hpp>

#include <carve/geom.hpp>

#include <vector>
#include <iostream>

#include <math.h>

namespace carve {
  namespace geom2d {

    typedef carve::geom::vector<2> P2;
    typedef carve::geom::ray<2> Ray2;
    typedef carve::geom::linesegment<2> LineSegment2;

    struct p2_adapt_ident {
      P2 &operator()(P2 &p) { return p; }
      const P2 &operator()(const P2 &p) { return p; }
    };

    typedef std::vector<P2> P2Vector;

    static inline double atan2(const P2 &p) {
      return ::atan2(p.y, p.x);
    }

    struct LineIntersectionInfo {
      LineIntersectionClass iclass;
      P2 ipoint;
      int p1, p2;

      LineIntersectionInfo(LineIntersectionClass _iclass,
                           P2 _ipoint = P2::ZERO(),
                           int _p1 = -1,
                           int _p2 = -1) :
        iclass(_iclass), ipoint(_ipoint), p1(_p1), p2(_p2) {
      }
    };

    struct PolyInclusionInfo {
      PointClass iclass;
      int iobjnum;

      PolyInclusionInfo(PointClass _iclass,
                        int _iobjnum = -1) :
        iclass(_iclass), iobjnum(_iobjnum) {
      }
    };

    struct PolyIntersectionInfo {
      IntersectionClass iclass;
      P2 ipoint;
      size_t iobjnum;

      PolyIntersectionInfo(IntersectionClass _iclass,
                           const P2 &_ipoint,
                           size_t _iobjnum) :
        iclass(_iclass), ipoint(_ipoint), iobjnum(_iobjnum) {
      }
    };

    LineIntersectionInfo lineSegmentIntersection(const P2 &l1v1, const P2 &l1v2, const P2 &l2v1, const P2 &l2v2);
    LineIntersectionInfo lineSegmentIntersection(const LineSegment2 &l1, const LineSegment2 &l2);

    double signedArea(const std::vector<P2> &points);

    bool pointInPolySimple(const std::vector<P2> &points, const P2 &p);

    PolyInclusionInfo pointInPoly(const std::vector<P2> &points, const P2 &p);

    int lineSegmentPolyIntersections(const std::vector<P2> &points,
                                     LineSegment2 line,
                                     std::vector<PolyInclusionInfo> &out);

    int sortedLineSegmentPolyIntersections(const std::vector<P2> &points,
                                           LineSegment2 line,
                                           std::vector<PolyInclusionInfo> &out);

    bool pickContainedPoint(const std::vector<P2> &poly, P2 &result);

    template<typename T, typename adapt_t>
    double signedArea(const std::vector<T> &points, adapt_t adapt) {
      P2Vector::size_type l = points.size();
      double A = 0.0;

      for (P2Vector::size_type i = 0; i < l - 1; i++) {
        A += (adapt(points[i + 1]).y + adapt(points[i]).y) * (adapt(points[i + 1]).x - adapt(points[i]).x);
      }
      A += (adapt(points[0]).y + adapt(points[l - 1]).y) * (adapt(points[0]).x - adapt(points[l - 1]).x);

      return A / 2.0;
    }



    template<typename T, typename adapt_t>
    bool pointInPolySimple(const std::vector<T> &points, adapt_t adapt, const P2 &p) {
      ASSERT(points.size() > 0);
      P2Vector::size_type l = points.size();
      double s = 0.0;
      double rp, r0, d;

      rp = r0 = atan2(adapt(points[0]) - p);

      for (P2Vector::size_type i = 1; i < l; i++) {
        double r = atan2(adapt(points[i]) - p);
        d = r - rp;
        if (d > M_PI) d -= M_TWOPI;
        if (d < -M_PI) d += M_TWOPI;
        s = s + d;
        rp = r;
      }

      d = r0 - rp;
      if (d > M_PI) d -= M_TWOPI;
      if (d < -M_PI) d += M_TWOPI;
      s = s + d;

      return !carve::math::ZERO(s);
    }



    template<typename T, typename adapt_t>
    PolyInclusionInfo pointInPoly(const std::vector<T> &points, adapt_t adapt, const P2 &p) {
      P2Vector::size_type l = points.size();
      for (unsigned i = 0; i < l; i++) {
        if (equal(adapt(points[i]), p)) return PolyInclusionInfo(POINT_VERTEX, i);
      }

      for (unsigned i = 0; i < l; i++) {
        unsigned j = (i + 1) % l;

        if (std::min(adapt(points[i]).x, adapt(points[j]).x) - EPSILON < p.x &&
            std::max(adapt(points[i]).x, adapt(points[j]).x) + EPSILON > p.x &&
            std::min(adapt(points[i]).y, adapt(points[j]).y) - EPSILON < p.y &&
            std::max(adapt(points[i]).y, adapt(points[j]).y) + EPSILON > p.y &&
            distance2(carve::geom::rayThrough(adapt(points[i]), adapt(points[j])), p) < EPSILON2) {
          return PolyInclusionInfo(POINT_EDGE, i);
        }
      }

      if (pointInPolySimple(points, adapt, p)) {
        return PolyInclusionInfo(POINT_IN);
      }

      return PolyInclusionInfo(POINT_OUT);
    }



    template<typename T, typename adapt_t>
    bool pickContainedPoint(const std::vector<T> &poly, adapt_t adapt, P2 &result) {
#if defined(DEBUG)
      std::cerr << "pickContainedPoint ";
      for (unsigned i = 0; i < poly.size(); ++i) std::cerr << " " << adapt(poly[i]);
      std::cerr << std::endl;
#endif

      const size_t S = poly.size();
      P2 a, b, c;
      for (unsigned i = 0; i < S; ++i) {
        a = adapt(poly[i]);
        b = adapt(poly[(i + 1) % S]);
        c = adapt(poly[(i + 2) % S]);

        if (cross(a - b, c - b) < 0) {
          P2 p = (a + b + c) / 3;
          if (pointInPolySimple(poly, adapt, p)) {
            result = p;
            return true;
          }
        }
      }
      return false;
    }

    /** 
     * \brief Return the orientation of c with respect to the ray defined by a->b.
     *
     * (Can be implemented exactly)
     * 
     * @param[in] a 
     * @param[in] b 
     * @param[in] c 
     * 
     * @return positive, if c to the left of a->b.
     *         zero, if c is colinear with a->b.
     *         negative, if c to the right of a->b.
     */
    inline double orient2d(const carve::geom2d::P2 &a, const carve::geom2d::P2 &b, const carve::geom2d::P2 &c) {
      double acx = a.x - c.x;
      double bcx = b.x - c.x;
      double acy = a.y - c.y;
      double bcy = b.y - c.y;
      return acx * bcy - acy * bcx;
    }

  }
}
