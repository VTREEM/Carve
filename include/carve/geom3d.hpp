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
#include <carve/geom.hpp>

#include <math.h>
#include <carve/math_constants.hpp>

#include <vector>
#include <list>
#include <map>

namespace carve {
  namespace geom3d {

    typedef carve::geom::plane<3> Plane;
    typedef carve::geom::ray<3> Ray;
    typedef carve::geom::linesegment<3> LineSegment;
    typedef carve::geom::vector<3> Vector;

    template<typename iter_t, typename adapt_t>
    bool fitPlane(iter_t begin, iter_t end, adapt_t adapt, Plane &plane) {
      Vector centroid;
      carve::geom::centroid(begin, end, adapt, centroid);
      iter_t i;

      Vector n = Vector::ZERO();
      Vector v, z;
      Vector p1, p2, p3, c1, c2;
      if (begin == end) return false;

      i = begin;
      p1 = c1 = adapt(*i++); if (i == end) return false;
      p2 = c2 = adapt(*i++); if (i == end) return false;

#if defined(CARVE_DEBUG)
      size_t N = 2;
#endif
      while (i != end) {
        p3 = adapt(*i++);
        v = cross(p3 - p2, p1 - p2);
        if (v.v[dominantAxis(v)]) v.negate();
        n += v;
        p1 = p2; p2 = p3;
#if defined(CARVE_DEBUG)
        ++N;
#endif
      }

      p1 = p2; p2 = p3; p3 = c1;
      v = cross(p3 - p2, p1 - p2);
      if (v.v[dominantAxis(v)]) v.negate();
      n += v;

      p1 = p2; p2 = p3; p3 = c2;
      v = cross(p3 - p2, p1 - p2);
      if (v.v[dominantAxis(v)]) v.negate();
      n += v;

      n.normalize();
      plane.N = n;
      plane.d = -dot(n, centroid);
#if defined(CARVE_DEBUG)
      if (N > 3) {
        std::cerr << "N = " << N << " fitted distance:";
        for (i = begin; i != end; ++i) {
          Vector p = adapt(*i);
          std::cerr << " {" << p << "} " << distance(plane, p);
        }
        std::cerr << std::endl;
      }
#endif
      return true;
    }

    bool planeIntersection(const Plane &a, const Plane &b, Ray &r);

    IntersectionClass rayPlaneIntersection(const Plane &p,
                                           const Vector &v1,
                                           const Vector &v2,
                                           Vector &v,
                                           double &t);

    IntersectionClass lineSegmentPlaneIntersection(const Plane &p,
                                                   const LineSegment &line,
                                                   Vector &v);

    RayIntersectionClass rayRayIntersection(const Ray &r1,
                                            const Ray &r2,
                                            Vector &v1,
                                            Vector &v2,
                                            double &mu1,
                                            double &mu2);



    // test whether point d is above, below or on the plane formed by the triangle a,b,c.
    // return: +ve = d is below a,b,c
    //         -ve = d is above a,b,c
    //           0 = d is on a,b,c
    static inline double orient3d(const carve::geom3d::Vector &a,
        const carve::geom3d::Vector &b,
        const carve::geom3d::Vector &c,
        const carve::geom3d::Vector &d) {
      return dotcross((a - d), (b - d), (c - d));
    }



    // Volume of a tetrahedron described by 4 points. Will be
    // positive if the anticlockwise normal of a,b,c is oriented out
    // of the tetrahedron.
    //
    // see: http://mathworld.wolfram.com/Tetrahedron.html
    static inline double tetrahedronVolume(const carve::geom3d::Vector &a,
        const carve::geom3d::Vector &b,
        const carve::geom3d::Vector &c,
        const carve::geom3d::Vector &d) {
      return dotcross((a - d), (b - d), (c - d)) / 6.0;
    }



    // The anticlockwise angle from vector "from" to vector "to", oriented around the vector "orient".
    static inline double antiClockwiseAngle(const Vector &from, const Vector &to, const Vector &orient) {
      double dp = dot(from, to);
      Vector cp = cross(from, to);
      if (cp.isZero()) {
        if (dp < 0) {
          return M_PI;
        } else {
          return 0.0;
        }
      } else {
        if (dot(cp, orient) > 0.0) {
          return acos(dp);
        } else {
          return M_TWOPI - acos(dp);
        }
      }
    }



    static inline double antiClockwiseOrdering(const Vector &from, const Vector &to, const Vector &orient) {
      double dp = dot(from, to);
      Vector cp = cross(from, to);
      if (cp.isZero()) {
        if (dp < 0) {
          return 2.0;
        } else {
          return 0.0;
        }
      } else {
        if (dot(cp, orient) > 0.0) {
          // 1..-1 -> 0..2
          return 1.0 - dp;
        } else {
          // -1..1 -> 2..4
          return dp + 1.0;
        }
      }
    }



  }
}
