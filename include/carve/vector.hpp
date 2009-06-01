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

#include <carve/math_constants.hpp>
#include <carve/geom.hpp>

#include <iostream>
#include <sstream>

#include <math.h>

namespace carve {
  namespace geom3d {

    typedef carve::geom::vector<3> Vector;

    static inline double tetrahedronVolume(const carve::geom3d::Vector &a,
                                           const carve::geom3d::Vector &b,
                                           const carve::geom3d::Vector &c,
                                           const carve::geom3d::Vector &d) {
      // Volume of a tetrahedron described by 4 points. Will be
      // positive if the anticlockwise normal of a,b,c is oriented out
      // of the tetrahedron.
      //
      // see: http://mathworld.wolfram.com/Tetrahedron.html
      return dotcross((a - d), (b - d), (c - d)) / 6.0;
    }



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



    struct hash_vector_ptr {
      size_t operator()(const Vector * const &v) const {
        return (size_t)v;
      }
      size_t operator()(const std::pair<const Vector *, const Vector *> &v) const {
        size_t r = (size_t)v.first;
        size_t s = (size_t)v.second;
        return r ^ ((s >> 16) | (s << 16));
      }
    };



    struct vec_adapt_ident {
      const Vector &operator()(const Vector &v) { return v; }
      Vector &operator()(Vector &v) { return v; }
    };



    struct vec_adapt_ptr {
      const Vector &operator()(const Vector * const &v) { return *v; }
      Vector &operator()(Vector *&v) { return *v; }
    };


  
    struct vec_adapt_pair_first {
      template<typename pair_t> const Vector &operator()(const pair_t &v) { return v.first; }
      template<typename pair_t> Vector &operator()(pair_t &v) { return v.first; }
    };



    struct vec_adapt_pair_second {
      template<typename pair_t> const Vector &operator()(const pair_t &v) { return v.second; }
      template<typename pair_t> Vector &operator()(pair_t &v) { return v.second; }
    };



    template<typename adapt_t>
    struct vec_cmp_lt_x {
      adapt_t adapt;
      vec_cmp_lt_x(adapt_t _adapt) : adapt(_adapt) {}
      template<typename input_t> bool operator()(input_t &a, input_t &b) { return adapt(a).x < adapt(b).x; }
    };
    template<typename adapt_t> vec_cmp_lt_x<adapt_t> vec_lt_x(adapt_t &adapt) { return vec_cmp_lt_x<adapt_t>(adapt); }



    template<typename adapt_t>
    struct vec_cmp_lt_y {
      adapt_t adapt;
      vec_cmp_lt_y(adapt_t _adapt) : adapt(_adapt) {}
      template<typename input_t> bool operator()(input_t &a, input_t &b) { return adapt(a).y < adapt(b).y; }
    };
    template<typename adapt_t> vec_cmp_lt_y<adapt_t> vec_lt_y(adapt_t &adapt) { return vec_cmp_lt_y<adapt_t>(adapt); }



    template<typename adapt_t>
    struct vec_cmp_lt_z {
      adapt_t adapt;
      vec_cmp_lt_z(adapt_t _adapt) : adapt(_adapt) {}
      template<typename input_t> bool operator()(input_t &a, input_t &b) { return adapt(a).z < adapt(b).z; }
    };
    template<typename adapt_t> vec_cmp_lt_z<adapt_t> vec_lt_z(adapt_t &adapt) { return vec_cmp_lt_z<adapt_t>(adapt); }



    template<typename adapt_t>
    struct vec_cmp_gt_x {
      adapt_t adapt;
      vec_cmp_gt_x(adapt_t _adapt) : adapt(_adapt) {}
      template<typename input_t> bool operator()(input_t &a, input_t &b) { return adapt(a).x > adapt(b).x; }
    };
    template<typename adapt_t> vec_cmp_gt_x<adapt_t> vec_gt_x(adapt_t &adapt) { return vec_cmp_gt_x<adapt_t>(adapt); }



    template<typename adapt_t>
    struct vec_cmp_gt_y {
      adapt_t adapt;
      vec_cmp_gt_y(adapt_t _adapt) : adapt(_adapt) {}
      template<typename input_t> bool operator()(input_t &a, input_t &b) { return adapt(a).y > adapt(b).y; }
    };
    template<typename adapt_t> vec_cmp_gt_y<adapt_t> vec_gt_y(adapt_t &adapt) { return vec_cmp_gt_y<adapt_t>(adapt); }



    template<typename adapt_t>
    struct vec_cmp_gt_z {
      adapt_t adapt;
      vec_cmp_gt_z(adapt_t _adapt) : adapt(_adapt) {}
      template<typename input_t> bool operator()(input_t &a, input_t &b) { return adapt(a).z > adapt(b).z; }
    };
    template<typename adapt_t> vec_cmp_gt_z<adapt_t> vec_gt_z(adapt_t &adapt) { return vec_cmp_gt_z<adapt_t>(adapt); }



    template<typename iter_t, typename adapt_t>
    void sortInDirectionOfRay(const Vector &ray_dir, iter_t begin, iter_t end, adapt_t adapt) {
      switch (carve::geom::dominantAxis(ray_dir)) {
      case 0:
        if (ray_dir.x > 0) {
          std::sort(begin, end, vec_lt_x(adapt));
        } else {
          std::sort(begin, end, vec_gt_x(adapt));
        }
        break;
      case 1:
        if (ray_dir.y > 0) {
          std::sort(begin, end, vec_lt_y(adapt));
        } else {
          std::sort(begin, end, vec_gt_y(adapt));
        }
        break;
      case 2:
        if (ray_dir.z > 0) {
          std::sort(begin, end, vec_lt_z(adapt));
        } else {
          std::sort(begin, end, vec_gt_z(adapt));
        }
        break;
      }
    }

  }
}
