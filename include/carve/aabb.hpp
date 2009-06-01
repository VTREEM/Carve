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

#include <carve/vector.hpp>
#include <carve/geom3d.hpp>

#include <carve/geom.hpp>

#include <vector>
#include <iostream>

namespace carve {
  namespace geom {
    // n-dimensional AABB
    template<unsigned ndim>
    struct aabb {
      typedef vector<ndim> vector_t;

      vector_t pos;     // the centre of the AABB
      vector_t extent;  // the extent of the AABB - the vector from the centre to the maximal vertex.

      void empty() {
        pos.setZero();
        extent.setZero();
      }

      bool isEmpty() const {
        return extent.exactlyZero();
      }

      void fit(const vector_t &v1) {
        pos = v1;
        extent.setZero();
      }

      void fit(const vector_t &v1, const vector_t &v2) {
        vector_t min, max;
        assign_op(min, v1, v2, carve::util::min_functor());
        assign_op(max, v1, v2, carve::util::max_functor());

        pos = (min + max) / 2.0;
        extent = max - pos;
      }

      void fit(const vector_t &v1, const vector_t &v2, const vector_t &v3) {
        vector_t min, max;
        min = max = v1;

        assign_op(min, min, v2, carve::util::min_functor());
        assign_op(max, max, v2, carve::util::max_functor());
        assign_op(min, min, v3, carve::util::min_functor());
        assign_op(max, max, v3, carve::util::max_functor());

        pos = (min + max) / 2.0;
        extent = max - pos;
      }

      template<typename iter_t, typename adapt_t>
      void fit(iter_t begin, iter_t end, adapt_t adapt) {
        vector_t min, max;
        bounds(begin, end, adapt, min, max);

        pos = (min + max) / 2.0;
        extent = max - pos;
      }

      template<typename iter_t>
      void fit(iter_t begin, iter_t end) {
        vector_t min, max;
        bounds(begin, end, min, max);

        pos = (min + max) / 2.0;
        extent = max - pos;
      }

      template<typename iter_t>
      void fitAABB(iter_t begin, iter_t end) {
        if (begin == end) {
          empty();
        } else {
          vector_t min, max;
          aabb<ndim> a = *begin++;
          min = a.min();
          max = a.max();
          while (begin != end) {
            aabb<ndim> a = *begin; ++begin;
            assign_op(min, min, a.min(), carve::util::min_functor());
            assign_op(max, max, a.max(), carve::util::max_functor());
          }

          pos = (min + max) / 2.0;
          extent = max - pos;
        }
      }

      void expand(double pad) {
        extent += pad;
      }

      void fitAABB(const aabb<ndim> &a, const aabb<ndim> &b) {
        vector_t min, max;
        assign_op(min, a.min(), b.min(), carve::util::min_functor());
        assign_op(max, a.max(), b.max(), carve::util::max_functor());

        pos = (min + max) / 2.0;
        extent = max - pos;
      }

      void unionAABB(const aabb<ndim> &a) {
        vector_t vmin, vmax;
        assign_op(vmin, min(), a.min(), carve::util::min_functor());
        assign_op(vmax, max(), a.max(), carve::util::max_functor());

        pos = (vmin + vmax) / 2.0;
        extent = vmax - pos;
      }

      aabb(const vector_t &_pos = vector_t::ZERO(),
           const vector_t &_extent = vector_t::ZERO()) : pos(_pos), extent(_extent) {
      }

      template<typename iter_t, typename adapt_t>
      aabb(iter_t begin, iter_t end, adapt_t adapt) {
        fit(begin, end, adapt);
      }

      template<typename iter_t, typename adapt_t>
      aabb(iter_t begin, iter_t end) {
        fit(begin, end);
      }

      aabb(const aabb<ndim> &a, const aabb<ndim> &b) {
        fit(a, b);
      }

      bool intersects(const aabb<ndim> &other) const {
        for (unsigned i = 0; i < ndim; ++i) {
          if (fabs(other.pos.v[i] - pos.v[i]) > (extent.v[i] + other.extent.v[i])) return false;
        }
        return true;
      }
 
      bool completelyContains(const aabb<ndim> &other) const {
        for (unsigned i = 0; i < ndim; ++i) {
          if (fabs(other.pos.v[i] - pos.v[i] + other.extent.v[i]) > extent.v[i]) return false;
        }
        return true;
      }

     bool containsPoint(const vector_t &v) const {
        for (unsigned i = 0; i < ndim; ++i) {
          if (fabs(v.v[i] - pos.v[i]) > extent.v[i]) return false;
        }
        return true;
      }

      bool intersectsSphere(const vector_t &v, double d) const {
        double r = 0.0;
        for (unsigned i = 0; i < ndim; ++i) {
          double t = fabs(v[i] - pos[i]) - extent[i]; if (t > 0.0) r += t*t;
        }
        return r <= d*d;
      }

      bool intersectsPlane(const plane<ndim> &plane) const {
        double d1 = fabs(distance(plane, pos));
        double d2 = dot(abs(plane.N), extent);
        return d1 <= d2;
      }

      bool intersectsRay(const ray<ndim> &ray) const;
      bool intersectsLineSegment(const vector_t &v1, const vector_t &v2) const;

      vector_t min() const { return pos - extent; }
      vector_t max() const { return pos + extent; }

      int compareAxis(const axis_pos &ap) const {
        double p = ap.pos - pos[ap.axis];
        if (p > extent[ap.axis]) return -1;
        if (p < -extent[ap.axis]) return +1;
        return 0;
      }

      void constrainMax(const axis_pos &ap) {
        if (pos[ap.axis] + extent[ap.axis] > ap.pos) {
          double min = std::min(ap.pos, pos[ap.axis] - extent[ap.axis]);
          pos[ap.axis] = (min + ap.pos) / 2.0;
          extent[ap.axis] = ap.pos - pos[ap.axis];
        }
      }

      void constrainMin(const axis_pos &ap) {
        if (pos[ap.axis] - extent[ap.axis] < ap.pos) {
          double max = std::max(ap.pos, pos[ap.axis] + extent[ap.axis]);
          pos[ap.axis] = (ap.pos + max) / 2.0;
          extent[ap.axis] = pos[ap.axis] - ap.pos;
        }
      }

    };

    template<unsigned ndim>
    bool operator==(const aabb<ndim> &a, const aabb<ndim> &b) {
      return a.pos == b.pos && a.extent == b.extent;
    }

    template<unsigned ndim>
    bool operator!=(const aabb<ndim> &a, const aabb<ndim> &b) {
      return a.pos != b.pos || a.extent != b.extent;
    }

    template<unsigned ndim>
    std::ostream &operator<<(std::ostream &o, const aabb<ndim> &a) {
      o << (a.pos - a.extent) << "--" << (a.pos + a.extent);
      return o;
    }

    template<>
    inline bool aabb<3>::intersectsRay(const ray<3> &ray) const {
      vector<3> t = pos - ray.v;
      vector<3> v;
      double r;

      //l.cross(x-axis)?
      r = extent.y * fabs(ray.D.z) + extent.z * fabs(ray.D.y);
      if (fabs(t.y * ray.D.z - t.z * ray.D.y) > r) return false;

      //ray.D.cross(y-axis)?
      r = extent.x * fabs(ray.D.z) + extent.z * fabs(ray.D.x);
      if (fabs(t.z * ray.D.x - t.x * ray.D.z) > r) return false;

      //ray.D.cross(z-axis)?
      r = extent.x*fabs(ray.D.y) + extent.y*fabs(ray.D.x);
      if (fabs(t.x * ray.D.y - t.y * ray.D.x) > r) return false;

      return true;
    }

    template<>
    inline bool aabb<3>::intersectsLineSegment(const vector<3> &v1, const vector<3> &v2) const {
      vector<3> half_length = 0.5 * (v2 - v1);
      vector<3> t = pos - half_length - v1;
      vector<3> v;
      double r;

      //do any of the principal axes form a separating axis?
      if(fabs(t.x) > extent.x + fabs(half_length.x)) return false;
      if(fabs(t.y) > extent.y + fabs(half_length.y)) return false;
      if(fabs(t.z) > extent.z + fabs(half_length.z)) return false;

      // NOTE: Since the separating axis is perpendicular to the line in
      // these last four cases, the line does not contribute to the
      // projection.

      //line.cross(x-axis)?
      r = extent.y * fabs(half_length.z) + extent.z * fabs(half_length.y);
      if (fabs(t.y * half_length.z - t.z * half_length.y) > r) return false;

      //half_length.cross(y-axis)?
      r = extent.x * fabs(half_length.z) + extent.z * fabs(half_length.x);
      if (fabs(t.z * half_length.x - t.x * half_length.z) > r) return false;

      //half_length.cross(z-axis)?
      r = extent.x*fabs(half_length.y) + extent.y*fabs(half_length.x);
      if (fabs(t.x * half_length.y - t.y * half_length.x) > r) return false;

      return true;
    }
  }
}

namespace carve {
  namespace geom3d {
    typedef carve::geom::aabb<3> AABB;
  }
}
