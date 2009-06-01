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
#include <carve/geom2d.hpp>
#include <carve/poly.hpp>
#include <carve/csg.hpp>

namespace carve {
  namespace interpolate {

    template<typename T, typename adapt_t>
    std::vector<double> polyInterpolate(const std::vector<T> &points, adapt_t adapt, const carve::geom2d::P2 &v) {
      // see hormann et al. 2006
      const size_t SZ = points.size();
      std::vector<carve::geom2d::P2> s;
      std::vector<double> r;
      std::vector<double> A;
      std::vector<double> D;

      std::vector<double> result;

      s.resize(SZ);
      r.resize(SZ);
      A.resize(SZ);
      D.resize(SZ);

      result.resize(SZ, 0.0);

      for (size_t i = 0; i < SZ; ++i) s[i] = adapt(points[i]) - v;

      for (size_t i = 0; i < SZ; ++i) {
        size_t i2 = (i + 1) % SZ;

        r[i] = sqrt(dot(s[i], s[i]));
        A[i] = cross(s[i], s[i2]) / 2.0;
        D[i] = dot(s[i], s[i2]);
        if (fabs(r[i]) < 1e-16) {
          result[i] = 1.0;
          return result;
        } else if (fabs(A[i]) < 1e-16 && D[i] < 0.0) {
          double r2 = sqrt(dot(s[i2], s[i2]));
          result[i2] = r[i] / (r[i] + r2);
          result[i] = r2 / (r[i] + r2);
          return result;
        }
      }

      double w_sum = 0.0;

      for (size_t i = 0; i < SZ; ++i) {
        size_t i_m = (i + SZ - 1) % SZ;
        size_t i_p = (i + 1) % SZ;

        double w = 0.0;
        if (fabs(A[i_m]) > 1e-16)
          w += (r[i_m] - D[i_m] / r[i]) / A[i_m];
        if (fabs(A[i]) > 1e-16)
          w += (r[i_p] - D[i] / r[i]) / A[i];

        result[i] = w;
        w_sum += w;
      }
  
      for (size_t i = 0; i < SZ; ++i) {
        result[i] /= w_sum;
      }

      carve::geom2d::P2 test;
      for (size_t i = 0; i < SZ; ++i) {
        test = test + result[i] * s[i];
      }
      // std::cerr << "test = " << test << std::endl;

      return result;
    }

    std::vector<double> polyInterpolate(const std::vector<carve::geom2d::P2> &points, const carve::geom2d::P2 &v) {
      return polyInterpolate(points, carve::geom2d::p2_adapt_ident(), v);
    }

    template<typename T,
             typename adapt_t,
             typename val_t,
             typename mod_t>
    val_t interp(const std::vector<T> &poly,
                 adapt_t adapt,
                 const std::vector<val_t> &vals,
                 double x,
                 double y,
                 mod_t mod = mod_t()) {
      std::vector<double> weight = polyInterpolate(poly, adapt, carve::geom::VECTOR(x, y));
      val_t v;
      for (size_t z = 0; z < poly.size(); z++) {
        v += weight[z] * vals[z];
      }

      return mod(v);
    }

    template<typename T,
             typename adapt_t,
             typename val_t>
    val_t interp(const std::vector<T> &poly,
                 adapt_t adapt,
                 const std::vector<val_t> &vals,
                 double x,
                 double y) {
      return interp(poly, adapt, vals, x, y, identity_t<val_t>());
    }

    template<typename val_t,
             typename mod_t>
    val_t interp(const std::vector<carve::geom2d::P2> &poly,
                 const std::vector<val_t> &vals,
                 double x,
                 double y,
                 mod_t mod = mod_t()) {
      return interp(poly, carve::geom2d::p2_adapt_ident(), vals, x, y, mod);
    }

    template<typename val_t>
    val_t interp(const std::vector<carve::geom2d::P2> &poly,
                 const std::vector<val_t> &vals,
                 double x,
                 double y) {
      return interp(poly, carve::geom2d::p2_adapt_ident(), vals, x, y, identity_t<val_t>());
    }



    class Interpolator {
    public:
      virtual void interpolate(const carve::poly::Face<3> *new_face,
                               const carve::poly::Face<3> *orig_face,
                               bool flipped) =0;

      Interpolator() {
      }

      virtual ~Interpolator() {
      }

      class Hook : public carve::csg::CSG::Hook {
        Interpolator *interpolator;
      public:
        virtual void resultFace(const carve::poly::Face<3> *new_face,
                                const carve::poly::Face<3> *orig_face,
                                bool flipped) {
          interpolator->interpolate(new_face, orig_face, flipped);
        }

        Hook(Interpolator *_interpolator) : interpolator(_interpolator) {
        }

        virtual ~Hook() {
        }
      };

      void installHooks(carve::csg::CSG &csg) {
        csg.hooks.registerHook(new Hook(this), carve::csg::CSG::Hooks::RESULT_FACE_BIT);
      }
    };

    template<typename attr_t>
    class FaceVertexAttr : public Interpolator {

    protected:
      struct fv_hash {
        size_t operator()(const std::pair<const carve::poly::Face<3> *, unsigned> &v) const {
          return size_t(v.first) ^ size_t(v.second);
        }
      };

      typedef std::unordered_map<const carve::poly::Vertex<3> *, attr_t, carve::poly::hash_vertex_ptr> attrvmap_t;
      typedef std::unordered_map<std::pair<const carve::poly::Face<3> *, unsigned>, attr_t, fv_hash> attrmap_t;

      attrmap_t attrs;

    public:
      bool hasAttribute(const carve::poly::Face<3> *f, unsigned v) {
        return attrs.find(std::make_pair(f, v)) != attrs.end();
      }

      attr_t getAttribute(const carve::poly::Face<3> *f, unsigned v, const attr_t &def = attr_t()) {
        typename attrmap_t::const_iterator fv = attrs.find(std::make_pair(f, v));
        if (fv != attrs.end()) {
          return (*fv).second;
        }
        return def;
      }

      void setAttribute(const carve::poly::Face<3> *f, unsigned v, const attr_t &attr) {
        attrs[std::make_pair(f, v)] = attr;
      }

      virtual void interpolate(const carve::poly::Face<3> *new_face,
                               const carve::poly::Face<3> *orig_face,
                               bool flipped) {
        std::vector<attr_t> vertex_attrs;
        attrvmap_t base_attrs;
        vertex_attrs.reserve(orig_face->vertices.size());

        for (size_t i = 0; i < orig_face->vertices.size(); ++i) {
          typename attrmap_t::const_iterator a = attrs.find(std::make_pair(orig_face, i));
          if (a == attrs.end()) return;
          vertex_attrs.push_back((*a).second);
          base_attrs[orig_face->vertices[i]] = vertex_attrs.back();
        }

        for (size_t i = 0; i < new_face->vertices.size(); ++i) {
          const carve::poly::Vertex<3> *vertex = new_face->vertices[i];
          typename attrvmap_t::const_iterator b = base_attrs.find(vertex);
          if (b != base_attrs.end()) {
            attrs[std::make_pair(new_face, i)] = (*b).second;
          } else {
            carve::geom2d::P2 p = carve::poly::face::project(orig_face, new_face->vertices[i]->v);
            attr_t attr = interp(orig_face->vertices,
                                 carve::poly::p2_adapt_project<3>(orig_face->project),
                                 vertex_attrs,
                                 p.x,
                                 p.y);
            attrs[std::make_pair(new_face, i)] = attr;
          }
        }
      }

      FaceVertexAttr() : Interpolator() {
      }

      virtual ~FaceVertexAttr() {
      }

    };


    template<typename attr_t>
    class FaceAttr : public Interpolator {

    protected:
      struct f_hash {
        size_t operator()(const carve::poly::Face<3> * const &f) const {
          return size_t(f);
        }
      };

      typedef std::unordered_map<const carve::poly::Face<3> *, attr_t, f_hash> attrmap_t;

      attrmap_t attrs;

    public:
      bool hasAttribute(const carve::poly::Face<3> *f) {
        return attrs.find(f) != attrs.end();
      }

      attr_t getAttribute(const carve::poly::Face<3> *f, const attr_t &def = attr_t()) {
        typename attrmap_t::const_iterator i = attrs.find(f);
        if (i != attrs.end()) {
          return (*i).second;
        }
        return def;
      }

      void setAttribute(const carve::poly::Face<3> *f, const attr_t &attr) {
        attrs[f] = attr;
      }

      virtual void interpolate(const carve::poly::Face<3> *new_face,
                               const carve::poly::Face<3> *orig_face,
                               bool flipped) {
        typename attrmap_t::const_iterator i = attrs.find(orig_face);
        if (i != attrs.end()) {
          attrs[new_face] = (*i).second;
        }
      }

      FaceAttr() : Interpolator() {
      }

      virtual ~FaceAttr() {
      }

    };

  }
}
