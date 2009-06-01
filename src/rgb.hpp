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

#undef RGB

#include <functional>

struct RGB {
  typedef float value_type;
  value_type r, g, b;

  RGB() : r(0), g(0), b(0) { }

  template <typename T>
  RGB(T _r, T _g, T _b) : r((value_type)_r), g((value_type)_g), b((value_type)_b) { }
};

struct RGBA {
  typedef float value_type;
  value_type r, g, b, a;

  RGBA() : r(0),g(0),b(0),a(1) { }
  template <typename T>
  RGBA(T _r, T _g, T _b, T _a = T(1)) : r((value_type)_r), g((value_type)_g), b((value_type)_b), a((value_type)_a) { }

  RGBA(const RGB &rgb) : r(rgb.r), g(rgb.g), b(rgb.b), a(1) { }
};

static inline RGB operator+(const RGB &a, const RGB &b) {
  return RGB(a.r + b.r, a.g + b.g, a.b + b.b);
}
static inline RGB &operator+=(RGB &a, const RGB &b) {
  a.r += b.r; a.g += b.g; a.b += b.b;
  return a;
}

static inline RGB operator*(double s, const RGB &a) {
  return RGB(s * a.r, s * a.g, s * a.b);
}

static inline RGBA operator+(const RGBA &a, const RGBA &b) {
  return RGBA(a.r + b.r, a.g + b.g, a.b + b.b, a.a + b.a);
}
static inline RGBA &operator+=(RGBA &a, const RGBA &b) {
  a.r += b.r; a.g += b.g; a.b += b.b; a.a += b.a;
  return a;
}

static inline RGBA operator*(double s, const RGBA &a) {
  return RGBA(s * a.r, s * a.g, s * a.b, s * a.a);
}

static inline RGB HSV2RGB(float H, float S, float V) {
  H = 6.0f * H;
  if (S < 5.0e-6) {
    RGB(V, V, V);
  } else {
    int i = (int)H;
    float f = H - i;
    float p1 = V * (1.0f - S);
    float p2 = V * (1.0f - S * f);
    float p3 = V * (1.0f - S * (1.0f - f));
    switch (i) {
    case 0: return RGB(V, p3, p1);
    case 1: return RGB(p2,  V, p1);
    case 2: return RGB(p1,  V, p3);
    case 3: return RGB(p1, p2,  V);
    case 4: return RGB(p3, p1,  V);
    case 5: return RGB(V, p1, p2);
    }
  }
  return RGB(0, 0, 0);
}

struct colour_clamp_t {
  RGB operator()(const RGB &c) const {
    return RGB(std::min(std::max(c.r, 0.0f), 1.0f),
               std::min(std::max(c.g, 0.0f), 1.0f),
               std::min(std::max(c.b, 0.0f), 1.0f));
  }
  RGBA operator()(const RGBA &c) const {
    return RGBA(std::min(std::max(c.r, 0.0f), 1.0f),
                std::min(std::max(c.g, 0.0f), 1.0f),
                std::min(std::max(c.b, 0.0f), 1.0f),
                std::min(std::max(c.a, 0.0f), 1.0f));
  }
};
