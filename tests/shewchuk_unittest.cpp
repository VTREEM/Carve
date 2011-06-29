// Begin License:
// Copyright (C) 2006-2011 Tobias Sargeant (tobias.sargeant@gmail.com).
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

#include <gtest/gtest.h>

#if defined(HAVE_CONFIG_H)
#  include <carve_config.h>
#endif

#include <carve/carve.hpp>
#include "../lib/shewchuk_predicates.cpp"

TEST(ExactTest, ExactTest) {
  double vala,valb,valc,vald;

  vala = 805306457.0;
  valb = 1610612741.0;
  valc = 3221225473.0;
  vald = 4294967291.0;

  double c;
  INEXACT double bvirt;
  INEXACT double abig;
  double avirt, bround, around;
  double ahi, alo, bhi, blo;
  double a0hi, a0lo, a1hi, a1lo;
  double err1, err2, err3;
  INEXACT double _i, _j, _k, _l, _m, _n;
  double _0, _1, _2;

  double aa1, aa0, as1, as0;
  double ab1, ab0, cd1, cd0;
  Two_Product(vala, valb, ab1, ab0);
  Two_Product(valc, vald, cd1, cd0);

  Two_Product(vala, vala, aa1, aa0);
  Square(vala, as1, as0);

  EXPECT_EQ(aa1, as1);
  EXPECT_EQ(aa0, as0);

  double a4[6], a4p[8];

  Two_Two_Product(aa1, aa0, aa1, aa0, a4p[7], a4p[6], a4p[5], a4p[4], a4p[3], a4p[2], a4p[1], a4p[0]);
  Two_Square(aa1, aa0, a4[5], a4[4], a4[3], a4[2], a4[1], a4[0]);

  std::cerr << "prod: " << a4p[7] << "," << a4p[6] << "," << a4p[5] << "," << a4p[4] << "," << a4p[3] << "," << a4p[2] << "," << a4p[1] << "," << a4p[0] << std::endl;
  std::cerr << "square: " << a4[5] << "," << a4[4] << "," << a4[3] << "," << a4[2] << "," << a4[1] << "," << a4[0] << std::endl;

  double x[8], y[8];
  std::fill(x+0,x+8,0.0);
  std::fill(y+0,y+8,0.0);
  Two_Two_Product(ab1, ab0, cd1, cd0, x[7], x[6], x[5], x[4], x[3], x[2], x[1], x[0]);

  std::cerr << "ab: " << ab1 << "," << ab0 << std::endl;
  std::cerr << "cd: " << cd1 << "," << cd0 << std::endl;
  std::cerr << "prod: " << x[7] << "," << x[6] << "," << x[5] << "," << x[4] << "," << x[3] << "," << x[2] << "," << x[1] << "," << x[0] << std::endl;

  int n = shewchuk::compress(4, x, y);

  std::cerr << "comp: " << n << " " << y[7] << "," << y[6] << "," << y[5] << "," << y[4] << "," << y[3] << "," << y[2] << "," << y[1] << "," << y[0] << std::endl;
}


TEST(ExactTest, FastExpansionSumZeroelim) {
  double d[8];
  std::fill(d, d+8, 0.0);

  double a[2];
  double b[2];
  a[0] = -376.0;          a[1] = 18158514394416283648.0;         
  b[0] = -543313364724.0; b[1] = 9748778911853561743662907392.0; 

  int dlen = shewchuk::fast_expansion_sum_zeroelim(2, a, 2, b, d);

  std::cerr << "sum: ";
  for (int i = 0; i < dlen; ++i) {
    std::cerr << " " << d[i];
  }
  std::cerr << std::endl;
}
