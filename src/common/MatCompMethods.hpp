/*
 * Copyright (C) 2006 Cold Spring Harbor Laboratory
 * Authors: Andrew D. Smith, Pavel Sumazin and Michael Q. Zhang
 *
 * This file is part of CREAD.
 *
 * CREAD is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * CREAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CREAD; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef MAT_COMP_METHODS_HPP
#define MAT_COMP_METHODS_HPP

/*!
  \file MatCompMethods.hpp

  \brief Organizes functions for comparing position-weight matrices.

*/

#include "Matrix.hpp"
#include "cread.hpp"

#ifdef HAVE_GSL
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#endif

class MatCompMethods {
public:
  static float divergence(const Matrix&, const Matrix&);
  static float chisquared(const Matrix&, const Matrix&);
  static float fishertest(const Matrix&, const Matrix&);
  
  static float sliding_divergence(const Matrix&, const Matrix&,
				  size_t max_overhang);
  static float sliding_divergence_offset(const Matrix&, const Matrix&,
					 int& offset, size_t max_overhang);
  static float sliding_chisquared(const Matrix&, const Matrix&, 
				  size_t max_overhang);
  static int sliding_chisquared_offset(const Matrix&, const Matrix&, 
				       size_t max_overhang);
  static float sliding_fishertest(const Matrix&, const Matrix&,
				  size_t max_overhang);
  static int sliding_fishertest_offset(const Matrix&, const Matrix&,
                                       size_t max_overhang);
private:

  struct column_struct {
    int top[alphabet_size];
    int bottom[alphabet_size];
    column_struct(int t[alphabet_size], int b[alphabet_size]) {
      std::copy(t, t + alphabet_size, top);
      std::copy(b, b + alphabet_size, bottom);
    }
  };
  MatCompMethods() {} // no objects
  static float sliding(const Matrix&, const Matrix&, 
		       float (*)(const float *a1[], const float *a2[],
				 const float *b1[], const float *b2[],
				 float (*)(const float [], const float [])),
		       float (*)(const float[], const float[]),
		       int &offset, size_t = 0);
  static float column_set_sum(const float *a1[], const float *a2[],
			      const float *b1[], const float *b2[],
			      float (*)(const float [], const float []));
  static float column_set_product(const float *a1[], const float *a2[],
				  const float *b1[], const float *b2[],
				  float (*)(const float [], const float []));
  static float column_divergence(const float a[], const float b[]);
  static float column_chisquared(const float a[], const float b[]);
  static float column_fishertest(const float a[], const float b[]);
  static float epsilon() {return 0.005;}
  static void build_all_tables(int column_total[], int row_total[], 
			       std::vector<column_struct>& column_matrices);
};
#endif
