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

#include "MatCompMethods.hpp"

using std::vector;

float
MatCompMethods::column_divergence(const float a[], const float b[]) {
  float divergence = 0;
  for (size_t i = 0; i < alphabet_size; i++) {
    float x = std::max(a[i], epsilon());
    float y = std::max(b[i], epsilon());
    divergence += (x - y) * log(x/y);
  }
  return divergence;
}

#ifdef HAVE_GSL
float
MatCompMethods::column_chisquared(const float a[], const float b[]) {
  // get row totals & total
  float row1_frac = std::accumulate(a, a + alphabet_size, 0.0);
  float row2_frac = std::accumulate(b, b + alphabet_size, 0.0);
  float total = row1_frac + row2_frac;
  row1_frac /= total;
  row2_frac /= total;
  
  // now do calculation
  float chi = 0;
  for (size_t j = 0; j < alphabet_size; ++j) {
    float expect_a = (a[j] + b[j])*row1_frac;
    if (expect_a)
      chi += (a[j] - expect_a)*(a[j] - expect_a)/expect_a;
    float expect_b = (a[j] + b[j])*row2_frac;
    if (expect_b)
      chi += (b[j] - expect_b)*(b[j] - expect_b)/expect_b;
  }
  return gsl_cdf_chisq_Q(chi, alphabet_size - 1);
}
#endif

void
MatCompMethods::build_all_tables(int column_total[], int row_total[], 
			     vector<column_struct>& column_matrices){
  int x[alphabet_size];//elements for row 0
  int y[alphabet_size];//elements for row 1
  int remain1[2],remain2[2];
  int lim1 = std::min(column_total[0], row_total[0]);
  for (x[0] = std::max(0, column_total[0] - row_total[1]); x[0] <= lim1; ++x[0]){ 
    y[0] = column_total[0] - x[0];
    remain1[0] = row_total[0] - x[0];
    remain1[1] = row_total[1] - y[0];
    int lim2 = std::min(column_total[1], remain1[0]);
    for (x[1] = std::max(0, column_total[1] - remain1[1]); x[1] <= lim2; ++x[1]) {
      y[1] = column_total[1] - x[1];
      remain2[0] = remain1[0] - x[1];
      remain2[1] = remain1[1] - y[1];
      int lim3 = std::min(column_total[2], remain2[0]);
      for (x[2] = std::max(0, column_total[2] - remain2[1]); x[2] <= lim3; ++x[2]) {
        y[2] = column_total[2] - x[2];
        x[3] = remain2[0] - x[2];
        y[3] = remain2[1] - y[2];
        column_matrices.push_back(column_struct(x,y));
      }
    }
  }
}


float
MatCompMethods::column_fishertest(const float a[], const float b[]){
  long double total_prob = 1;
  int column_totals[alphabet_size] = {0,0,0,0};
  int row_totals[2] = {0,0};
  int total = 0;

  // get row totals & total
  for (size_t j = 0; j < alphabet_size; ++j){
    row_totals[0] += static_cast<int>(a[j]);
    row_totals[1] += static_cast<int>(b[j]);
    total += (static_cast<int>(a[j]) + static_cast<int>(b[j]));
  }

  // get column totals
  for (size_t j = 0; j < alphabet_size; ++j)
    column_totals[j] += (static_cast<int>(a[j]) + static_cast<int>(b[j]));

  // make worse tables
  vector<column_struct> column_matrices;
  build_all_tables(column_totals, row_totals, column_matrices);

  // do calculation
  long double const_prob = lgamma(row_totals[0]+1) + lgamma(row_totals[1]+1) - lgamma(total + 1);
  for (size_t i = 0; i < alphabet_size; ++i)
    const_prob += lgamma(column_totals[i]+1);
  long double prob = 0;
  for (size_t i = 0; i < alphabet_size; ++i)
    prob += lgamma(a[i]+1) + lgamma(b[i]+1);
  long double column_prob = 0;
  for (vector<column_struct>::iterator it = column_matrices.begin(); it != column_matrices.end(); ++it){
    long double this_prob = 0;
    for (size_t i = 0; i < alphabet_size; ++i)
      this_prob += lgamma(it->top[i]+1) + lgamma(it->bottom[i]+1);
    if (this_prob >= prob)
      column_prob += exp(const_prob - this_prob);
  }
  total_prob *= column_prob;
  return total_prob;
}

float
MatCompMethods::column_set_sum(const float *a1[], const float *a2[], 
			       const float *b1[], const float *b2[],
			       float (*cmpr)(const float [], const float [])) {
  size_t n_col = a2 - a1;
  float measure = 0;
  while (a1 < a2 && b1 < b2)
    measure += (*cmpr)(*a1++, *b1++);
  return measure/n_col;
}

float
MatCompMethods::column_set_product(const float *a1[], const float *a2[], 
			       const float *b1[], const float *b2[],
			       float (*cmpr)(const float [], const float [])) {
  size_t n_col = a2 - a1;
  float measure = 1;
  while (a1 < a2 && b1 < b2)
    measure *= (*cmpr)(*a1++, *b1++);
  return (pow((double) measure, 1.0/(float)n_col));
}

float
MatCompMethods::sliding(const Matrix& matrix1, const Matrix& matrix2,
			float (*cmpr_set)(const float *a1[], const float *a2[], 
					  const float *b1[], const float *b2[],
					  float (*cmpr)(const float [],
							const float [])),
			float (*cmpr)(const float[], const float[]),
			int &offset, size_t max_overhang) {
  const Matrix *longer = (matrix1.get_width() > matrix2.get_width()) ?
    &matrix1 : &matrix2;
  const Matrix *shorter = (matrix1.get_width() > matrix2.get_width()) ?
    &matrix2 : &matrix1;
  // float minimum = 0;//std::numeric_limits<float>::max();
  offset = 0;
  
  float best_score;
  if (cmpr == column_divergence)
    best_score = std::numeric_limits<float>::max();
  else
    best_score = 0;

  int sliding_limit = static_cast<int>(longer->get_width() + max_overhang -
				       shorter->get_width() + 1);
    
  for (int i = -max_overhang; i < sliding_limit; ++i) {
    int short_start = std::max(0, -i);
    int long_start = std::max(0, i);
    int short_end = std::min(shorter->get_width(),
			     longer->get_width() - long_start);
    int long_end = std::min(i + shorter->get_width(),
			    longer->get_width());
    float temp_dist = (*cmpr_set)(longer->at(long_start),
				  longer->at(long_end),
				  shorter->at(short_start),
				  shorter->at(short_end),
				  cmpr);

    if (cmpr == column_divergence) {
      if (temp_dist < best_score) {
        best_score = temp_dist;
        offset = i;
      }
    }
    else {
      if (temp_dist > best_score){
        best_score = temp_dist;
        offset = i;
      }
    }
  }
  return best_score;
}

float
MatCompMethods::divergence(const Matrix &matrix1, const Matrix &matrix2) {
  int dummy; // because we don't care what the offset is
  return sliding(matrix1, matrix2, &column_set_sum,
		 &column_divergence, dummy, 0);
}

#ifdef HAVE_GSL
float
MatCompMethods::chisquared(const Matrix &matrix1, const Matrix &matrix2) {
  int dummy;
  return sliding(matrix1, matrix2, &column_set_product,
		 &column_chisquared, dummy, 0);
}
#endif

float
MatCompMethods::fishertest(const Matrix &matrix1, const Matrix &matrix2) {
  int dummy;
  return sliding(matrix1, matrix2, &column_set_product,
                 &column_fishertest, dummy, 0);
}
float
MatCompMethods::sliding_divergence(const Matrix &matrix1, const Matrix &matrix2,
				   size_t max_overhang) {
  int dummy;
  return sliding(matrix1, matrix2, &column_set_sum, &column_divergence,
		 dummy, max_overhang);
}

float
MatCompMethods::sliding_divergence_offset(const Matrix &matrix1, const Matrix &matrix2,
					  int &offset, size_t max_overhang) {
  return sliding(matrix1, matrix2, &column_set_sum, &column_divergence,
		 offset, max_overhang);
}

#ifdef HAVE_GSL
float
MatCompMethods::sliding_chisquared(const Matrix &matrix1, const Matrix &matrix2,
			       size_t max_overhang) {
  int dummy;
  return sliding(matrix1, matrix2, &column_set_product, &column_chisquared,
		 dummy, max_overhang);
}
#endif

#ifdef HAVE_GSL
int
MatCompMethods::sliding_chisquared_offset(const Matrix &matrix1, const Matrix &matrix2,
				      size_t max_overhang) {
  int offset;
  sliding(matrix1, matrix2, &column_set_product, &column_chisquared,
	  offset, max_overhang);
  return offset;
}
#endif

float
MatCompMethods::sliding_fishertest(const Matrix &matrix1, const Matrix &matrix2,
                               size_t max_overhang) {
  int dummy;
  return sliding(matrix1, matrix2, &column_set_product, &column_fishertest,
                 dummy, max_overhang);
}
int
MatCompMethods::sliding_fishertest_offset(const Matrix &matrix1, const Matrix &matrix2,
                                      size_t max_overhang) {
  int offset;
  sliding(matrix1, matrix2, &column_set_product, &column_fishertest,
          offset, max_overhang);
  return offset;
}
