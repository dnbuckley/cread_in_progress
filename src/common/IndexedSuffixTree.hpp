/*
 * Copyright (C) 2006 Cold Spring Harbor Laboratory
 * Authors: Andrew D. Smith and Michael Q. Zhang
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

#ifndef SUFFIXTREE_H
#define SUFFIXTREE_H

#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <utility>
#include <iostream>
#include <iterator>
#include <limits>
#include <algorithm>
#include <queue>

#include "cread.hpp"
#include "ScoringMatrix.hpp"

typedef std::pair<size_t, float> pos_score;

class IndexedSuffixNode;

class IndexedSuffixTree {
  
  std::string sequence;
  IndexedSuffixNode* root;
  
  static bool valid(char b) {
    return b == 'A' || b == 'C' || b == 'G' || b == 'T';
  }
  
  IndexedSuffixTree();
  IndexedSuffixTree(const IndexedSuffixTree&);
  IndexedSuffixTree& operator=(const IndexedSuffixTree&);
  static float** PrepareScoringMatrix(const ScoringMatrix&, float &);
  static void destroy_st_scoring_matrix(float **sm, const size_t w);
  
public:
  IndexedSuffixTree(const std::string&);
  ~IndexedSuffixTree();
  
  std::string tostring() const;
  
  void scores_greater_indices(const ScoringMatrix&, const float cutoff,
			      std::vector<pos_score> &) const;
  void top_scores_indices(const ScoringMatrix&, const size_t,
			  std::vector<pos_score>&) const;
};

#endif
