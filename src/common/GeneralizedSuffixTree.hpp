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

#ifndef GENERALIZED_SUFFIXTREE_H
#define GENERALIZED_SUFFIXTREE_H

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

#include "ScoringMatrix.hpp"
#include "cread.hpp"

typedef std::pair<std::pair<size_t, size_t>, float> seq_pos_score;

class GeneralizedSuffixNode;

class GeneralizedSuffixTree {
  
  std::string sequence;
  size_t depth;
  GeneralizedSuffixNode* root;
  std::vector<size_t> offset;
  
  std::pair<size_t, size_t> index2seq_offset(size_t) const;
  bool ValidOccurrenceLocation(size_t, size_t) const;
  static bool valid(char b) {
    return b == 'A' || b == 'C' || b == 'G' || b == 'T';
  }
  
  static float** PrepareScoringMatrix(const ScoringMatrix&, float &);
  
  GeneralizedSuffixTree();
  GeneralizedSuffixTree(const GeneralizedSuffixTree&);
  GeneralizedSuffixTree& operator=(const GeneralizedSuffixTree&);

public:
  GeneralizedSuffixTree(const std::vector<std::string>&, int = -1);
  GeneralizedSuffixTree(const std::string&, int = -1);
  ~GeneralizedSuffixTree();
  
  std::string tostring() const;
  // Get scores, but don't return the index
  
  /* The final argument, which defaults to '-1', allows the user to
     specify a particular sequence. If no sequence is specified, then
     all sequences are searched */
  // Get the top score, option to specify sequence
  float top_score(const ScoringMatrix&, int = -1) const;
  float top_score(float**, size_t, int = -1) const;
  // get top score within a window
  float window_top_score(float **, size_t, int = -1) const;
  
  /* These next two are the same as above, except they can take a range */
  float top_score(const ScoringMatrix&, size_t, size_t, int = -1) const;

  // Get the top k scores, option to specify sequence
  void top_scores(std::vector<float>&, size_t, const ScoringMatrix&, int = -1) const;
  // Get scores greater than a threshold, option to specify sequence
  void scores_greater(std::vector<float>&, float, 
		      const ScoringMatrix&, int = -1) const;
  
  static bool no_overlaps(std::vector<std::pair<size_t, size_t> > &,
			  size_t, size_t);
  void window_scores_greater(std::vector<float>&, float,
			     const ScoringMatrix&, int, int) const;
  void window_scores_greater(std::vector<float>&, float, float**, size_t,
			     size_t, size_t) const;

  // Get scores and indices corresponding to those scores
  
  // Get the occurrence and its score, option to specify sequence
  seq_pos_score top_score_index(const ScoringMatrix&, int = -1) const;

  seq_pos_score top_score_index(const ScoringMatrix&, 
				size_t, size_t, int = -1) const;
  
  // Get the top k scores with indices, option to specify sequence
  void top_scores_indices(std::vector<seq_pos_score>&, size_t, 
			  const ScoringMatrix&, int = -1) const;

  // Get scores and indices with scores greater than threshold, option
  // to specify sequence
  void scores_greater_indices(std::vector<seq_pos_score>&, float, 
			      const ScoringMatrix&, int = -1) const;

  // get scores (and indices) above threshold for sequence within a window
  void window_scores_greater_indices(std::vector<seq_pos_score>&, float, 
				     const ScoringMatrix&, size_t, size_t, size_t) const;
  void window_scores_greater_indices(std::vector<seq_pos_score>&, float, 
				     const ScoringMatrix&, size_t, size_t, size_t, 
				     std::vector<std::pair<size_t, size_t> >&) const;
  void window_scores_greater_indices(std::vector<seq_pos_score>&, float, 
				     const ScoringMatrix&, const ScoringMatrix&, 
				     size_t, size_t, size_t) const;
  
  /************************************************/
  // Get locations of substrings of width in a particular window,
  // and number of copies greater than a particular threshold
  void GetRepeatLocations(std::vector<std::pair<size_t, size_t> >& V, 
			  size_t threshold_number, size_t min_width) const;
  void get_kmer_counts(std::vector<std::pair<std::string, size_t> > &,
		       size_t) const;
  size_t get_count(std::string w) const;
  void get_locations(std::string w, 
		     std::vector<std::pair<size_t, size_t> > &L) const;
  
};

#endif
