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

#ifndef SUFFIXTREE_H
#define SUFFIXTREE_H

/*!
  \file SuffixTree.hpp
  \brief This header file contains the class definition for SuffixTree.
*/

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

class SuffixNode;

typedef std::pair<float, SuffixNode*> valnode;
typedef std::priority_queue<valnode, 
			    std::vector<valnode>, 
			    std::greater<valnode> > hit_queue;

/*!
  \brief A suffix tree class.
  \class SuffixTree

  This class represents an edge compressed suffix tree over the DNA
  alphabet (ACGTN), and provides functionality to query the tree for
  sites matching scoring matrices.
  
  For many of the querying functions, there is an option to specify
  the sequence in which to search. This is always the final argument,
  and is a number (specifying the index of the sequence in the
  original set, i.e. vector, passed to construct the suffix tree). If
  no sequence is specified, then all sequences are searched.
  
*/
class SuffixTree {
public:

  /*!
    \brief Constructor to make a tree from a single sequence.
    
    This constructor has the option to specify the depth of the tree,
    which defaults to the maximum depth (length of the sequence).
    
    \param the_sequence sequence from which to build the tree
    \param maximum_tree_depth maximum depth of the tree (in terms of
    characters, not nodes).
  */
  SuffixTree(const std::string& the_sequence, 
	     int maximum_tree_depth = FULL_DEPTH);
  /*!
    \brief Constructor to make a tree from a single sequence.
    
    This constructor has the option to specify the depth of the tree,
    which defaults to the maximum depth (length of the sequence).

    \param the_sequence sequence from which to build the tree
    \param maximum_tree_depth maximum depth of the tree (in terms of
    characters, not nodes).
  */
  SuffixTree(const std::vector<std::string>& the_sequences, 
	     int maximum_tree_depth = FULL_DEPTH);

  /*!
    \brief Destructor recursively deletes subtrees.
  */
  ~SuffixTree();
  
  /*!
    \brief Get a string representation.

    This function is mainly used for debugging. For a sequence of
    length greater than roughly 10, this will be very difficult to
    read and interpret.
    
    \return string representation of the suffix tree

  */
  std::string tostring() const;
  
  /*!
    \brief Get the length of the sequence(s) encoded in the suffix
    tree.
  */
  size_t seqlen() const {return sequence.length() - 1;}
  
  /*!
    \brief Get the top score, option to specify sequence.

    This function returns the score of the top matching substring
    relative to the specified scoring matrix.

    \param sm 
    the scoring matrix relative to which the score will be calculated

    \param sequence_id
    (optional) The index of the sequence in which to search. This defaults
    to searching all sequences
    
    \return the value of the top scoring match between the specified
    scoring matrix and any substring of the encoded sequence having
    the same width.
  */
  float top_score(const ScoringMatrix& sm, 
		  int sequence_id = ALL_SEQUENCES) const;
  /*!
    \brief Get the top score, within specified range, option to
    specify sequence.
    
    This function returns the score of the top matching substring
    relative to the specified scoring matrix. This function also has
    parameters to specify a range (a window) in which to search within
    the sequence.

    \param sm 
    the scoring matrix relative to which the score will be calculated
    \param min_offset
    Search for matches at positions greater than or equal to this.
    \param max_offset
    Search for matches at positions smaller than this.
    \param sequence_id
    (optional) The index of the sequence in which to search. This defaults
    to searching all sequences
    
    \return the value of the top scoring match between the specified
    scoring matrix and any substring of the encoded sequence having
    the same width.
  */
  float top_score(const ScoringMatrix& sm,
		  size_t min_offset, size_t max_offset,
		  int sequence_id = ALL_SEQUENCES) const;
  
  /*!
    \brief Get the top k scores, option to specify sequence

    \param sm 
    the scoring matrix relative to which the score will be calculated

    \param sequence_id
    (optional) The index of the sequence in which to search. This defaults
    to searching all sequences

  */
  void top_scores(std::vector<float>&, 
		  size_t, 
		  const ScoringMatrix& sm, 
		  int sequence_id = ALL_SEQUENCES) const;
  
  /*!
    \brief Get scores greater than a threshold, option to specify
    sequence

    \param sm 
    the scoring matrix relative to which the score will be calculated

    \param sequence_id
    (optional) The index of the sequence in which to search. This defaults
    to searching all sequences

  */
  void scores_greater(std::vector<float>&, float, 
		      const ScoringMatrix& sm,
		      int sequence_id = ALL_SEQUENCES) const;
  
  /*!
    \brief Get the top score location and score, option to specify
    sequence.

    \param sm 
    the scoring matrix relative to which the score will be calculated

    \param sequence_id
    (optional) The index of the sequence in which to search. This defaults
    to searching all sequences

    \return an aggregate of the location of the top match and the
    corresponding score
  */
  seq_pos_score top_score_index(const ScoringMatrix& sm, 
				int sequence_id = ALL_SEQUENCES) const;

  /*!  
    \brief Get location and score of top scoring match within a range,
    option to specify sequence.
    
    \param sm 
    the scoring matrix relative to which the score will be calculated
    \param min_offset
    Search for matches at positions greater than or equal to this.
    \param max_offset
    Search for matches at positions smaller than this.
    \param sequence_id
    (optional) The index of the sequence in which to search. This defaults
    to searching all sequences

  */
  seq_pos_score top_score_index(const ScoringMatrix& sm,
				size_t min_offset,
				size_t max_offset,
				int the_sequence_id = ALL_SEQUENCES) const;

  /*!
    \brief Get the top k scores with indices, option to specify sequence
    
    \param sm 
    the scoring matrix relative to which the score will be calculated
    
    \param sequence_id
    (optional) The index of the sequence in which to search. This defaults
    to searching all sequences
  */
  void top_scores_indices(std::vector<seq_pos_score>&,
			  size_t,
			  const ScoringMatrix &sm,
			  int sequence_id = ALL_SEQUENCES) const;
  
  /*!
    \brief Get scores and indices with scores greater than a cutoff,
    option to specify sequence

    \param sm 
    the scoring matrix relative to which the score will be calculated
    \param cutoff
    score above which to find matches
    \param sequence_id
    (optional) The index of the sequence in which to search. This defaults
    to searching all sequences

  */
  void scores_greater_indices(float cutoff, 
			      const ScoringMatrix &sm,
			      std::vector<seq_pos_score> &hits,
			      int sequence_id = ALL_SEQUENCES) const;
  
  // get scores (and indices) above threshold for sequence within a window
  void window_scores_greater_indices(std::vector<seq_pos_score>&, float, 
				     const ScoringMatrix& sm,
				     size_t, size_t, size_t) const;
  void window_scores_greater_indices(std::vector<seq_pos_score>&,
				     float,
				     const ScoringMatrix& sm,
				     size_t,
				     size_t,
				     size_t,
				     std::vector<std::pair<size_t, size_t> >&) const;
  void window_scores_greater_indices(std::vector<seq_pos_score>&,
				     float,
				     const ScoringMatrix&,
				     const ScoringMatrix&,
				     size_t,
				     size_t,
				     size_t) const;

  /*!
    \warning DEPRICATED
    
    return True of there are no overlaps in the specified sites
  */
  static bool no_overlaps(std::vector<std::pair<size_t, size_t> > &,
 			  size_t,
 			  size_t);
  
  void window_scores_greater(std::vector<float>&, float,
			     const ScoringMatrix&,
			     int,
			     int) const;
  
  /************************************************/
  // Get locations of substrings of width in a particular window,
  // and number of copies greater than a particular threshold
  void get_kmer_counts(std::vector<std::pair<std::string, size_t> > &,
		       size_t) const;
  size_t get_count(std::string w) const;
  void get_locations(std::string w,
		     std::vector<std::pair<size_t, size_t> > &L) const;
private:
  std::string sequence;
  size_t depth;
  SuffixNode* root;
  std::vector<size_t> offset;
  
  std::pair<size_t, size_t> index2seq_offset(size_t) const;
  static float** PrepareScoringMatrix(const ScoringMatrix&, float &);
  
  SuffixTree(const SuffixTree&);
  SuffixTree& operator=(const SuffixTree&);
  
  static const int FULL_DEPTH = -1;
  static const int ALL_SEQUENCES = -1;
};

#endif
