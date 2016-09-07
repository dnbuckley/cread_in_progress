/*
 * Copyright (C) 2006 Cold Spring Harbor Laboratory
 * Authors: Andrew D. Smith, Dustin Schones, Pavel Sumazin and Michael Q. Zhang
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

#include "cread.hpp"
#include "Motif.hpp"
#include "ScoringMatrix.hpp"
#include "IndexedSuffixTree.hpp"
#include "GeneralizedSuffixTree.hpp"
#include "smithlab_utils.hpp"
#include "FastaFile.hpp"
#include "WordTable.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::greater;
using std::priority_queue;
using std::min;
using std::cout;
using std::cerr;
using std::endl;
using std::numeric_limits;
using std::ofstream;
using std::ptr_fun;
using std::ostream;
using std::ostream_iterator;
using std::mem_fun_ref;
using std::pair;
using std::make_pair;
using std::copy;


string threshold_label;

/* This is the suffix of files that are supposed to be the FASTA files
 * in the directory of sequences to search (assuming the user gives a
 * directory). This variable also acts like a flag: if it is set, we
 * assume the user intends the '-s' value to be a
 * directory. Otherwise, the value specified with '-s' is a single
 * FASTA file.
 */ 
string fasta_suffix;

bool VERBOSE = false;
int buffer_size = 2500000;
int n_top = 0;
int core_size = 15;
float threshold_param = numeric_limits<float>::max();
bool functional_depth = false;
bool p_value = false;
bool handleties = false;
bool bysequence = false;
bool dme_matrix = false;
bool max_score = false;
bool no_preprocessing = false;
bool single_strand = false;

/* The threshold_tolerance value, defined below, is subtracted from
   thresholds to avoid problems due to rounding. Since, when
   thresholds are used, we are interested in occurrences scoring
   greater than _or_equal_to_ the threshold, it is important that we
   subtract some small amount in case the threshold has been rounded
   up slightly.  The tolerance is also required when trying to get top
   scores and allowing ties. Some rounding happens when calculating
   scores.  */
const float threshold_tolerance = 0.0001; 
const float base_comp_tolerance = 0.0001;

/* The "contains_valid_site(string, size_t)" function just tests
 * whether or not the sequence has a substring, without any 'N's, that
 * is long enough to be an occurrence of a motif with the specified
 * width.
 */
bool contains_valid_site(const string& sequence, const size_t width) {
  size_t count = width;
  for (size_t i = 0; i < sequence.length(); ++i)
    if (!valid_base(sequence[i]))
      count = width;
    else if (--count == 0) 
      return true;
  return false;
}

/* The "contains_valid_site(vector<string>, size_t)" makes sure at
 * least one of the strings in the vector contains a valid site (using
 * the "contains_valid_site(string, size_t) function).
 */
bool contains_valid_site(const vector<string>& sequences, const size_t width) {
  for (size_t i = 0; i < sequences.size(); ++i)
    if (contains_valid_site(sequences[i], width)) 
      return true;
  return false;
}

// TODO: get rid of this hit thing, and just use a "Site"
struct Hit {
  int seq, pos;
  float score;
  char strand;
  Hit(const size_t s, const size_t p, const float sc, const char st) :
    seq(s), pos(p), score(sc), strand(st) {}
  Hit(const int s, const pos_score& ps, const char st) : 
    seq(s), pos(ps.first), score(ps.second), strand(st) {}
  Hit(const pos_score& ps, char st) :
    pos(ps.first), score(ps.second), strand(st) {}
  bool operator>(const Hit &h) const {return score > h.score;}
  static bool seq_pos_order(const Hit& h1, const Hit& h2) {
    return h1.seq < h2.seq || (h1.seq == h2.seq && h1.pos < h2.pos);
  }
  string tostring() const {
    std::ostringstream ss;
    ss << seq << "\t" << pos << "\t" << score << "\t" << strand;
    return ss.str();
  }
};
typedef priority_queue<Hit, vector<Hit>, greater<Hit> > Hit_queue;


void AddDirection(const vector<pos_score>& ps, 
		  vector<Hit>& hits, const char dir) {
  for (size_t i = 0; i < ps.size(); ++i) 
    hits.push_back(Hit(ps[i], dir));
}


void AddDirection(const vector<seq_pos_score>& sps, 
		  vector<Hit>& hits, const char dir) {
  for (size_t i = 0; i < sps.size(); ++i) 
    hits.push_back(Hit(sps[i].first.first, sps[i].first.second, 
		       sps[i].second, dir));
}


void AddDirectionAndSequence(const vector<pos_score>& ps, 
			     vector<Hit>& hits, const char dir, const int sn) {
  for (size_t i = 0; i < ps.size(); i++)
    hits.push_back(Hit(sn, ps[i], dir));
}

// TODO: test these inlines to make sure they help
inline bool
valid_subsequence(const int* offset, const size_t matwidth) {
  size_t k = 0;
  while (k < matwidth && valid_base_id(offset[k]))
    ++k;
  return k == matwidth;
}

inline float
match_matrix(const ScoringMatrix &sm, const int* offset, const size_t width) {
  float score = 0;
  for (size_t i = 0; i < width; ++i)
    score += sm[i][offset[i]];
  return score;
}

void 
QuerySequenceNoPreprocessing(const string sequence, const vector<ScoringMatrix>& sm,
			     const vector<ScoringMatrix>& smrc,
			     vector<vector<Hit> >& occ,
			     const vector<float>& threshold, size_t sn) {
//  typedef pair<pair<size_t, size_t>, float> seq_pos_score;
//  typedef pair<size_t, size_t> seq_pos;
  vector<int> helper(sequence.length());
  transform(sequence.begin(), sequence.end(), 
	    helper.begin(), &base2int);
  if (occ.empty()) occ = vector<vector<Hit> >(sm.size());
  for (size_t i = 0; i < sm.size(); ++i) {
    const size_t matwidth = sm[i].get_width();
    const size_t lim = std::max(static_cast<int>(sequence.length()) - 
				static_cast<int>(matwidth) + 1, 0);
    for (size_t j = 0; j < lim; ++j) {
      const int* offset = &helper[j];
      if (valid_subsequence(offset, matwidth)) {
	float temp_score = match_matrix(sm[i], offset, matwidth);
	if (temp_score >= threshold[i])
	  occ[i].push_back(Hit(sn, pos_score(j, temp_score), 'p'));
	if (!single_strand) {
	  temp_score = match_matrix(smrc[i], offset, matwidth);
	  if (temp_score >= threshold[i])
	    occ[i].push_back(Hit(sn, pos_score(j, temp_score), 'n'));
	}
      }
    }
    sort(occ[i].begin(), occ[i].end(), greater<Hit>());
    sort(occ[i].begin(), occ[i].end(), ptr_fun(&Hit::seq_pos_order));
  }
}

void 
QuerySequenceNoPreprocessing(const string& sequence, 
			     const vector<ScoringMatrix>& sm, 
			     const vector<ScoringMatrix>& smrc, 
			     vector<vector<Hit> >& occ, 
			     size_t n_top, size_t sn) {
//  typedef pair<pair<size_t, size_t>, float> seq_pos_score;
  vector<int> helper(sequence.length());
  transform(sequence.begin(), sequence.end(), 
	    helper.begin(), &base2int);
  for (size_t i = 0; i < sm.size(); ++i) {
    
    const size_t matwidth = sm[i].get_width();
    float cutoff = -numeric_limits<float>::max();
    Hit_queue hq;
    const size_t lim = std::max(static_cast<int>(sequence.length()) - 
				static_cast<int>(matwidth) + 1, 0);
    for (size_t j = 0; j < lim; ++j) {
      const int* offset = &helper[j];
      if (valid_subsequence(offset, matwidth)) {
	float temp_score = match_matrix(sm[i], offset, matwidth);
	if (temp_score >= cutoff) {
	  if (hq.size() == n_top) 
	    hq.pop();
	  hq.push(Hit(sn, j, temp_score, 'p'));
	  if (hq.size() == n_top) 
	    cutoff = hq.top().score;
	}
	if (!single_strand) {
	  temp_score = match_matrix(smrc[i], offset, matwidth);
	  if (temp_score >= cutoff) {
	    if (hq.size() == n_top) hq.pop();
	    hq.push(Hit(sn, j, temp_score, 'n'));
	    if (hq.size() == n_top) cutoff = hq.top().score;
	  }
	}
      }
    }
    while (hq.size() > 0) {
      occ[i].push_back(hq.top());
      hq.pop();
    }
    sort(occ[i].begin(), occ[i].end(), greater<Hit>());
    sort(occ[i].begin(), occ[i].end(), ptr_fun(&Hit::seq_pos_order));
  }
}

typedef priority_queue<float, vector<float>, greater<float> > greater_pq;
inline void
update_tied_queue(float value, greater_pq &q, size_t max_size) {
  if (value >= q.top() || q.size() < max_size) {
    if (q.size() == max_size) q.pop();
    q.push(value);
  }
}


void 
QuerySequenceSetNoPreprocessingTies(const vector<string>& sequences,
				    const vector<ScoringMatrix>& sm,
				    const vector<ScoringMatrix>& smrc,
				    vector<vector<Hit> >& occ, size_t n_top) {
  vector<vector<int> > helpers(sequences.size());
  for (size_t i = 0; i < sequences.size(); ++i)
    transform(sequences[i].begin(), sequences[i].end(),
	      back_inserter(helpers[i]), &base2int);
  vector<greater_pq> sq(sm.size());
  for (size_t s = 0; s < sequences.size(); ++s)
    for (size_t i = 0; i < sm.size(); ++i) {
      const size_t matwidth = sm[i].get_width();
      const size_t lim = std::max(static_cast<int>(sequences[s].length()) - 
				  static_cast<int>(matwidth) + 1, 0);
      for (size_t j = 0; j < lim; ++j) {
	const int* offset = &helpers[s][j];
	if (valid_subsequence(offset, matwidth)) {
	  float temp_score = match_matrix(sm[i], offset, matwidth);
	  update_tied_queue(temp_score, sq[i], n_top);
	  if (!single_strand) {
	    temp_score = match_matrix(smrc[i], offset, matwidth);
	    update_tied_queue(temp_score, sq[i], n_top);
	  }
	}
      }
    }
  vector<float> cutoff;
  for (size_t i = 0; i < sq.size(); ++i)
    cutoff.push_back(sq[i].top() - threshold_tolerance);
  for (size_t s = 0; s < sequences.size(); ++s)
    for (size_t i = 0; i < sm.size(); ++i) {
      const size_t matwidth = sm[i].get_width();
      const size_t lim = std::max(static_cast<int>(sequences[s].length()) -
				  static_cast<int>(matwidth) + 1, 0);
      for (size_t j = 0; j < lim; ++j) {
	const int* offset = &helpers[s][j];
	if (valid_subsequence(offset, matwidth)) {
	  float temp_score = match_matrix(sm[i], offset, matwidth);
	  if (temp_score >= cutoff[i])
	    occ[i].push_back(Hit(s, j, temp_score, 'p'));
	  if (!single_strand) {
	    temp_score = match_matrix(smrc[i], offset, matwidth);
	    if (temp_score >= cutoff[i])
	      occ[i].push_back(Hit(s, j, temp_score, 'n'));
	  }
	}
      }
    }
  for (size_t i = 0; i < sm.size(); ++i) {
    sort(occ[i].begin(), occ[i].end(), greater<Hit>());
    sort(occ[i].begin(), occ[i].end(), ptr_fun(&Hit::seq_pos_order));
  }
}


void 
QuerySequenceNoPreprocessingTies(const string& sequence, 
				 const vector<ScoringMatrix>& sm,
				 const vector<ScoringMatrix>& smrc,
				 vector<vector<Hit> >& occ, 
				 size_t n_top, size_t sn) {
  vector<int> helper(sequence.length());
  transform(sequence.begin(), sequence.end(),
	    helper.begin(), &base2int);
  for (size_t i = 0; i < sm.size(); ++i) {
    size_t matwidth = sm[i].get_width();
    // Find the k-th top score
    const size_t lim = std::max(static_cast<int>(sequence.length()) -
				static_cast<int>(matwidth) + 1, 0);
    greater_pq sq;
    for (size_t j = 0; j < lim; ++j) {
      const int* offset = &helper[j];
      if (valid_subsequence(offset, matwidth)) {
	float temp_score = match_matrix(sm[i], offset, matwidth);
	update_tied_queue(temp_score, sq, n_top);
	if (!single_strand) {
	  temp_score = match_matrix(smrc[i], offset, matwidth);
	  update_tied_queue(temp_score, sq, n_top);
	}
      }
    }
    float cutoff = sq.top() - threshold_tolerance;
    for (size_t j = 0; j < lim; ++j) {
      const int* offset = &helper[j];
      if (valid_subsequence(offset, matwidth)) {
	float temp_score = match_matrix(sm[i], offset, matwidth);
	if (temp_score >= cutoff)
	  occ[i].push_back(Hit(sn, j, temp_score, 'p'));
	if (!single_strand) {
	  temp_score = match_matrix(smrc[i], offset, matwidth);
	  if (temp_score >= cutoff)
	    occ[i].push_back(Hit(sn, j, temp_score, 'n'));
	}
      }
    }
    sort(occ[i].begin(), occ[i].end(), greater<Hit>());
    sort(occ[i].begin(), occ[i].end(), ptr_fun(&Hit::seq_pos_order));
  }
}



void 
QuerySequenceSetNoPreprocessing(const vector<string>& sequences, 
				const vector<ScoringMatrix>& sm,
				const vector<ScoringMatrix>& smrc, 
				vector<vector<Hit> >& occ, size_t n_top) {
  vector<Hit_queue> hq(sm.size());
  vector<float> cutoff(sm.size());
  fill_n(cutoff.begin(), sm.size(), -numeric_limits<float>::max());
  for (size_t s = 0; s < sequences.size(); ++s) {
    vector<int> helper(sequences[s].length());
    transform(sequences[s].begin(), sequences[s].end(), 
	      helper.begin(), &base2int);
    for (size_t i = 0; i < sm.size(); ++i) {
      const size_t matwidth = sm[i].get_width();
      const size_t lim = std::max(static_cast<int>(sequences[s].length()) - 
				  static_cast<int>(matwidth) + 1, 0);
      for (size_t j = 0; j < lim; ++j) {
	const int* offset = &helper[j];
	if (valid_subsequence(offset, matwidth)) {
	  float temp_score = match_matrix(sm[i], offset, matwidth);
	  if (temp_score >= cutoff[i]) {
	    if (hq[i].size() == n_top)
	      hq[i].pop();
	    hq[i].push(Hit(s, j, temp_score, 'p'));
	    if (hq[i].size() == n_top)
	      cutoff[i] = hq[i].top().score;
	  }
	  if (!single_strand) {
	    temp_score = match_matrix(smrc[i], offset, matwidth);
	    if (temp_score >= cutoff[i]) {
	      if (hq[i].size() == n_top) hq[i].pop();
	      hq[i].push(Hit(s, j, temp_score, 'n'));
	      if (hq[i].size() == n_top) cutoff[i] = hq[i].top().score;
	    }
	  }
	}
      }
    }
  }
  for (size_t i = 0; i < sm.size(); ++i) {
    while (hq[i].size() > 0) {
      occ[i].push_back(hq[i].top());
      hq[i].pop();
    }
    sort(occ[i].begin(), occ[i].end(), greater<Hit>());
    sort(occ[i].begin(), occ[i].end(), ptr_fun(&Hit::seq_pos_order));
  }
}


void 
QuerySequenceByThreshold(const string& sequence,
			 const vector<ScoringMatrix>& sm,
			 const vector<ScoringMatrix>& smrc,
			 vector<vector<Hit> >& occ,
			 vector<float>& threshold, size_t sn) {
  if (no_preprocessing)
    QuerySequenceNoPreprocessing(sequence, sm, smrc, occ, threshold, sn);
  else {
    IndexedSuffixTree tree(sequence);
    for (size_t i = 0; i < sm.size(); ++i) {
      if (contains_valid_site(sequence, sm[i].get_width())) {
	vector<pos_score> temp_hits;
	tree.scores_greater_indices(sm[i], threshold[i], temp_hits);
	vector<Hit> hits;
	AddDirectionAndSequence(temp_hits, hits, 'p', sn);
	if (!single_strand) {
	  temp_hits.clear();
	  tree.scores_greater_indices(smrc[i], threshold[i], temp_hits);
	  AddDirectionAndSequence(temp_hits, hits, 'n', sn);
	}
	sort(hits.begin(), hits.end(), ptr_fun(&Hit::seq_pos_order));
	copy(hits.begin(), hits.end(), back_inserter(occ[i]));
      }
    }
  }
}

/* The parameter "sn" is an index for the sequence (sequence name)
 * which allows us to map the occurrences back to the sequences.
 * For reasons of efficiency, SuffixTree doesn't store the names,
 * just numbers.
 */
void
QuerySequenceByCount(const string& sequence,
		     const vector<ScoringMatrix>& sm,
		     const vector<ScoringMatrix>& smrc,
		     vector<vector<Hit> >& occ,
		     size_t n_top, size_t sn) {
  if (no_preprocessing)
    if (handleties)
      QuerySequenceNoPreprocessingTies(sequence, sm, smrc, occ, n_top, sn);
    else QuerySequenceNoPreprocessing(sequence, sm, smrc, occ, n_top, sn);
  else {
    IndexedSuffixTree tree(sequence);
    for (size_t i = 0; i < sm.size(); i++) {
      if (contains_valid_site(sequence, sm[i].get_width())) {
	vector<pos_score> temp_hits;
	
	tree.top_scores_indices(sm[i], n_top, temp_hits);
	vector<Hit> hits;
	AddDirectionAndSequence(temp_hits, hits, 'p', sn);
	
	if (!single_strand) {
	  temp_hits.clear();
	  tree.top_scores_indices(smrc[i], n_top, temp_hits);
	  AddDirectionAndSequence(temp_hits, hits, 'n', sn);
	}
      
	sort(hits.begin(), hits.end(), greater<Hit>());

	if (handleties) {
	  size_t index = min(n_top - 1, hits.size() - 1);
	  float threshold = hits[index].score - threshold_tolerance;
	
	  hits.clear();
	  temp_hits.clear();
	  tree.scores_greater_indices(sm[i], threshold, temp_hits);
	  AddDirectionAndSequence(temp_hits, hits, 'p', sn);
	
	  if (!single_strand) {
	    temp_hits.clear();
	    tree.scores_greater_indices(smrc[i], threshold, temp_hits);
	    AddDirectionAndSequence(temp_hits, hits, 'n', sn);
	  }
	  sort(hits.begin(), hits.end(), greater<Hit>());
	  sort(hits.begin(), hits.end(), ptr_fun(&Hit::seq_pos_order));
	  copy(hits.begin(), hits.end(), back_inserter(occ[i]));
	}
	else {
	  sort(hits.begin(), hits.begin() + min(n_top, hits.size()), 
	       ptr_fun(&Hit::seq_pos_order));
	  copy(hits.begin(), hits.begin() + min(n_top, hits.size()),
	       back_inserter(occ[i]));
	}
      }
    }
  }
}



void 
QuerySequenceSetByCount(const vector<string>& sequences,
			const vector<ScoringMatrix>& sm,
			const vector<ScoringMatrix>& smrc,
			vector<vector<Hit> >& occ,
			size_t n_top) {
  if (no_preprocessing)
    if (handleties)
      QuerySequenceSetNoPreprocessingTies(sequences, sm, smrc, occ, n_top);
    else QuerySequenceSetNoPreprocessing(sequences, sm, smrc, occ, n_top);
  else {
    GeneralizedSuffixTree tree(sequences);
    for (size_t i = 0; i < sm.size(); ++i) {
      if (contains_valid_site(sequences, sm[i].get_width())) {
	vector<seq_pos_score> temp_hits;
	tree.top_scores_indices(temp_hits, n_top, sm[i]);
	
	vector<Hit> hits;
	AddDirection(temp_hits, hits, 'p');
	if (!single_strand) {
	  temp_hits.clear();
	  tree.top_scores_indices(temp_hits, n_top, smrc[i]);
	  AddDirection(temp_hits, hits, 'n');
	}
	sort(hits.begin(), hits.end(), greater<Hit>());
	if (handleties) {
	  const float threshold = ((n_top < hits.size()) ? 
				   hits[n_top - 1].score :
				   hits.back().score) - threshold_tolerance;
	  hits.clear();
	  temp_hits.clear();
	  tree.scores_greater_indices(temp_hits, threshold, sm[i]);
	  AddDirection(temp_hits, hits, 'p');
	  if (!single_strand) {
	    temp_hits.clear();
	    tree.scores_greater_indices(temp_hits, threshold, smrc[i]);
	    AddDirection(temp_hits, hits, 'n');
	  }
	  sort(hits.begin(), hits.end(), greater<Hit>());
	  sort(hits.begin(), hits.end(), ptr_fun(&Hit::seq_pos_order));
	  copy(hits.begin(), hits.end(), back_inserter(occ[i]));
	}
	else {
	  sort(hits.begin(), hits.begin() + min(n_top, hits.size()), 
	       ptr_fun(&Hit::seq_pos_order));
	  copy(hits.begin(), hits.begin() + min(n_top, hits.size()),
	       back_inserter(occ[i]));
	}
      }
    }
  }
}


void 
get_motif_thresholds(vector<float>& threshold, vector<Motif>& motifs) {
  if (threshold_param == numeric_limits<float>::max() && threshold_label.c_str())
    for (vector<Motif>::iterator i = motifs.begin(); i != motifs.end(); ++i) {
      string threshold_str(i->get_attribute(threshold_label.c_str()));
      if (!threshold_str.empty())
	threshold.push_back(atof(threshold_str.c_str()));
      // TODO: better handle the case where threshold label doesn't exist
      else {
	cerr << "ERROR: no attribute with label: " 
	     << threshold_label << " to be used as threshold." << endl;
	exit(EXIT_FAILURE);
      }
    }
  else std::fill_n(back_inserter(threshold), motifs.size(), threshold_param);
}

void 
get_scoring_matrices(const vector<Matrix> &matrices, 
		     const vector<float>& base_comp,
		     vector<ScoringMatrix> &sm,
		     vector<ScoringMatrix> &smrc, 
		     size_t &max_width) {
  max_width = 0;
  vector<Matrix>::const_iterator i(matrices.begin());
  for (; i != matrices.end(); ++i) {
    max_width = std::max(max_width, i->get_width());
    if (!i->is_count_mat() || dme_matrix)
      sm.push_back(ScoringMatrix(*i, base_comp));
    else 
      sm.push_back(ScoringMatrix::StormoScoringMatrix(*i, base_comp));
    if (!single_strand) 
      smrc.push_back(sm.back().revcomp());
  }
}



void 
parse_base_comp_str(const char *bc_str, vector<float>& base_comp) {
  vector<string> parts(cread::split(bc_str, ","));
  if (parts.size() != alphabet_size) {
    cerr << "ERROR: incorrect base composition: " << bc_str << endl;
    exit(EXIT_FAILURE);
  }
  float total = 0;
  for (size_t i = 0; i < alphabet_size; ++i) {
    base_comp[i] = atof(parts[i].c_str());
    total += base_comp[i];
  }
  if (std::abs(1.0 - total) > base_comp_tolerance) {
    cerr << "ERROR: base composition must sums to 1.0, not " << total << endl;
    exit(EXIT_FAILURE);
  }
}


struct NameScoreIndex {
  string name;
  float score;
  size_t index;
  NameScoreIndex() : score(0), index(0) {}
  NameScoreIndex(string n, float s, size_t i) : 
    name(n), score(s), index(i) {}
  bool operator>(const NameScoreIndex& rhs) const {
    return ((name > rhs.name) || (score > rhs.score && name == rhs.name));
  }
};

void
determine_which_to_keep_with_ties(const vector<NameScoreIndex>& nm_scr_idx,
				  const size_t k_to_keep, vector<bool>& keep) {
  size_t in_same_seq = 0;
  bool kept_prev = false;
  for (size_t j = 0; j < nm_scr_idx.size(); ++j) {
    if (j == 0 || nm_scr_idx[j].name != nm_scr_idx[j - 1].name)
      in_same_seq = 0;
    else ++in_same_seq;
    if (in_same_seq < k_to_keep || 
	(kept_prev && nm_scr_idx[j].score == nm_scr_idx[j - 1].score)) {
      keep[nm_scr_idx[j].index] = true;
      kept_prev = true;
    }
    else kept_prev = false;
  }
}


void
determine_which_to_keep(const vector<NameScoreIndex>& nm_scr_idx, 
			const size_t k_to_keep, vector<bool>& keep) {
  size_t in_same_seq = 0;
  for (size_t i = 0; i < nm_scr_idx.size(); ++i) {
    if (i == 0 || nm_scr_idx[i].name != nm_scr_idx[i - 1].name)
      in_same_seq = 0;
    else ++in_same_seq;
    if (in_same_seq < k_to_keep)
      keep[nm_scr_idx[i].index] = true;
  }
}


/* Remove all but the top k sites for a motif in each sequence. */
void
top_k_per_sequence(vector<vector<MotifSite> >& sites,
		   const size_t k_to_keep, const bool allow_ties) {
  // iterate over each motif
  for (size_t i = 0; i < sites.size(); ++i) {
    // get the hits for each sequence and sort them
    vector<NameScoreIndex> nm_scr_idx(sites[i].size());
    for (size_t j = 0; j < sites[i].size(); ++j)
      nm_scr_idx[j] = NameScoreIndex(sites[i][j].get_seq_name(),
				     sites[i][j].get_score(), j);
    sort(nm_scr_idx.begin(), nm_scr_idx.end(), 
	 greater<NameScoreIndex>());
    
    // vector of indices of sites that we will keep
    vector<bool> keep(sites[i].size(), false);
    if (allow_ties)
      determine_which_to_keep_with_ties(nm_scr_idx, k_to_keep, keep);
    else determine_which_to_keep(nm_scr_idx, k_to_keep, keep);

    vector<MotifSite> temp_sites;
    for (size_t j = 0; j < keep.size(); ++j)
      if (keep[j])
	temp_sites.push_back(sites[i][j]);
    sites[i].swap(temp_sites);
  }
}


/* Function to remove all but the top k sites for a motif in entire
 * sequence set
 */
struct ScoreIndex {
  float score;
  size_t index;
  ScoreIndex(float s = 0, size_t i = 0) : score(s), index(i) {}
  bool operator>(const ScoreIndex& rhs) const {return score > rhs.score;}
};
void
top_k_overall(vector<vector<MotifSite> >& sites,
	      const size_t k_to_keep, const bool allow_ties) {
  // iterate over each motif
  for (size_t i = 0; i < sites.size(); ++i) {
    // get the hits for each sequence and sort them
    vector<ScoreIndex> scr_idx(sites[i].size());
    for (size_t j = 0; j < sites[i].size(); ++j)
      scr_idx[j] = ScoreIndex(sites[i][j].get_score(), j);
    sort(scr_idx.begin(), scr_idx.end(), greater<ScoreIndex>());
    // get a new vector including only the top k hits
    vector<MotifSite> temp_sites;
    const size_t loop_lim = min(k_to_keep, scr_idx.size());
    for (size_t j = 0; j < loop_lim || 
	   (allow_ties && 
	    scr_idx[j].score == scr_idx[j - 1].score && 
	    j < scr_idx.size()); ++j)
      temp_sites.push_back(sites[i][scr_idx[j].index]);
    sites[i].swap(temp_sites);
  }
}


struct MotifSiteSeqnameLess {
  bool operator()(const MotifSite &lhs, const MotifSite &rhs) {
    return lhs.get_seq_name() < rhs.get_seq_name();
  }
};


struct CoreScoreMat {
  Matrix mat;
  ScoringMatrix sm;
  size_t upstream;
  size_t downstream;
  CoreScoreMat(const Matrix& mat, const ScoringMatrix& sm, 
	       const size_t core_size, const WordTable &wt);
  void convert_score(vector<vector<MotifSite> >& sites,
		     const vector<ScoringMatrix>& sm,
		     const bool convert_to_fd,
		     const bool convert_to_pval,
		     const WordTable& wt) const;
};


CoreScoreMat::CoreScoreMat(const Matrix& matrix,
			   const ScoringMatrix& scoremat,
			   const size_t core_size, 
			   const WordTable &wt) {
  vector<float> info_profile(wt.info_profile(matrix));
  float max_info = 0.0;
  for (size_t i = 0; i <= matrix.get_width() - core_size; ++i) {
    const float sum = accumulate(info_profile.begin() + i,
				 info_profile.begin() + i + core_size, 0.0);
    if (sum > max_info) {
      max_info = sum;
      upstream = i;
    }
  }
  downstream = matrix.get_width() - (core_size + upstream);
  float *newmat[core_size], *newsm[core_size];
  for (size_t i = 0; i < core_size; ++i) {
    newmat[i] = new float[alphabet_size];
    std::copy(matrix[i + upstream], matrix[i + upstream] + alphabet_size, newmat[i]);
    newsm[i] = new float[alphabet_size];
    std::copy(scoremat[i + upstream], scoremat[i + upstream] + alphabet_size, newsm[i]);
  }
  mat = Matrix(newmat, core_size);
  sm = ScoringMatrix(newsm, core_size);
}


void
convert_scores_to_fds(vector<vector<MotifSite> >& sites,
		      const vector<ScoringMatrix>& sm) {
  for (size_t i = 0; i < sites.size(); ++i)
    for (size_t j = 0; j < sites[i].size(); ++j)
      sites[i][j].set_score(sm[i].functional_depth(sites[i][j].get_score()));
}

void 
update_motif_sites(const vector<vector<Hit> >& occ,
		   const vector<string>& sequences, 
		   const vector<string>& seqnames,
		   const vector<size_t>& full_motif_widths,
		   const size_t first_seq_offset,
		   vector<vector<MotifSite> >& sites) {
  typedef vector<Hit>::const_iterator hititer;
  for (size_t i = 0; i < sites.size(); ++i) {
    for (hititer j = occ[i].begin(); j != occ[i].end(); ++j) {
      const string temp(sequences[j->seq].substr(j->pos, full_motif_widths[i]));
      const string siteseq((j->strand == 'p') ? temp :
			   revcomp(temp));
      const size_t true_offset = j->pos + ((j->seq == 0) ?
					   first_seq_offset : 0);
      sites[i].push_back(MotifSite(siteseq, seqnames[j->seq],
				   true_offset, full_motif_widths[i], "",
				   j->strand, j->score));
    }
  }
}


/* Right now the function below takes waaay too long because each site
 * has to have a separate calculation. This needs to be fixed, and
 * also must be able to handle the Staden p-values, which would be
 * much quicker and easier, and should be done first.
 */
// void
// convert_scores_to_pvals(vector<vector<MotifSite> >& sites,
//                         const WordTable& wt,
//                         const vector<CoreScoreMat>& cores) {
//   for (size_t i = 0; i < sites.size(); ++i) {
//     if (VERBOSE) {
//       cerr << i + 1 << "/" << sites.size() << endl;
//     }
//     size_t hits = 0, totals = 0;
//     std::map<float, float> cached_pvals;
//     for (size_t j = 0; j < sites[i].size(); ++j) {
//       std::map<float, float>::iterator k(cached_pvals.find(sites[i][j].get_score()));
//       if (k == cached_pvals.end()) {
//         float pval = wt.cutoff2pval(sites[i][j].get_score(),
//                                     cores[i].mat, cores[i].sm,
//                                     !single_strand);
//         cached_pvals[sites[i][j].get_score()] = pval;
//         sites[i][j].set_score(pval);
//       }
//       else {
//         sites[i][j].set_score(k->second);
//         hits++;
//       }
//       totals++;
//     }
//   }
// }


void 
update_motif_sites(const vector<vector<Hit> >& occ,
		   const vector<string>& sequences, 
		   const vector<string>& seqnames,
		   const vector<size_t>& widths,
		   const size_t first_seq_offset,
		   const vector<CoreScoreMat>& cores,
		   vector<vector<MotifSite> >& sites) {
  typedef vector<Hit>::const_iterator hititer;
  for (size_t i = 0; i < sites.size(); ++i) {
    const size_t downstream = cores[i].downstream;
    const size_t upstream = cores[i].upstream;
    for (hititer j = occ[i].begin(); j != occ[i].end(); ++j) {
      const int position = j->pos - ((j->strand == 'p') ?
				     upstream : downstream);
      if (position >= 0 && position + widths[i] - 1 < sequences[j->seq].length()) {
	const string temp(sequences[j->seq].substr(position, widths[i]));
	const string siteseq((j->strand == 'p') ? temp :
			     revcomp(temp));
	size_t true_offset = position + ((j->seq == 0) ?
					 first_seq_offset : 0);
	sites[i].push_back(MotifSite(siteseq, seqnames[j->seq],
				     true_offset, widths[i], "",
				     j->strand, j->score));
      }
    }
  }
}


/* Functions to remove duplicate sites, which can happen (rarely) if
 * motifs are short enough and have matches in parts that are shared
 * between consecutive chunks of the FASTA file.
 */
MotifSite* get_ptr(MotifSite& s) {return &s;}
MotifSite get_obj(MotifSite* s) {return *s;}
struct MotifSitePtrGreaterScore {
  bool operator()(const MotifSite* lhs, const MotifSite* rhs) {
    return lhs->get_score() > rhs->get_score();}
};
struct MotifSitePtrLess {
  bool operator()(const MotifSite* lhs, const MotifSite* rhs) {
    return (*lhs) < (*rhs);}
};
struct MotifSitePtrEqual {
  bool operator()(const MotifSite* lhs, const MotifSite* rhs) {
    return (*lhs) == (*rhs);}
};
void
remove_duplicate_sites(const string progress_prefix, 
		       vector<vector<MotifSite> >& sites) {
  for (size_t i = 0; i < sites.size(); ++i) {
    if (VERBOSE) 
      cerr << "\r" << progress_prefix << "\t" 
	   << static_cast<size_t>((100.0*i)/sites.size()) << "%";
    vector<MotifSite*> site_ptrs(sites[i].size());
    transform(sites[i].begin(), sites[i].end(), 
	      site_ptrs.begin(), ptr_fun(get_ptr));
    sort(site_ptrs.begin(), site_ptrs.end(), MotifSitePtrGreaterScore());
    sort(site_ptrs.begin(), site_ptrs.end(), MotifSitePtrLess());
    vector<MotifSite*>::iterator j = unique(site_ptrs.begin(), 
					    site_ptrs.end(),
					    MotifSitePtrEqual());
    site_ptrs.erase(j, site_ptrs.end());
    vector<MotifSite> site_objs;
    transform(site_ptrs.begin(), site_ptrs.end(),
	      back_inserter(site_objs), ptr_fun(get_obj));
    sites[i].swap(site_objs);
  }
  if (VERBOSE) 
    cerr << "\r" << progress_prefix << "\t100%" << endl;
}



size_t 
get_total_file_sizes(const vector<string> &seqfiles) {
  size_t total = 0;
  for (size_t i = 0; i < seqfiles.size(); ++i) {
    std::ifstream f(seqfiles[i].c_str());
    size_t begin_pos = 0, end_pos = 0;
    if (f.good()) {
      begin_pos = f.tellg();
      f.seekg(0, std::ios_base::end);
      end_pos = f.tellg();
      f.close();
    }
    total += (end_pos - begin_pos);
  }
  return total;
}



string
file_progress_string(const string& filename,
		     const size_t file_id, const size_t max_file_id,
		     const size_t file_offset, const size_t max_file_offset) {
  std::ostringstream ss;
  ss << "\r" << "file " << filename << "\t"
     << static_cast<size_t>((100.0*file_offset)/max_file_offset) << "%\t"
     << "(" << static_cast<size_t>((100.0*file_id)/max_file_id) << "% total)";
  return ss.str();
}

int main(int argc, const char **argv) {
  
/* INPUT PARAMETERS */
  string seqfile;           // file containing foreground
       				     // sequences
  string outfile;           // file in which to print output
  string motif_file;        // file containing the motifs
  string base_comp_str;     // string of base composition
  string word_table;
  
  try {
  
    /***************** GET COMMAND LINE ARGUMENTS *******************/
    OptionParser opt_parse(strip_path(argv[0]),
                "evaluates cross-species conervation of motif occurrences",
			              "[Motif]");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                       false, outfile);
    opt_parse.add_opt("fasta-suffix", '\0', "Suffix for FASTA files;"
                      " assumes input is directory", false, fasta_suffix);
    opt_parse.add_opt("verbose", 'v', "Print extra info to standard error",
                       false, VERBOSE);
    opt_parse.add_opt("buffer-size", '\0', "Size for input buffer",
                       false, buffer_size);
    opt_parse.add_opt("hit-count", 'n', "Find this many top occurrences",
                       false, n_top);
    opt_parse.add_opt("sequences", 's', "File of sequences in which to search",
                       false, seqfile);
    opt_parse.add_opt("by-sequence", 'q', "Top matchs are per sequences (needs -n)",
                       false, bysequence);
    opt_parse.add_opt("handle-ties", 'h', "Report equivalent hits in case of ties",
                       false, handleties);
    opt_parse.add_opt("ax-score", 'm', "Alias for the parameters: \"-hq -n 1\"",
                       false, max_score);
    opt_parse.add_opt("dme-scoring", 'd', "Use scoring matrix used by dme program",
                       false, dme_matrix);
    opt_parse.add_opt("no-preprocessing", 'N', "No suffix tree preprocessing of sequences",
                       false, no_preprocessing);
    opt_parse.add_opt("single-strand", 'S', "Search only one strand",
                       false, single_strand);
    opt_parse.add_opt("base-comp", 'C', "comma separated base freqs (order:A,C,G,T)",
                       false, base_comp_str);
    opt_parse.add_opt("threshold", 't', "threshold score defining occurrences",
                       false, threshold_param);
    opt_parse.add_opt("func-depth", 'f',"thresholds given as functional depths",
                       false, functional_depth);
    opt_parse.add_opt("p-value", 'p', "threshold given as p-values",
                       false, p_value);
    opt_parse.add_opt("core-size", 'c', "core size for calculating p-values",
                       false, core_size);
    opt_parse.add_opt("label", 'l', "thresholds specifieds as motif attribute",
                       false, threshold_label);
    opt_parse.add_opt("word-table", 'H', "word-table file", false, word_table);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const vector<string> input_filenames(leftover_args);
    // TODO: add a test somewhere to determine if the functional depths
    // specified inside a motif file are reasonable.
    /**********************************************************************/
    motif_file = input_filenames[0];

    // get the names of all the sequence files
    vector<string> seqfiles;
    if (!fasta_suffix.empty()) {
      seqfiles = FastaFile::read_seqs_dir(seqfile.c_str(), fasta_suffix.c_str());
      if (VERBOSE)
	cerr << "found " << seqfiles.size() 
	     << " sequence files in " << seqfile.c_str() << endl;
    }
    else seqfiles.push_back(seqfile);

    /* Get the base composition, either from the sequences supplied,
     * the command line, or the WordTable file. If WordTable is used,
     * this is where we read it in.
     */
    vector<float> base_comp(alphabet_size);

    WordTable wt;
    if (!word_table.empty()) { // base comp from WordTable (and read in hit table)
      wt = WordTable::ReadWordTable(word_table.c_str());
      vector<float> bc = wt.extract_base_comp();
      copy(bc.begin(), bc.end(), base_comp.begin());
    }
    else{ // no WordTable
      if (!base_comp_str.empty()) // base comp from command line
	parse_base_comp_str(base_comp_str.c_str(), base_comp);
      else  // base comp from the input file(s)
	FastaFile::base_comp_from_files(seqfiles, base_comp);
      if (VERBOSE) {
	cerr << "base comp: ";
	copy(base_comp.begin(), base_comp.end() - 1,
	     ostream_iterator<float>(cerr, ","));
	cerr << base_comp.back() << endl;
      }
    }

    /* Now we deal with reading in the motifs, turning them into
     * scoring matrices, and obtaining their thresholds.  If the 'top
     * k' matches for a motif is what was requested, then nothing is
     * done about the thresholds. If thresholds are needed, they are
     * specified by:
     * (1) Supplied on the command line (optionally with the '-f' flag
     * to indicate functional depth
     * (2) Specified as an attribute of the motif inside the motifs
     * file, in which case the '-l' argument will specify the 'key'
     * whose value should be the threshold.
     *
     * The thresholds can be specified as an absolute score value, a
     * functional depth value, or a p-value. For p-values, a WordTable
     * must also be given.
     */
    
    // read in the motifs
    vector<Motif> motifs(Motif::ReadMotifVector(motif_file.c_str()));
    
    // get the score cutoff associated with each motif
    vector<float> threshold;
    if (threshold_param != numeric_limits<float>::max() || !threshold_label.empty())
      get_motif_thresholds(threshold, motifs);

    // extract the matrix from each motif
    vector<Matrix> matrices;
    transform(motifs.begin(), motifs.end(), back_inserter(matrices),
	      mem_fun_ref(&Motif::get_matrix));

    // get a vector of motif widths
    vector<size_t> full_motif_widths;
    transform(motifs.begin(), motifs.end(), back_inserter(full_motif_widths),
	      mem_fun_ref(&Motif::get_width));
    
    // make vector of scoring matrices (and their reverse complements
    // if needed).
    vector<ScoringMatrix> sm, smrc;
    size_t max_width;
    get_scoring_matrices(matrices, base_comp, sm, smrc, max_width);

    // CORES STUFF
    /* If we are calculating p-values, then sometimes we can't do that
     * for very long motifs, and so we get some shorter 'cores', which
     * have a default width of 15
     */
    vector<CoreScoreMat> cores;
    if (!word_table.empty()) {
      for (size_t i = 0; i < matrices.size(); ++i) {
	cores.push_back(CoreScoreMat(matrices[i], sm[i],
				     min(static_cast<size_t>(core_size),
					 matrices[i].get_width()), wt));
	sm[i] = cores[i].sm;
	smrc[i] = sm[i].revcomp();
      }
    }
    
    StadenPValue *spv = (p_value && word_table.empty()) ?
      new StadenPValue(base_comp, 1000) : 0;
    
    if (!threshold.empty()) {
      // if thresholds were specified as functional depths, get
      // corresponding scores
      if (functional_depth) 
	for (size_t i = 0; i < sm.size(); ++i)
	  threshold[i] = sm[i].functional_depth_to_score(threshold[i]);
      else 
	// if they were specified as p-values, get the corresponding
	// scores
	if (p_value) {
	  for (size_t i = 0; i < sm.size(); ++i) {
	    if (VERBOSE)
	      cerr << "\r" << "calculating score cutoff\t"
		   << static_cast<size_t>((100.0*i)/sm.size()) << "%";
	    if (!word_table.c_str())
	      threshold[i] = wt.pval2cutoff(threshold[i], cores[i].mat,
					    cores[i].sm, !single_strand);
	    else{
	      threshold[i] = spv->get_score(sm[i], threshold[i]);
            }
	  }
	  if (VERBOSE)
	    cerr << "\r" << "calculating score cutoff\t100%" << endl;
	}
      // correct the threshold to overcome roundoff errors... The
      // only thing more annoying than having to do this would be
      // taking the time to figure out the proper way to fix it
      if (!threshold.empty())
	transform(threshold.begin(), threshold.end(), threshold.begin(),
		  bind2nd(std::minus<float>(), threshold_tolerance));
    }
    
    // erase all sites in the motifs -- we only want the new ones
    // and memory will grow, so lets delete them now instead of later
    for_each(motifs.begin(), motifs.end(), mem_fun_ref(&Motif::clear_sites));

    // The 'sites' vector is where the new sites will be stored as
    // they are obtained. It would be better not to have to store full
    // sites, and only the 'Hit' objects, but converting them to sites
    // later would not work because the actual sequence are needed,
    // and they might not still be available.
    vector<vector<MotifSite> > sites(motifs.size());

    // The variable below keeps track of how much of the total amount
    // of sequence (which might be divided across files) has been
    // processed.
    size_t current_progress = 0;
    const size_t total_file_sizes = (VERBOSE) ?
      get_total_file_sizes(seqfiles) : 0;
    
    // amount by which to increment the file offset for each chunk:
    const size_t file_offset_increment = 
      static_cast<size_t>(buffer_size) - max_width;
    
    // iterate over all the sequence files
    for (size_t i = 0; i < seqfiles.size(); ++i) {
      
      // read in the sequences and their names
      BigFastaFile faa(seqfiles[i]);//, static_cast<size_t>(buffer_size));
      current_progress += faa.get_filesize();
      
      // iterate over the chunks of the file, with 'j' always
      // indicating the offset of the current chunk within the file
      for (size_t j = 0; j < faa.get_filesize(); j += file_offset_increment) {
	
	if (VERBOSE)
	  cerr << file_progress_string(seqfiles[i], 
				       current_progress, total_file_sizes,
				       j, faa.get_filesize());

	// get the sequences in this chunk and their corresponding names
	const vector<string> sequences(faa.get_sequences(j, buffer_size));
	const vector<string> seqnames(faa.get_names(j, buffer_size));
	
	/* 'occ' is a 2D vector that holds the occurrences, outer
	 * vector is for motifs, inner vector is for occurrences of
	 * that motif
	 */
	vector<vector<Hit> > occ(sm.size());
	
	/* This is where we decide which type of search we want, and
	 * call the appropriate function
	 */
	if (!threshold.empty() || bysequence)
	  for (size_t k = 0; k < sequences.size(); ++k) {
	    if (n_top > 0)
	      // one sequence at a time, get top k matches
	      QuerySequenceByCount(sequences[k], sm, smrc, occ, n_top, k);
	    else
	      // one sequence at a time, get everything above a cutoff
	      QuerySequenceByThreshold(sequences[k], sm, smrc, occ, threshold, k);
	  }
	else // do all sequences at once -- we want top k matches in entire file
	  QuerySequenceSetByCount(sequences, sm, smrc, occ, n_top);
	
	// The 'first_seq_offset' is the offset of the first sequence
	// in the current chunk from the start of that entire entire
	// sequence (of which it is a suffix) in the FASTA file.
	const size_t first_seq_offset = faa.get_first_sequence_offset(j);
	
	
	// populate the sites vector from the occ vector
	if (!word_table.c_str())
	  update_motif_sites(occ, sequences, seqnames, full_motif_widths, 
			     first_seq_offset, cores, sites);
	else
	  update_motif_sites(occ, sequences, seqnames, full_motif_widths, 
			     first_seq_offset, sites);
      }
      if (VERBOSE)
	cerr << file_progress_string(seqfiles[i], 1, 1, 1, 1) << endl;
    }
    

    /******************************************************************
     * Now because we went chunk by chunk, we need to figure out which
     * sites to keep, and which ones to throw out. We might have to
     * throw out sites when we wanted the top k sites in a files, but
     * had to get the top k sites in each chunk.  We might also have
     * duplicate sites for motifs that are smaller than the maximum
     * motif size, and therefore can have occurrences in the parts of
     * sequences that are read in consecutive chunks. Those duplicates
     * must also be removed.
     ******************************************************************/
    // First we get rid of duplicates
    remove_duplicate_sites("removing duplicate sites", sites);
    
    // now we make sure we only have the top 'k':
    if (n_top) {
      if (bysequence) // in each sequence
	top_k_per_sequence(sites, n_top, handleties);
      else  // in the whole set
	top_k_overall(sites, n_top, handleties);
    }
    
    /* Since we have our stable set of sites (to be output), we can
     * convert the scores into the desired format (func. depth or
     * p-value):
     */
    if (functional_depth) {
      if (VERBOSE)
	cerr << "converting site scores";
      convert_scores_to_fds(sites, sm);
      if (VERBOSE)
	cerr << "\t" << "done" << endl;
    }
    // !!!TODO: need to make the code able to convert back from scores
    // to p-values; not an easy thing. Would be easier for Staden
    // p-values to start with.
    //     else if (word_table)
    //       convert_scores_to_pvals(sites, wt, cores);
    
    // now we add all the sites to the corresponding motifs
    if (VERBOSE)
      cerr << "adding new sites to motifs";
    for (size_t i = 0; i < sites.size(); ++i)
      for (size_t j = 0; j < sites[i].size(); ++j)
	motifs[i].add_site(sites[i][j]);
    if (VERBOSE)
      cerr << "\t" << "done" << endl;    
    
    // output the motifs (which have updated binding sites)
    ostream* out = (!outfile.empty()) ? new ofstream(outfile.c_str()) : &cout;
    copy(motifs.begin(), motifs.end(), ostream_iterator<Motif>(*out, "\n"));
    if (out != &cout) delete out;
  }
  catch (CREADException &excpt) {
    cerr << "ERROR:\t" << excpt.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
