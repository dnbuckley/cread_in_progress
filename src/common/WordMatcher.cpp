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

#include "WordMatcher.hpp"

using std::vector;
using std::string;
using std::sort;
using std::min;

struct WordMatch {
  size_t mismatches;
  size_t seqid;
  size_t offset;
  WordMatch(const size_t m, const size_t s, const size_t o) : 
    mismatches(m), seqid(s), offset(o) {}
  bool operator<(const WordMatch& rhs) const {
    return mismatches < rhs.mismatches;
  }
};

static void
get_word_max_mismatch_sequence_set(const string& word,
				   const vector<string>& seqs,
				   const size_t max_mismatches,
				   vector<WordMatch>& matches) {
  for (size_t i = 0; i < seqs.size(); ++i) {
    const size_t n_substrings = seqs[i].length() - word.length() + 1;
    for (size_t j = 0; j < n_substrings; ++j) {
      size_t mismatches = 0;
      for (size_t k = 0; k < word.length() && mismatches <= max_mismatches; ++k)
	mismatches += (seqs[i][j + k] != word[k]);
      if (mismatches <= max_mismatches)
	matches.push_back(WordMatch(mismatches, i, j));
    }
  }
}

void
get_matches(const Word& word, const vector<string>& seqs,
	    const size_t min_matching_positions,
	    vector<string>& matches, size_t max_matches) {
  vector<WordMatch> hits;
  const size_t max_mismatches = word.get_width() - min_matching_positions;
  get_word_max_mismatch_sequence_set(word.get_word(), seqs, max_mismatches, hits);
  sort(hits.begin(), hits.end());
  for (size_t i = 0; i < min(hits.size(), max_matches); ++i)
    matches.push_back(seqs[hits[i].seqid].substr(hits[i].offset,
						 word.get_width()));
}

static bool
sequence_has_site_max_mismatch(const string& word,
			       const string& seq,
			       const size_t max_mismatches) {
  const size_t limit = seq.length() - word.length() + 1;
  for (size_t i = 0; i < limit; ++i) {
    size_t mismatches = 0;
    for (size_t j = 0; j < word.length() && mismatches <= max_mismatches; ++j)
      mismatches += (seq[i + j] != word[j]);
    if (mismatches <= max_mismatches)
      return true;
  }
  return false;
}

size_t
count_sequences_with_match(const Word& word, 
			   const std::vector<std::string>& seqs,
			   const size_t min_matching_positions) {
  const string word_repr = word.get_word();
  const size_t max_mismatches = word.get_width() - min_matching_positions;
  size_t sequences_with_match = 0;
  for (size_t i = 0; i < seqs.size(); ++i)
    sequences_with_match += sequence_has_site_max_mismatch(word_repr,
							   seqs[i],
							   max_mismatches);
  return sequences_with_match;
}
