/*
 * Copyright (C) 2006 Cold Spring Harbor Laboratory
 * Authors: Dustin E. Schones, Andrew D. Smith and Michael Q. Zhang
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

#ifndef WORD_TABLE
#define WORD_TABLE

#include "ScoringMatrix.hpp"
#include "FastaFile.hpp"
#include "Motif.hpp"
#include "cread.hpp"

class WordTable {
private:
  static const size_t maxiter = 50;
  static const size_t buffer_size = 10000;
  static const float tolerance;
  static const float increment;

  size_t nseqs;
  size_t mersize;
  size_t wordsize;
  size_t nwords;
  size_t maxgap;
  size_t gapstart;
  std::vector<float> base_comp;
  std::vector<std::vector<float> > from_left;
  std::vector<std::vector<float> > from_right;
  std::vector<float> mer_freq;
  std::vector<size_t> nposs;
  std::vector<std::vector<size_t> > hits;
  
  size_t count_mers(std::vector<std::string> &seqs, size_t curr);
  size_t count_mers(std::string filename, size_t curr);
  float count_words_above(float **sm, char *prefix, float curr, 
			  float cutoff, size_t depth, size_t maxdepth,
			  size_t gapsize, size_t offset) const;
  float column_info(const float a[]) const;
  std::pair<size_t, size_t> compress(Matrix mat) const;
  float get_cutoff(ScoringMatrix &sm, float pval, size_t gap,
		   size_t offset) const;
  float get_pval(const ScoringMatrix &sm, const float cutoff, 
		 const size_t gap, const size_t offset) const;
  float expected_freq(char *w, size_t wlen, size_t gapsize, size_t offset) const;
  float only_using_freq(char *w, size_t len, size_t gapsize) const;
  static size_t mer2i(const char *s, size_t n);
  size_t mer2i(const char *s, size_t gapsize, size_t n) const;
  static void i2mer(char *s, size_t n, size_t index);
  size_t mer2i_rc(const char *s, size_t gapsize, size_t n) const;
  static size_t onright(size_t word, size_t c, size_t ws);
  static size_t onleft(size_t word, size_t c, size_t ws);
  static size_t find_start(char *s, size_t len, size_t profsize);
  static size_t skip_ahead(char *s, size_t len, size_t profsize, size_t offset);

  static size_t find_next_start(char *s, size_t len, size_t profsize, size_t loc);
  static void base_comp_from_file(std::string filename, float *bc);
  static size_t count_seqs(std::string filename);
  static size_t kmer_counts_from_file(std::string filename,
				      std::vector<size_t> &counts, size_t k);
  
public:
  WordTable() {}
  WordTable(std::vector<std::string> &seqs,
	   size_t wordsize, size_t gapstart, size_t maxgap);
  WordTable(std::string filename,
	   size_t wordsize, size_t gapstart, size_t maxgap);
  static WordTable ReadWordTable(std::string filename);
  
  void WriteWordTable(const std::string &filename) const;
  void WriteWordTableText(std::string filename) const;
  
  size_t get_wordsize() const {return wordsize;}
  size_t get_maxgap() const {return maxgap;}
  size_t get_gapstart() const {return gapstart;}
  float cutoff2pval(const float cutoff, const Matrix& matrix, 
		    const ScoringMatrix& sm, bool both_strands = false) const;
  float pval2cutoff(const float pval, const Matrix& matrix, 
		    const ScoringMatrix& sm, bool both_strands = false) const;
  std::vector<float> extract_base_comp() const {return base_comp;}
  std::string tostring() const;
  std::vector<float> info_profile(Matrix matrix) const;
};

class WordTableException : public CREADException {
public:
  WordTableException(std::string m = std::string()) : CREADException(m) {}
}; 

#endif
