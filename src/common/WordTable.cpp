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

#include "WordTable.hpp"
#include "smithlab_utils.hpp"

using std::string;
using std::ostream;
using std::ofstream;
using std::ostream_iterator;
using std::vector;
using std::pair;
using std::fill_n;
using std::copy;
using std::map;
using std::pow;
using std::cout;
using std::endl;
using std::count_if;
using std::ostringstream;
using std::min;
using std::cerr;

const float WordTable::tolerance = 0.0001;
const float WordTable::increment = 0.05;

void
inline get_base_comp(const std::vector<std::string>& sequences, float *base_comp) {
  std::fill(base_comp, base_comp + smithlab::alphabet_size, 0.0);
  float total = 0;
  for (std::vector<std::string>::const_iterator i = sequences.begin();
       i != sequences.end(); ++i)
    for (std::string::const_iterator j = i->begin(); j != i->end(); ++j)
      if (valid_base(*j)) {
        base_comp[base2int(*j)]++;
        total++;
      }
  transform(base_comp, base_comp + smithlab::alphabet_size, base_comp,
            std::bind(std::divides<float>(), std::placeholders::_1, total));
}

WordTable::WordTable(std::vector<string> &seqs, size_t wordsize, 
		     size_t gapstart, size_t maxgap) :
  nseqs(seqs.size()), wordsize(wordsize), 
  maxgap(maxgap), gapstart(gapstart) {
  
  base_comp = vector<float>(alphabet_size);
  get_base_comp(seqs, &base_comp.front());
  
  mersize = std::min(gapstart, wordsize - gapstart);
  vector<size_t> counts;
  float total = kmer_counts(seqs, counts, mersize);
  
  std::transform(counts.begin(), counts.end(), back_inserter(mer_freq),
		 bind2nd(std::divides<float>(), total));
  nwords = static_cast<size_t>(pow(static_cast<float>(alphabet_size), 
				   static_cast<int>(wordsize)));
  hits = vector<vector<size_t> >(nwords, vector<size_t>(maxgap + 1));
  for (size_t i = 0; i < maxgap + 1; ++i)
    nposs.push_back(count_mers(seqs, i));
}

size_t 
WordTable::count_seqs(string filename) {
  std::ifstream fin(filename.c_str());
  size_t count = 0;
  char c;
  while (fin.get(c)) if (c == '>') ++count;
  fin.close();
  return count;
}

void
WordTable::base_comp_from_file(string filename, float *bc) {
  std::ifstream fin(filename.c_str());
  bool reading_sequence = true;
  char buffer[buffer_size + 1];
  int total = 0;
  int tempbc[alphabet_size];
  std::fill(tempbc, tempbc + alphabet_size, 0);
  while (!fin.eof()) {
    fin.get(buffer, buffer_size + 1, '\0');
    for (int j = 0; j < fin.gcount(); ++j)
      if (buffer[j] == '>') reading_sequence = false;
      else if (!reading_sequence && buffer[j] == '\n') 
	reading_sequence = true;
      else if (reading_sequence && buffer[j] != '\n' && 
	       valid_base(buffer[j])) {
	tempbc[base2int(buffer[j])]++;
	++total;
      }
  }
  std::transform(tempbc, tempbc + alphabet_size, bc,
		 bind2nd(std::divides<double>(), total));
  fin.close();
}

size_t 
WordTable::kmer_counts_from_file(string filename, 
				vector<size_t> &counts, size_t k) {
  counts.clear();
  size_t nwords = static_cast<size_t>(pow(static_cast<float>(alphabet_size), 
					  static_cast<int>(k)));
  counts.resize(nwords, 0);
  size_t total = 0;
  char buffer[buffer_size + 1];
  std::ifstream fin(filename.c_str());
  bool good = true;
  bool reading_first = true;
  char seq[buffer_size + k + 1];
  //std::fill_n(seq, buffer_size + k, '\0');
  while (!fin.eof()) {
    fin.get(buffer, buffer_size + 1, '\0');
    size_t prefix_len = 0, len = fin.gcount();
    if (!reading_first && len > 0) {
      prefix_len = k - 1;
      const size_t ll = strlen(seq) - prefix_len;
      for (size_t j = 0; j < prefix_len; ++j) seq[j] = seq[ll + j];
      seq[prefix_len] = '\0';
    }
    else reading_first = false;
    std::fill_n(seq + prefix_len, len, '\0');
    for (size_t j = 0; j < len; ++j)
      if (buffer[j] == '>') good = false;
      else if (!good && buffer[j] == '\n') {
	seq[prefix_len++] = 'N';
	good = true;
      }
      else if (good && buffer[j] != '\n') seq[prefix_len++] = buffer[j];
    const size_t frames = (strlen(seq) >= k) ? strlen(seq) - k + 1 : 0;
    for (size_t j = 0; j < frames; ++j)
      if (std::count_if(seq + j, seq + j + k, &valid_base) == 
	  static_cast<int>(k)) {
	counts[mer2index(seq + j, k)]++;
	++total;
      }
  }
  fin.close();
  return total;
}

WordTable::WordTable(string filename, size_t wordsize, 
		   size_t gapstart, size_t maxgap) :
  wordsize(wordsize), maxgap(maxgap), gapstart(gapstart) {
  
  nseqs = count_seqs(filename);
  base_comp = vector<float>(alphabet_size);
  base_comp_from_file(filename, &base_comp.front());
  
  mersize = std::min(gapstart, wordsize - gapstart);
  vector<size_t> counts;
  const float total = kmer_counts_from_file(filename, counts, mersize);
  std::transform(counts.begin(), counts.end(), back_inserter(mer_freq),
		 bind2nd(std::divides<float>(), total));
  nwords = static_cast<size_t>(pow(static_cast<float>(alphabet_size), 
				   static_cast<int>(wordsize)));
  hits = vector<vector<size_t> >(nwords, vector<size_t>(maxgap + 1));
  for (size_t i = 0; i < maxgap + 1; ++i)
    nposs.push_back(count_mers(filename, i));
}

size_t
WordTable::find_next_start(char *s, size_t len, size_t profsize, size_t loc){
  size_t good_count = 0;
  while (good_count < profsize && loc <= len){
    if (valid_base(s[loc++]))
      ++good_count;
    else
      good_count = 0;
  }
  return loc - profsize;
}

size_t
WordTable::count_mers(vector<string> &seqs, size_t curr) {
  const size_t profsize = wordsize + curr;
  size_t possible = 0;
  for (size_t i = 0; i < seqs.size(); ++i) {
    char seq[seqs[i].length() + 1];
    seq[seqs[i].length()] = '\0';
    copy(seqs[i].begin(), seqs[i].end(), seq);

    const size_t frames = (seqs[i].length() >= profsize) ? 
      seqs[i].length() - profsize + 1 : 0;

    for (size_t j = find_next_start(seq, seqs[i].length(), profsize, 0); 
	 j < frames; ) {
      if (valid_base(seq[j + profsize - 1])) {
	hits[mer2i(seq + j, curr, profsize)][curr]++;
	++possible;
	++j;
      }
      else j = find_next_start(seq, seqs[i].length(), profsize, j);
    }
  }
  return possible;
}

size_t
WordTable::count_mers(string filename, size_t curr) {
  const size_t profsize = wordsize + curr;
  size_t possible = 0;
  char buffer[buffer_size + 1];
  std::ifstream fin(filename.c_str());
  bool good = true;
  bool reading_first = true;
  char seq[buffer_size + profsize + 1];
  std::fill_n(seq, buffer_size + profsize, '\0');
  while (!fin.eof()) {
    fin.get(buffer, buffer_size + 1, '\0');
    size_t prefix_len = 0, len = fin.gcount();
    if (!reading_first && len > 0) {
      prefix_len = profsize - 1;
      size_t ll = strlen(seq) - prefix_len;
      for (size_t j = 0; j < prefix_len; ++j)
	seq[j] = seq[ll + j];
      seq[prefix_len] = '\0';
    }
    else reading_first = false;
    std::fill_n(seq + prefix_len, len, '\0');
    for (size_t j = 0; j < len; ++j)
      if (buffer[j] == '>') good = false;
      else if (!good && buffer[j] == '\n') {
	seq[prefix_len++] = 'N';
	good = true;
      }
      else if (good && buffer[j] != '\n') seq[prefix_len++] = buffer[j];
    
    const size_t frames = (strlen(seq) >= profsize) ? 
      strlen(seq) - profsize + 1 : 0;
    
    for (size_t j = find_next_start(seq, strlen(seq), profsize, 0);
 	 j < frames; ) {
      if (valid_base(seq[j + profsize - 1])) {
	hits[mer2i(seq + j, curr, profsize)][curr]++;
	++possible;
	++j;
      }
      else j = find_next_start(seq, strlen(seq), profsize, j);
    }
  }
  fin.close();
  return possible;
}

string
WordTable::tostring() const {
  ostringstream ss;
  ss << wordsize << endl << gapstart << endl
     << maxgap << endl << nseqs << endl;
  copy(base_comp.begin(), base_comp.end(),
       ostream_iterator<float>(ss, "\n"));
  ss << mersize << endl;
  copy(mer_freq.begin(), mer_freq.end(),
       ostream_iterator<float>(ss, "\n"));
  copy(nposs.begin(), nposs.end(),
       ostream_iterator<size_t>(ss, "\n"));
  size_t count = 0;
  for (vector<vector<size_t> >::const_iterator i = hits.begin(); 
       i != hits.end(); ++i) {
    char mymer[100];
    i2mer(mymer, wordsize, count);
    mymer[wordsize] = '\0';
    ss << mymer << "\t";
    copy(i->begin(), i->end(), ostream_iterator<size_t>(ss, "\t"));
    ss << endl;
    ++count;
  }
  return ss.str();
}

void
WordTable::WriteWordTableText(string filename = "") const {
  ofstream out(filename.c_str(), std::ios::out );
  out << wordsize << " " << sizeof(wordsize) << endl 
      << gapstart << " " << sizeof(gapstart) << endl
      << maxgap << " " << sizeof(maxgap) <<  endl 
      << nseqs << " " << sizeof(nseqs) << endl;
  copy(base_comp.begin(), base_comp.end(),
       ostream_iterator<float>(out, "\n"));
  out << alphabet_size*sizeof(float) << endl;
  out << mersize << " " << sizeof(mersize) << endl;
  copy(mer_freq.begin(), mer_freq.end(),
       ostream_iterator<float>(out, "\n"));
  out << mer_freq.size()*sizeof(float) << endl;
  copy(nposs.begin(), nposs.end(),
       ostream_iterator<size_t>(out, "\t")); 
  out << nposs.size()*sizeof(size_t) << endl;
  vector<vector<size_t> >::const_iterator i = hits.begin();
  for (; i != hits.end(); ++i) {
    copy(i->begin(), i->end(), ostream_iterator<size_t>(out, "\t"));
    out << i->size()*sizeof(size_t);
    out << endl;
  }
  out.close();
}


void
WordTable::WriteWordTable(const string &filename) const {
  ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
  const int int_wordsize = static_cast<int>(wordsize);
  out.write((char*)&int_wordsize, sizeof(int_wordsize));
  const int int_gapstart = static_cast<int>(gapstart);
  out.write((char*)&int_gapstart, sizeof(int_gapstart));
  const int int_maxgap = static_cast<int>(maxgap);
  out.write((char*)&int_maxgap, sizeof(int_maxgap));
  const int int_nseqs = static_cast<int>(nseqs);
  out.write((char*)&int_nseqs, sizeof(int_nseqs));
  out.write((char*)&base_comp.front(), alphabet_size*sizeof(float));
  const int int_mersize = static_cast<int>(mersize);
  out.write((char*)&int_mersize, sizeof(int_mersize));
  out.write((char*)&mer_freq.front(), mer_freq.size()*sizeof(float));

  for (vector<size_t>::const_iterator i = nposs.begin(); 
       i != nposs.end(); ++i) {
    const int int_nposs = static_cast<int>(*i);
    out.write((char*)&int_nposs, sizeof(int_nposs));
  }
  for (vector<vector<size_t> >::const_iterator i = hits.begin(); 
       i != hits.end(); ++i) {
    for (vector<size_t>::const_iterator j = i->begin(); j != i->end(); ++j) {
      const int int_hits = static_cast<int>(*j);
      out.write((char*)&int_hits, sizeof(int_hits));
    }
  }
  out.close();
}


WordTable
WordTable::ReadWordTable(string filename) {
  std::ifstream fin(filename.c_str(), std::ios::binary);
  if (!fin) {
    string message("ERROR in WordTable:\n"
		   "cannot open table file: " + filename);
    throw WordTableException(message);
  }
  
  WordTable wt;
  // First read the word size
  int int_wordsize;
  fin.read((char*)&int_wordsize, sizeof(int_wordsize));
  wt.wordsize = static_cast<size_t>(int_wordsize);
  
  // Now get the number of words
  wt.nwords = static_cast<size_t>(pow(static_cast<float>(alphabet_size), 
				      int_wordsize));
  
  // Read start of the gap
  int int_gapstart;
  fin.read((char*)&int_gapstart, sizeof(int_gapstart));
  wt.gapstart = static_cast<size_t>(int_gapstart);
  
  // Read the max gaps size
  int int_maxgap;
  fin.read((char*)&int_maxgap, sizeof(int_maxgap));
  wt.maxgap = static_cast<size_t>(int_maxgap);
  
  // Read the number of seqs
  int int_nseqs;
  fin.read((char*)&int_nseqs, sizeof(int_nseqs));
  wt.nseqs = static_cast<size_t>(int_nseqs);
  
  // Next read the base composition
  wt.base_comp = vector<float>(alphabet_size);
  fin.read((char*)&wt.base_comp[0], alphabet_size*sizeof(float));
  
  // Now read the short mer composition
  int int_mersize;
  fin.read((char*)&int_mersize, sizeof(int_mersize));
  wt.mersize = static_cast<size_t>(int_mersize);
  
  // calculate the number of mers
  const size_t n_mers = 
    static_cast<size_t>(pow(static_cast<float>(alphabet_size), int_mersize));
  
  wt.mer_freq = vector<float>(n_mers);
  fin.read((char*)&wt.mer_freq.front(), n_mers*sizeof(float));
  
  size_t nshorter = n_mers/alphabet_size;
  wt.from_left =
    vector<vector<float> >(nshorter, vector<float>(alphabet_size, 0.0));
  wt.from_right = 
    vector<vector<float> >(nshorter, vector<float>(alphabet_size, 0.0));
  
  for (size_t i = 0; i < nshorter; ++i) {
    float total_left = 0, total_right = 0;
    for (size_t j = 0; j < alphabet_size; ++j) {
      wt.from_left[i][j] = wt.mer_freq[onright(i, j, wt.mersize)];
      total_left += wt.from_left[i][j];
      wt.from_right[i][j] = wt.mer_freq[onleft(i, j, wt.mersize)];
      total_right += wt.from_right[i][j];
    }
    for (size_t j = 0; j < alphabet_size; ++j) {
      wt.from_left[i][j] /= total_left;
      wt.from_right[i][j] /= total_right;
    }
  }
  
  // Next read the total number of possible hits for each gap size
  vector<int> int_nposs(wt.maxgap + 1);
  fin.read((char*)&int_nposs.front(), 
	   (wt.maxgap + 1)*sizeof(int_nposs.front()));
  // and convert it to size_t
  wt.nposs = vector<size_t>(wt.maxgap + 1);
  vector<size_t>::iterator j = wt.nposs.begin();
  for (vector<int>::iterator i = int_nposs.begin();
       i != int_nposs.end(); ++i, ++j)
    *j = static_cast<size_t>(*i);
  
  // Now read the words' frequencies
  wt.hits = vector<vector<size_t> >(wt.nwords, vector<size_t>(wt.maxgap + 1));
  vector<int> int_hits(wt.maxgap + 1);
  for (size_t i = 0; i < wt.nwords; ++i) {
    fin.read((char*)&int_hits.front(), 
	     (wt.maxgap + 1)*sizeof(int_hits.front()));
    for (size_t j = 0; j <= wt.maxgap; ++j)
      wt.hits[i][j] = static_cast<size_t>(int_hits[j]);
  }
  return wt;
}

float
WordTable::count_words_above(float **sm, char *prefix, float curr, 
			    float cutoff, size_t depth, size_t maxdepth,
			    size_t gapsize, size_t offset) const {
  if (depth == maxdepth)
    return expected_freq(prefix, maxdepth, gapsize, offset);
  //return only_using_freq(prefix, maxdepth, gapsize);
  float total = 0;
  for (size_t i = 0; i < alphabet_size; ++i)
    if (curr - sm[depth][i] >= cutoff) {
      prefix[depth] = int2base(i);
      total += count_words_above(sm, prefix, curr - sm[depth][i], cutoff,
				 depth + 1, maxdepth, gapsize, offset);
    }
  return total;
}


float
WordTable::expected_freq(char *w, size_t wlen, 
			size_t gapsize, size_t offset) const {
  const size_t part2end = offset + gapstart + gapsize;
  float freq = 1.0;
  for (int j = offset - 1; j >= 0; --j) // Freq for part 1
    freq *= from_right[mer2i(w + j + 1, 
			     mersize - 1)][base2int(w[j])];
  for (size_t j = offset + gapstart; j < part2end; ++j) // Freq for part 2
    freq *= from_left[mer2i(w + j - mersize + 1, 
			    mersize - 1)][base2int(w[j])];
  for (size_t j = offset + gapsize + wordsize; j < wlen; ++j) // Freq for part 3
    freq *= from_left[mer2i(w + j - mersize + 1, 
			    mersize - 1)][base2int(w[j])];
  return hits[mer2i(w + offset, gapsize, 
		    wordsize + gapsize)][gapsize]*freq;
}


float
WordTable::only_using_freq(char *w, size_t len, size_t gapsize) const {
  float freq = mer_freq[mer2i(w, mersize)];
  for (size_t j = mersize; j < len; ++j)
    freq *= from_left[mer2i(w + j - mersize + 1,
			    mersize - 1)][base2int(w[j])];
  return nposs[gapsize]*freq;
}

float
WordTable::get_pval(const ScoringMatrix &sm, const float cutoff, 
		   const size_t gapsize, const size_t offset) const {
  float max_score = 0;
  const size_t matsize = std::max(wordsize, sm.get_width());
  float **score_matrix = new float *[matsize];
  
  for (size_t i = 0; i < matsize; ++i) {
    float tempmax = 0;
    if (i < sm.get_width())
      for (size_t j = 0; j < alphabet_size; j++)
	tempmax = std::max(tempmax, sm[i][j]);
    max_score += tempmax;
    score_matrix[i] = new float[alphabet_size];
    if (i < sm.get_width())
      for (size_t j = 0; j < alphabet_size; ++j)
	score_matrix[i][j] = tempmax - sm[i][j];
    else fill_n(score_matrix[i], alphabet_size, 0);
  }
  vector<char> prefix(matsize + 1, '\0');
  const float total = count_words_above(score_matrix, &prefix.front(), 
					max_score, min(max_score, cutoff), 0, 
					matsize, gapsize, offset);
  for (size_t i = 0; i < matsize; ++i)
    delete[] score_matrix[i];
  delete[] score_matrix;
  
  return total/nposs[gapsize];
}



float
WordTable::cutoff2pval(const float cutoff, const Matrix& mat, 
		      const ScoringMatrix& sm, bool both_strands) const {
  pair<size_t, size_t> gap_off = compress(mat);
  const size_t matsize = std::max(wordsize, sm.get_width());
  float pval = get_pval(sm, cutoff, gap_off.first, gap_off.second);
  if (both_strands) {
    ScoringMatrix smrc(sm.revcomp());
    pval += get_pval(smrc, cutoff, gap_off.first, 
		     matsize - (gap_off.second + gap_off.first + wordsize));
    pval /= 2;
  }
  return pval;
}



float
WordTable::pval2cutoff(const float pval, const Matrix& mat, 
		      const ScoringMatrix& sm, bool both_strands) const {
  const pair<size_t, size_t> gap_off(compress(mat));
  if (get_pval(sm, sm.functional_depth_to_score(1.0), 
	       gap_off.first, gap_off.second) > pval)
    return std::numeric_limits<float>::max();
  
  const ScoringMatrix smrc(sm.revcomp());
  float curr = 1.0, top = 1.0, bottom = 1.0;
  const size_t matsize = std::max(wordsize, sm.get_width());
  do {
    top = bottom;
    bottom -= increment;
    curr = get_pval(sm, sm.functional_depth_to_score(bottom), 
		    gap_off.first, gap_off.second);
    if (both_strands) {
      curr += get_pval(smrc, sm.functional_depth_to_score(bottom), 
		       gap_off.first, 
		       matsize - (gap_off.second + gap_off.first + wordsize));
      curr /= 2;
    }
  } while (curr < pval);
  size_t iterations = 0;
  while (std::abs(curr - pval)/std::max(curr, pval) > tolerance && 
	 iterations++ < maxiter) {
    const float mid = (bottom + top)/2;
    curr = get_pval(sm, sm.functional_depth_to_score(mid), 
		    gap_off.first, gap_off.second);
    if (both_strands) {
      curr += get_pval(smrc, sm.functional_depth_to_score(mid), 
		       gap_off.first,
		       matsize - (gap_off.second + gap_off.first + wordsize));
      curr /= 2;
    }
    if (curr > pval) bottom = mid;
    else top = mid;
  }
  return sm.functional_depth_to_score(bottom);
}



float 
WordTable::column_info(const float a[]) const {
  const float sum = std::accumulate(a, a + alphabet_size, 0.0);
  float info = 0;
  for (size_t i = 0; i < alphabet_size; ++i)
    if (a[i] > 0) {
      float freq = a[i]/sum;
      info += freq*(cread::log2(freq) - cread::log2(base_comp[i]));
    }
  return info;
}

vector<float> 
WordTable::info_profile(Matrix matrix) const {
  vector<float> profile;
  for (size_t i = 0; i < matrix.get_width(); ++i)
    profile.push_back(column_info(matrix[i]));
  return profile;
}


pair<size_t, size_t>
WordTable::compress(Matrix mat) const {
  const vector<float> prof = info_profile(mat.freqmat());
  float max_info = 0;
  size_t offset = 0, gap = 0;
  const size_t lim = (wordsize + 1 > mat.get_width()) ? 0 : 
    mat.get_width() - wordsize + 1;
  for (size_t i = 0; i < lim; ++i)
    for (size_t j = 0; i + wordsize + j < mat.get_width() && j <= maxgap; ++j) {
      float temp_info = 0;
      for (size_t k = 0; k < gapstart; ++k) temp_info += prof[i + k];
      for (size_t k = 0; k < wordsize - gapstart; ++k)
	temp_info += prof[i + gapstart + j + k];
      if (temp_info > max_info) {
	max_info = temp_info;
	offset = i;
	gap = j;
      }
    }
  return pair<size_t, size_t>(gap, offset);
}

size_t
WordTable::mer2i(const char *s, size_t n) {
  size_t multiplier = 1, index = 0;
  do {
    --n;
    index += base2int(s[n])*multiplier;
    multiplier *= alphabet_size;
  } while (n > 0);
  return index;
}

size_t
WordTable::mer2i(const char *s, size_t gapsize, size_t n) const {
  size_t multiplier = 1, index = 0;
  const size_t gapend = wordsize - gapstart + gapsize;
  do {
    --n;
    index += base2int(s[n])*multiplier;
    multiplier *= alphabet_size;
  } while (n > gapend);
  n -= gapsize;
  do {
    --n;
    index += base2int(s[n])*multiplier;
    multiplier *= alphabet_size;
  } while (n > 0);
  return index;
}

size_t
WordTable::mer2i_rc(const char *s, size_t gapsize, size_t n) const {
  size_t multiplier = 1, index = 0, k = 0;
  do {
    index += (alphabet_size - 1 - base2int(s[k]))*multiplier;
    multiplier *= alphabet_size;
  } while (++k < gapstart);
  k += gapsize;
  do {
    index += (alphabet_size - 1 - base2int(s[k]))*multiplier;
    multiplier *= alphabet_size;
  } while (++k < n);
  return index;
}

void
WordTable::i2mer(char *s, size_t n, size_t index) {
  do {
    --n;
    s[n] = int2base(index % alphabet_size);
    index /= alphabet_size;
  } while (n > 0);
}

size_t 
WordTable::onright(size_t word, size_t c, size_t ws) {
  return word*alphabet_size + c;
}

size_t 
WordTable::onleft(size_t word, size_t c, size_t ws) {
  const size_t p = 
    static_cast<size_t>(pow(static_cast<float>(alphabet_size), 
			    static_cast<int>(ws - 1)));
  return p*c + word;
}
