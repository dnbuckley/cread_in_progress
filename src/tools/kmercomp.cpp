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

// TODO: Make this program report kmer composition by sequence (like
// TODO: -q in storm) - ADS
// TODO: Report kmers in motif format (until a word format is available)

#include "cread.hpp"
#include "Word.hpp"
//#include "Alphabet.hpp"
#include "FastaFile.hpp"
#include "SuffixTree.hpp"
#include "WordMatcher.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"

using std::ofstream;
using std::ostringstream;
using std::ostream_iterator;
using std::numeric_limits;
using std::string;
using std::vector;
using std::pair;
using std::ptr_fun;
using std::copy;
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;

const char *kmer_count_label = "COUNT";
const char *kmer_freq_label = "FREQUENCY";
const char *kmer_freq_lod_label = "FREQUENCY_LOD";
const char *kmer_freq_rel_lod_label = "RELATIVE_FREQUENCY_LOD";

float
expected_freq(const string &s, const vector<float>& base_comp) {
  float freq = 1;
  for (string::const_iterator i = s.begin(); i != s.end(); ++i)
    freq *= base_comp[base2int(*i)];
  return freq;
}

Word 
make_kmer(string w, string a, size_t count, float total,
	  const vector<float>& base_comp) {
  Word word(w, a);
  word.set_attribute(kmer_count_label, count);
  word.set_attribute(kmer_freq_label, count/total);
  word.set_attribute(kmer_freq_lod_label, cread::log2(count/total) - 
		     cread::log2(expected_freq(w, base_comp)));
  return word;
}

struct Greater2ndKey {
  template <class T, class U> 
  bool operator()(const pair<T, U> &a, const pair<T, U> &b) const {
    return a.second > b.second;
  }
};

void
i2mer(char *s, size_t n, size_t index) {
  do {
    --n;
    s[n] = int2base(index % alphabet_size);
    index /= alphabet_size;
  } while (n > 0);
}

float
max_valid_matches(const vector<string>& sequences, 
		  const size_t width) {
  float total = 0;
  for (size_t i = 0; i < sequences.size(); ++i)
    total += sequences[i].length() - width + 1;
    //     for (size_t j = 0; j < sequences[i].length(); ++j) {
    //       total += this substring if it is valid;
    //     }
  return total;
}

int main(int argc, const char **argv) {

  // Parameter variables
  int kmer = 1;
  int ntop = numeric_limits<int>::max();
  string outfile;
  string seqfile;
  string kmersfile;

  string prefix = "KMER";
  try {
    /********************** COMMAND LINE OPTIONS **************************/
    OptionParser opt_parse(strip_path(argv[0]), " ",
                           " ");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                       false, outfile);
    opt_parse.add_opt("prefix", 'p',  "accession prefix for generated k-mers",
                       false, prefix);
    opt_parse.add_opt(" ", 'k', "\"k\" in k-mer", false, kmer);
    opt_parse.add_opt(" ", 'n', "number of top k-mers to output", false, ntop);
    opt_parse.add_opt("words", 'w', "file of k-mers to query", false, kmersfile);

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
    /**********************************************************************/
    // read sequences
    seqfile = input_filenames[0];
    const vector<string> sequences(FastaFile(seqfile).get_sequences());
    vector<float> base_comp;
    get_base_comp(sequences, base_comp);

    vector<Word> words;
    if (!kmersfile.empty()) {
      words = Word::ReadWordVector(kmersfile.c_str());
      for (size_t i = 0; i < words.size(); ++i) {
	const float total = max_valid_matches(sequences, words[i].get_width());
	vector<string> matches;
	get_matches(words[i], sequences, words[i].get_width(), matches);
	words[i].set_attribute(kmer_count_label, matches.size());
	words[i].set_attribute(kmer_freq_label, matches.size()/total);
	words[i].set_attribute(kmer_freq_lod_label, 
			       cread::log2(matches.size()/total) - 
			       cread::log2(expected_freq(words[i].get_word(), 
							 base_comp)));
      }
      ostream* output = (!outfile.empty()) ? new ofstream(outfile.c_str()) : &cout;
      copy(words.begin(), words.end(), ostream_iterator<Word>(*output, "\n"));
      if (output != &cout) delete output;
    }
    else {
      vector<size_t> counts;
      const float total = kmer_counts(sequences, counts, kmer);

      ostream* output = (!outfile.empty()) ? new ofstream(outfile.c_str()) : &cout;
      for (size_t i = 0; i < counts.size(); ++i) {
	char char_kmer[kmer + 1];
	char_kmer[kmer] = '\0';
	i2mer(char_kmer, kmer, i);
	const string the_kmer(char_kmer);
 	*output << make_kmer(the_kmer, the_kmer,
			     counts[i], total, base_comp) << endl;
      }
      if (output != &cout) delete output;
    }
  }
  catch (CREADException &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
