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

#include "cread.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "FastaFile.hpp"
#include "OptionParser.hpp"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::ofstream;
using std::ostream;
using std::transform;
using std::endl;
using std::cout;
using std::cerr;
using std::endl;
using std::copy;
using std::bind2nd;
using std::divides;
using std::max;
using std::multiplies;
using std::fill;
using std::accumulate;
using std::ostream_iterator;

float bias = 0.5;
float epsilon = 0.00001;  

float adenine = 0.0;
float cytosine = 0.0;
float guanine = 0.0;
float thymine = 0.0;

int random_number_generator_seed = std::numeric_limits<int>::max();

static char i2b(int i) {
  static const char *itob = "ACGT";
  return (i < 0 || i > 3) ? -1 : itob[i];
}

static int b2i(char b) {
  static int btoi[20] = {
    //A, b, C, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, T
    0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3
  };
  return (::toupper(b) > 'T' || ::toupper(b) < 'A') ? 
    -1 : btoi[::toupper(b) - 'A'];
}

size_t
count_valid_bases(const std::string& s) {
  return count_if(s.begin(), s.end(), &valid_base);
}

size_t
count_valid_bases(const std::vector<std::string>& s) {
  size_t n_valid = 0;
  for (std::vector<std::string>::const_iterator i = s.begin(); i != s.end(); ++i)
    n_valid += count_valid_bases(*i);
  return n_valid;
}

void
get_base_comp(const std::vector<std::string>& sequences, float *base_comp) {
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

size_t random_nat(size_t upper) {
  return static_cast<size_t>(upper*(rand()/(RAND_MAX + 1.0)));
}

void modify_basecomp(vector<string> &sequences, float *target) {
  float bc[alphabet_size];
  get_base_comp(sequences, bc);
  size_t seq_len = count_valid_bases(sequences);
  
  size_t counts[alphabet_size];
  for (size_t i = 0; i < alphabet_size; ++i)
    counts[i] = static_cast<size_t>(seq_len*bc[i]);
  
  size_t target_counts[alphabet_size];
  for (size_t i = 0; i < alphabet_size; i++) 
    target_counts[i] = static_cast<size_t>(seq_len*target[i]);
  
  int changes[alphabet_size];
  for (size_t i = 0; i < alphabet_size; i++)
    changes[i] = target_counts[i] - counts[i];

  float donate_fraction[alphabet_size];
  for (size_t i = 0; i < alphabet_size; i++)
    donate_fraction[i] = max(0, -changes[i]);
  
  float total = std::accumulate(donate_fraction, 
				donate_fraction + alphabet_size, 0.0);
  transform(donate_fraction, donate_fraction + alphabet_size, 
	    donate_fraction, bind2nd(divides<float>(), total));

  // determine who needs what from whom
  typedef std::pair<size_t, vector<size_t> > mutation_set;
  map<char, map<char, mutation_set> > mutations;
  for (size_t i = 0; i < alphabet_size; i++)
    if (changes[i] > 0) {
      mutations[i2b(i)] = map<char, mutation_set>();
      for (size_t j = 0; j < alphabet_size; j++)
	if (changes[j] < 0) {
	  mutations[i2b(i)][i2b(j)] =
	    mutation_set(static_cast<size_t>(changes[i]*donate_fraction[j]),
			 vector<size_t>());
	}
    }

  vector<set<size_t> > used(alphabet_size);
  for (map<char, map<char, mutation_set> >::iterator i = mutations.begin(); 
       i != mutations.end(); i++) {
    for (map<char, mutation_set>::iterator j = i->second.begin(); 
	 j != i->second.end(); ++j) {
      while (j->second.second.size() < j->second.first) {
	size_t selected;
	do {
	  selected = random_nat(counts[b2i(j->first)]);
	}
	while (used[b2i(j->first)].count(selected) != 0);
	used[b2i(j->first)].insert(selected);
	j->second.second.push_back(selected);
      }
      sort(j->second.second.begin(), j->second.second.end());
    }
  }
  
  // perform mutations
  vector<string> sequences_static = sequences;
  for (map<char, map<char, mutation_set> >::iterator i = mutations.begin(); 
       i != mutations.end(); i++)
    for (map<char, mutation_set>::iterator j = i->second.begin(); 
	 j != i->second.end(); ++j) {
      size_t count = 0;
      vector<size_t>::iterator k = j->second.second.begin();
      for (size_t l = 0; l < sequences_static.size(); ++l)
	for (size_t m = 0; m < sequences_static[l].length(); ++m)
	  if (sequences_static[l][m] == j->first && 
	      k != j->second.second.end()) {
	    if (count == *k) {
	      sequences[l][m] = i->first;
	      ++k;
	    }
	    ++count;
	  }
    }
}

void target_basecomp(float fg_basecomp[alphabet_size], 
			  float target[alphabet_size]) {
  fill(target, target + alphabet_size, 0);
  if (adenine > 0) target[0] = adenine;
  if (cytosine > 0) target[1] = cytosine;
  if (guanine > 0) target[2] = guanine;
  if (thymine > 0) target[3] = thymine;

  float sum = accumulate(target, target + alphabet_size, 0.0);
  if (sum > 1 + epsilon) {
    cerr << "ERROR: target base compositions "
	 << " do not sum to 1\t" << sum << endl;
    exit(EXIT_FAILURE);
  }
  if (adenine <= 0) 
    target[0] = fg_basecomp[0] * (1 - sum);
  if (cytosine <= 0) 
    target[1] = fg_basecomp[1] * (1 - sum);
  if (guanine <= 0) 
    target[2] = fg_basecomp[2] * (1 - sum);
  if (thymine <= 0) 
    target[3] = fg_basecomp[3] * (1 - sum);
  sum = accumulate(target, target + alphabet_size, 0.0);
  if (sum > 1 + epsilon) {
    cerr << "ERROR: target base compositions "
	 << "do not sum to 1\t" << sum << endl;
    exit(EXIT_FAILURE);
  }
}


int main(int argc, const char **argv) {

  try {

    string fgfile;        // Name of input file containing foreground sequences
    string bgfile;        // Name of input file containing foreground sequences
    string fg_outfile;                 // Name of output file for foreground sequences
    string bg_outfile;                 // Name of output file for foreground sequences

    /*********************** COMMAND LINE OPTIONS ************************************/
    OptionParser opt_parse(strip_path(argv[0]), " ",
                           " ");   
    opt_parse.add_opt("bias", 'i', "fraction representing contribution of"
                    " foreground to target base composition", false, bias);
    opt_parse.add_opt(" ", 'A', "target adenine composition", false, adenine);
    opt_parse.add_opt(" ", 'C', "target cytosine composition", false, cytosine);
    opt_parse.add_opt(" ", 'G', "target guanine composition", false, guanine);
    opt_parse.add_opt(" ", 'T',"target thymine composition", false, thymine);
    opt_parse.add_opt("background", 'b', "name of background file", 
                       false, bgfile);
    opt_parse.add_opt("fg-output", '\0', "foreground output file", 
                       false, fg_outfile);
    opt_parse.add_opt("bg-output", '\0', "background output file", 
                       false, bg_outfile);
    opt_parse.add_opt("rng-seed", '\0', "seed for random number generator", 
                       false, random_number_generator_seed);

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
    const vector<string> input_files(leftover_args);

  /**********************************************************************/
    // Seed the random number generator
    srand(random_number_generator_seed);
    
    fgfile = input_files[0];
    FastaFile ff(fgfile); 
    vector<string> fgnames = ff.get_names();
    vector<string> foreground = ff.get_sequences(); 

    float fg_basecomp[alphabet_size];
    get_base_comp(foreground, fg_basecomp);
  
    float target[alphabet_size];
    if (!bgfile.empty()) {
      FastaFile bb(bgfile);
      vector<string> fgnames = bb.get_names();
      vector<string> background = bb.get_sequences();
      
      float bg_basecomp[alphabet_size];
      get_base_comp(background, bg_basecomp); 

      // Calculate the target base composition

      for (size_t i = 0; i < alphabet_size; i++)
	target[i] = bias * fg_basecomp[i] + (1 - bias) * bg_basecomp[i];
      
      modify_basecomp(background, target);

      ostream* seqout = (!bg_outfile.empty()) ? new ofstream(bg_outfile.c_str()) : &cout;
      for (size_t i = 0; i < background.size(); ++i)
	*seqout << ">" << fgnames[i] << endl << background[i] << endl;
      if (seqout != &cout) delete seqout;
      
    }
    else target_basecomp(fg_basecomp, target);
    
    modify_basecomp(foreground, target);
    
    ostream* seqout = (!fg_outfile.empty()) ? new ofstream(fg_outfile.c_str()) : &cout;
    for (size_t i = 0; i < foreground.size(); ++i)
      *seqout << ">" << fgnames[i] << endl << foreground[i] << endl;
    if (seqout != &cout) delete seqout;
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
