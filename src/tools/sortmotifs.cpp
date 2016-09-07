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
#include "Pattern.hpp"
#include "PatternFactory.hpp"
#include "Motif.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::ostream_iterator;
using std::ofstream;
using std::ostream;
using std::min;
using std::copy;
using std::cerr;
using std::endl;
using std::cout;


string key;
const char *rank_suffix = "_RANK";

void AssignRanks(vector<Pattern*>& patterns) {
  string key_rank = string(key) + string(rank_suffix);
  for (size_t i = 0, rank = 1; i < patterns.size(); ++i) {
    patterns[i]->set_attribute(key_rank, rank);
    if (i == patterns.size() - 1 ||
	patterns[i]->get_attribute(key) != patterns[i + 1]->get_attribute(key))
      rank = i + 2;
  }
}

int main(int argc, const char **argv) {
  
  // COMMAND LINE PARAMETERS
  string outfile;
  string patterns_file;     
  string cutoff;
  static int max_patterns_to_output_param = -1;
  bool reverse_order = false;
  bool numeric_sort_order = false;
  bool processing_motifs = false;
  bool assign_rank = false;

  try {
    /***************** COMMAND LINE OPTIONS *******************/
    OptionParser opt_parse(strip_path(argv[0]), "sorts patterns in CREAD"
                           " pattern by a key value", "  ");
    opt_parse.add_opt("reverse", 'r', "sort in reverse order", 
                         false, reverse_order);
    opt_parse.add_opt("rank", 'a', "set a rank attribute",
                       false, assign_rank);
    opt_parse.add_opt("numeric", 'n',"sort as though key has numeric value",
                         false, numeric_sort_order);
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                         false, outfile);
    opt_parse.add_opt("key", 'k', "name of attribute on which to sort",
                         false, key);
    opt_parse.add_opt("cutoff", 'c', "cutoff value for patterns to return", 
                       false, cutoff);
    opt_parse.add_opt("ntop", 'm', "only print this many of top outputs", 
                       false, max_patterns_to_output_param);
    opt_parse.add_opt("motifs", '\0', "input is motifs"
                      " (for old style motif files)",false, processing_motifs);
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
  /****************************************************************/
    vector<Pattern*> patterns;
    vector<Motif> motifs;
    patterns_file = input_filenames[0];
    if (processing_motifs) {
      motifs = Motif::ReadMotifVector(patterns_file.c_str());
      for (size_t i = 0; i < motifs.size(); ++i)
	patterns.push_back(&motifs[i]);
    }
    else patterns = PatternFactory::ReadPatternVector(patterns_file.c_str());
    
    /* HERE IS THE SORT FOR MOTIFS */
    PatternOrder po(key.c_str(), reverse_order, numeric_sort_order);
    std::stable_sort(patterns.begin(), patterns.end(), po);

    if (!cutoff.empty()) {
      PatternCutoff pc(key.c_str(), cutoff.c_str(), !reverse_order, numeric_sort_order);
      patterns.erase(std::find_if(patterns.begin(), patterns.end(), pc),
		     patterns.end());
    }
    
    if (assign_rank) AssignRanks(patterns);
    
    ostream* patternout;
    if (!outfile.empty()) patternout = new ofstream(outfile.c_str());
    else patternout = &cout;
    
    const size_t max_patterns_to_output = 
      (max_patterns_to_output_param < 1) ? 
      patterns.size() : static_cast<size_t>(max_patterns_to_output_param);
    
    for (size_t i = 0; i < min(max_patterns_to_output, patterns.size()); ++i)
      *patternout << *(patterns[i]) << endl;
    if (patternout != &cout) delete patternout;
  }
  catch (CREADException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
