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
#include "Motif.hpp"
#include "FastaFile.hpp"
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::fill_n;
using std::copy;
using std::map;
using std::ofstream;
using std::ostream;
using std::endl;
using std::cout;
using std::cerr;
using std::numeric_limits;

int main(int argc, const char **argv) {    
  string seqfile;
  string motifsfile;
  string outfile;
  bool mask = false;
  bool invert = false;
  size_t context_size = 0;
  size_t context_start = numeric_limits<size_t>::max();

  // TODO: make this program work for words and modules
  try {
  /****************** COMMAND LINE OPTIONS***************************/
    OptionParser opt_parse(strip_path(argv[0]), " ",
                           " ");
    opt_parse.add_opt("sequences", 'f', "name of sequences file", false, seqfile);
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", false, outfile);
    opt_parse.add_opt("context", 'c', "number of letters of context", false, context_size);
    //opt_parse.add_opt("start", 's', "upstream starting position of context", context_start);
    opt_parse.add_opt("mask", 'm', "mask the occurrences", false, mask);
    opt_parse.add_opt("invert", 'i', "get all but the context", false, invert);

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
  /*******************END COMMAND LINE OPTIONS***************************/
   
    // READ SEQUENCES
    motifsfile = input_filenames[0];
    FastaFile faa(seqfile.c_str());
    vector<string> names(faa.get_names());
    vector<string> sequences(faa.get_sequences());
    
    map<string, size_t> seq_index;
    for (size_t i = 0; i < names.size(); ++i)
      seq_index[names[i]] = i;
    
    // Read in the motifs
    vector<Motif> motifs = Motif::ReadMotifVector(motifsfile.c_str());
    vector<string> new_sequences;
    vector<string> new_names;
    map<string, size_t> used_names;
    for (size_t i = 0; i < motifs.size(); ++i) {
      vector<MotifSite> bs = motifs[i].get_sites();
      for (size_t j = 0; j < bs.size(); ++j)
	if (seq_index.find(bs[j].get_seq_name()) != seq_index.end()) {
	  size_t index = seq_index.find(bs[j].get_seq_name())->second;
	  size_t length = bs[j].get_length();
	  size_t offset = bs[j].get_start();
	  size_t true_offset = (offset > context_start) ?
	    offset - context_start : 0;
	  size_t true_size = std::min(context_size,
				      sequences[index].length() -
				      true_offset) + length;
	  if (invert) {
	    new_sequences.push_back(sequences[index]);
	    fill_n(new_sequences.back().begin() + true_offset, true_size, 'N');
	  }
	  else
	    new_sequences.push_back(sequences[index].substr(true_offset, 
							    true_size));
	  if (mask) {
	    size_t mask_start = std::min(context_start, offset);
	    fill_n(new_sequences.back().begin() + mask_start, length, 'N');
	  }
	  string new_name = names[index] + string(":") +
	    motifs[i].get_accession();
	  if (used_names.find(new_name) != used_names.end()) {
	    used_names[new_name]++;
	    new_name += string(":") + cread::toa(used_names[new_name]);
	  }
	  else used_names[new_name] = 0;
	  new_names.push_back(new_name);
	}
    }
    
    ostream* seqout = ((!outfile.empty()) ? new ofstream(outfile.c_str()) : &cout);
    for (size_t i = 0; i < new_sequences.size(); ++i)
      *seqout << ">" << new_names[i] << endl << new_sequences[i] << endl;
    if (seqout != &cout) delete seqout;
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
