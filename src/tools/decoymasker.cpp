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
#include "FastaFile.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"

using std::string;
using std::vector;
using std::ostream;
using std::ofstream;
using std::ostringstream;
using std::pair;
using std::cout;
using std::cerr;
using std::endl;


const char masked_base = 'N';

string mask(vector<string> &seqs, vector<string> &names, size_t width, 
	    size_t max_order) {
  ostringstream log;
  for (size_t i = 1; i <= max_order; ++i)
    for (size_t j = 0; j < seqs.size(); ++j) {
      size_t rep = i;
      for (size_t k = i; k < seqs[j].length(); ++k)
	if (seqs[j][k] == seqs[j][k - i] && seqs[j][k] != masked_base &&
	    k < seqs[j].length() - 1)
	  rep++;
	else {
	  if (rep > width) {
	    log << names[j] << ":\t" << seqs[j].substr(k - rep, rep) << endl;
	    fill_n(seqs[j].begin() + k - rep, rep, masked_base);
	  }
	  rep = i;
	}
    }
  return log.str();
}

int main(int argc, const char **argv) {

  string seqfile;
  string logname;
  string outfile;

  size_t width = 10;
  size_t max_order = 2;

  try {      
    /******************** COMMAND LINE OPTIONS ************************/

    OptionParser opt_parse(strip_path(argv[0]), " ",
                           " ");
    opt_parse.add_opt("width", 'w', "minimum width of decoys",
                       false, width);
    opt_parse.add_opt("max-order", 'r', "maximum order of repeats",
                       false, max_order);
    opt_parse.add_opt("log-file", 'l', "name of file in which to log masks",
                       false, logname);
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                       false, outfile);

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
    seqfile= input_filenames[0];
    FastaFile faa = FastaFile(seqfile.c_str());
    vector<string> names(faa.get_names());
    vector<string> sequences(faa.get_sequences());
    
    string log = mask(sequences, names, width, max_order);
    
    ostream* seqout = (!outfile.empty()) ? new ofstream(outfile.c_str()) : &cout;
    for (size_t i = 0; i < sequences.size(); ++i)
      *seqout << ">" << names[i] << endl << sequences[i] << endl;
    if (seqout != &cout) delete seqout;
    
    if (logname.c_str()) {
      ofstream logfile(logname.c_str());
      logfile << log;
    }
  }
  catch (CREADException &e) {
    cerr << "ERROR\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
} 
