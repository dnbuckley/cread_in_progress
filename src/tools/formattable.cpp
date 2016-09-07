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
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>

#include "cread.hpp"
#include "FastaFile.hpp"
#include "WordTable.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"

using std::string;
using std::vector;
using std::map;
using std::cerr;
using std::endl;
using std::count_if;

int main(int argc, const char *argv[]) {
  
/* INPUT PARAMETERS */
string outfile;
string seqfile; // file containing sequences

/* DEFAULT PARAMETERS */
size_t wordsize = 6;
size_t gapsize = 5;

  try {
    /***************** COMMAND LINE OPTIONS  *******************/

    OptionParser opt_parse(strip_path(argv[0]), " ",
                           " ");
    opt_parse.add_opt("wordsize", 'w', "size of words to tabulate (MUST BE EVEN)",
                      false, wordsize);
    opt_parse.add_opt("gapsize", 'g', "max size of gap",
                      false, gapsize);
    opt_parse.add_opt("output", 'o', "output file name",
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
    // READ IN THE SEQUENCES AND THEIR NAMES
    if (wordsize % 2 != 0)
      throw FastaFileException("ERROR: wordsize must be even");
    size_t gapstart = wordsize/2;
   
    seqfile = input_filenames[0]; 
    WordTable wt(seqfile.c_str(), wordsize, gapstart, gapsize);
    wt.WriteWordTable(outfile.c_str());
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
