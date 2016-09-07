/*
 * Copyright (C) 2006-2015 Cold Spring Harbor Laboratory
 *                         University of Southern California
 *                         Authors: Andrew D. Smith, Pavel Sumazin
 *                                  and Michael Q. Zhang
 *
 * Authors: Andrew D. Smith and Sarah S. Ma
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
#include "OptionParser.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;


static void
compute_base_comp(const string &infile, vector<double> &basecomp) {
  
  std::ifstream in(infile.c_str());
  if (!in)
    throw SMITHLABException("cannot open file: " + infile);
  
  string line;
  while (getline(in, line)) {
    if (line[0] != '>')
      for (string::const_iterator i(line.begin()); i != line.end(); ++i)
        if (isvalid(*i))
          ++basecomp[base2int(*i)];
  }
  
  const double total = 
    std::accumulate(basecomp.begin(), basecomp.end(), 0.0);
  
  std::transform(basecomp.begin(), basecomp.end(), basecomp.begin(),
                 std::bind2nd(std::divides<double>(), total));
}



static void
abbreviate_precision(const double decimal_places, vector<double> &basecomp) {
  
  const double shifting_factor = std::pow(10.0, decimal_places);

  transform(basecomp.begin(), basecomp.end(), basecomp.begin(),
            std::bind2nd(std::multiplies<double>(), shifting_factor));
  
  transform(basecomp.begin(), basecomp.end(), basecomp.begin(), 
            std::ptr_fun(round));

  transform(basecomp.begin(), basecomp.end(), basecomp.begin(),
            std::bind2nd(std::divides<double>(), shifting_factor));
}

static string
format_basecomp(const bool FULL_PRECISION, 
                const size_t decimal_places, vector<double> &basecomp) {
  if (FULL_PRECISION) 
    abbreviate_precision(decimal_places, basecomp);
  std::ostringstream oss;
  copy(basecomp.begin(), basecomp.end(), 
       std::ostream_iterator<double>(oss, "\t"));
  return oss.str();
}


int
main(int argc, const char **argv) {
  
  try {

    static const size_t decimal_places = 3;
    
    string outfile;

    // run mode flags
    bool VERBOSE = false;
    bool FILES_SEPARATELY = false;
    bool FULL_PRECISION = false;
    /// bool SEQUENCES_SEPARATELY = false; 
    /// bool SINGLE_STRAND = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "base composition of sequences",
                           "fasta1 [fasta2 ...]");
    opt_parse.add_opt("outfile", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("sep", 's', "compute files separately (default: together)",
                      false, FILES_SEPARATELY);
    opt_parse.add_opt("precise", 'p', "output full precision",
                      false, FULL_PRECISION);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    
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
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    vector<double> basecomp(smithlab::alphabet_size, 0.0);
    
    for (size_t i = 0; i < input_filenames.size(); ++i) {
      
      compute_base_comp(input_filenames[i], basecomp);
      
      if (FILES_SEPARATELY) {
        out << strip_path(input_filenames[i]) << '\t'
            << format_basecomp(FULL_PRECISION, decimal_places, basecomp) << endl;
        basecomp = vector<double>(smithlab::alphabet_size, 0.0);
      }
    }
    
    if (!FILES_SEPARATELY)
      out << format_basecomp(FULL_PRECISION, decimal_places, basecomp) << endl;
    
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
