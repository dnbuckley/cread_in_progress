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
#include "MatCompMethods.hpp"
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::ofstream;
using std::ostream;
using std::ostream_iterator;
using std::numeric_limits;
using std::endl;
using std::cout;
using std::cerr;
using std::min;

float threshold = 1.0;
size_t max_overhang = 0;

void combine_matrices(Matrix& A, Matrix& B, Matrix& C) {
  const Matrix *longer = (A.get_width() > B.get_width()) ? &A : &B;
  const Matrix *shorter = (A.get_width() > B.get_width()) ? &B : &A;
  int offset = 0;
  MatCompMethods::sliding_divergence_offset(*longer, *shorter,
					    offset, max_overhang);
  if (offset > 0) C = Matrix::combine(*longer, *shorter, offset);
  else C = Matrix::combine(*longer, *shorter, -offset);
}

int main(int argc, const char *argv[]) {
  
  string outfile;
  string motifsfile;     
  try {

    OptionParser opt_parse(strip_path(argv[0]), " ",
                           " ");
    opt_parse.add_opt("output",'o',"output file (default:stdout)", 
                       false, outfile);
    opt_parse.add_opt("overhang",'h',"greatest overhang when comparing motifs",
                       false, max_overhang);
    opt_parse.add_opt("threshold",'t',"threshold divergence below" 
                       "which motifs may be joined", false, threshold);

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

    vector<Motif> motifs = Motif::ReadMotifVector(motifsfile.c_str());
    vector<Matrix> matrix, rcmatrix;
    for (vector<Motif>::iterator i = motifs.begin(); i != motifs.end(); ++i) {
      matrix.push_back(i->const_get_matrix().freqmat());
      rcmatrix.push_back(matrix.back().revcomp());
    }
    while (true) {
      float min_diverg = numeric_limits<float>::max();
      size_t matrix1 = 0, matrix2 = 0;
      bool forward = true;
      for (size_t i = 0; i < matrix.size() - 1; ++i)
	for (size_t j = i + 1; j < matrix.size(); ++j) {
	  float temp_diverg =
	    MatCompMethods::sliding_divergence(matrix[i], matrix[j], max_overhang);
	  if (temp_diverg < threshold && temp_diverg < min_diverg) {
	    min_diverg = temp_diverg;
	    matrix1 = i;
	    matrix2 = j;
	    forward = true;
	  }
	  temp_diverg =
	    MatCompMethods::sliding_divergence(matrix[i], matrix[j], max_overhang);
	  if (temp_diverg < threshold && temp_diverg < min_diverg) {
	    min_diverg = temp_diverg;
	    matrix1 = i;
	    matrix2 = j;
	    forward = false;
	  }
	}
      if (min_diverg == numeric_limits<float>::max())
	break;
      Matrix newmat;
      if (forward)
	combine_matrices(matrix[matrix1], matrix[matrix2], newmat);
      else combine_matrices(matrix[matrix1], rcmatrix[matrix2], newmat);
      // Add the combined motif
      matrix.push_back(newmat);
      rcmatrix.push_back(matrix.back().revcomp());
      // Erase the motifs that were combined from the list to be
      // combined
      matrix.erase(matrix.begin() + matrix2);
      matrix.erase(matrix.begin() + matrix1);
      rcmatrix.erase(rcmatrix.begin() + matrix2);
      rcmatrix.erase(rcmatrix.begin() + matrix1);
    }
    vector<Motif> newmotifs;
    for (size_t i = 0; i < matrix.size(); ++i)
      newmotifs.push_back(Motif(matrix[i]));
    
    // output the new motifs
    ostream* output = (outfile.c_str()) ? new ofstream(outfile.c_str()) : &cout;
    copy(newmotifs.begin(), newmotifs.end(), 
	 ostream_iterator<Motif>(*output, "\n"));
    if (output != &cout) delete output;
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
