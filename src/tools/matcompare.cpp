/*
 * Copyright (C) 2006 Cold Spring Harbor Laboratory
 * Authors: Andrew D. Smith, Pavel Sumazin and Michael Q. Zhang
 * Copyright (C) 2006 Cold Spring Harbor Laboratory
 * Authors: Dustin E. Schones, Andrew D. Smith, Pavel Sumazin and Michael Q. Zhang
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
#include "Module.hpp"
#include "MatCompMethods.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"

using std::string;
using std::vector;
using std::pair;
using std::ostringstream;
using std::ostream_iterator;
using std::ofstream;
using std::ostream;
using std::endl;
using std::cout;
using std::cerr;
using std::min;
using std::max;
using std::priority_queue;
using std::numeric_limits;
using std::setw;
using std::ios_base;
using std::greater;
using std::less;

size_t n_top = 1;
int max_overhang = 0;
float threshold = 1.0;

string alignmentfile;

const char *null_factor_name = "NONE";
const char *null_accession = "NONE";
const char *match_label = "MATCH";
const char *bad_name_characters = ": \t";
const char replacement = '_';
const char *separator = ":";

int use_fisher = false;
int use_chisquared = false;

const char *null_name = "NONE";

static const float prior = 0.0005;                
static const size_t n_degenerate_nucleotides = 15;
static const float fixed_matrix[n_degenerate_nucleotides][alphabet_size] = {
  { 1.000, 0.000, 0.000, 0.000 },  // A
  { 0.000, 1.000, 0.000, 0.000 },  // C
  { 0.000, 0.000, 1.000, 0.000 },  // G
  { 0.000, 0.000, 0.000, 1.000 },  // T
  { 0.500, 0.500, 0.000, 0.000 },  // M
  { 0.500, 0.000, 0.500, 0.000 },  // R
  { 0.500, 0.000, 0.000, 0.500 },  // W
  { 0.000, 0.500, 0.500, 0.000 },  // S
  { 0.000, 0.500, 0.000, 0.500 },  // Y
  { 0.000, 0.000, 0.500, 0.500 },  // K
  { 0.333, 0.333, 0.333, 0.000 },  // V
  { 0.333, 0.333, 0.000, 0.333 },  // H
  { 0.333, 0.000, 0.333, 0.333 },  // D
  { 0.000, 0.333, 0.333, 0.333 },  // B
  { 0.250, 0.250, 0.250, 0.250 }   // N
};

static const char fixed_matrix_alphabet[] =
  {'A','C','G','T', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N'};



struct match {
  size_t index;
  float divergence;
  bool rc;
  int offset;
  match(size_t &i, float &d, bool &r, int &o) :
    index(i), divergence(d), rc(r), offset(o) {}
  bool operator<(const match &b) const {
    return divergence < b.divergence;
  }
  bool operator>(const match &b) const {
    return divergence > b.divergence;
  }
};

/**********************************************************************
 * OBTAIN THE IUPAC DEGENERATE NUCLEOTIDE SYMBOL THAT IS CLOSEST TO A
 * PARTICULAR COLUMN TYPE.
 **********************************************************************/
string get_consensus(Matrix matrix) {
  string r;
  for (size_t i = 0; i < matrix.get_width(); ++i) {
    size_t best = 0;
    float score = numeric_limits<float>::max();
    for (size_t j = 0; j < n_degenerate_nucleotides; ++j) {
      float temp_score = 0;
      for (size_t k = 0; k < alphabet_size; k++) 
	if (fixed_matrix[j][k] > 0) {
	  if (matrix[i][k] > 0) 
	    temp_score += (fixed_matrix[j][k] - matrix[i][k]) * 
	      log(fixed_matrix[j][k]/matrix[i][k]);
	  else temp_score += (fixed_matrix[i][j] - prior) * 
		 log(fixed_matrix[i][j]/prior);
	}
	else if (matrix[i][k] > 0) 
	  temp_score += (prior - matrix[i][k]) * log(prior/matrix[i][k]);
      if (temp_score < score) {
	best = j;
	score = temp_score;
      }
    }
    r += fixed_matrix_alphabet[best];
  }
  return r;
}


struct row_data {
  static const size_t n_fields = 6;
  static const size_t spacer_width = 2;
  static const char *spacer;
  static const char *field_names[];
  vector<string> fields;
  row_data(string s) {
    fields = vector<string>(n_fields);
    fields.back() = s;
  }
  row_data(vector<string> &f) : fields(f) {}
  string& operator[](int i) {return fields[i];}
  string format(vector<size_t> fw) {
    ostringstream row;
    for (size_t i = 0; i < fields.size(); ++i)
      row << setw(fw[i]) << fields[i] << spacer;
    return row.str();
  }
  static string get_headers(vector<size_t> fw) {
    ostringstream header;
    for (size_t i = 0; i < n_fields; ++i)
      header << setw(fw[i]) << field_names[i] << spacer;
    return header.str();
  }
  static vector<size_t> get_fw(vector<row_data> &rows) {
    vector<size_t> r;
    for (size_t i = 0; i < n_fields; ++i)
      r.push_back(strlen(field_names[i]));
    for (size_t i = 0; i < rows.size(); ++i)
      for (size_t j = 0; j < n_fields; ++j)
	r[j] = max(r[j], rows[i][j].length());
    return r;
  }
  static string get_separator(vector<size_t> fw) {
    return string(std::accumulate(fw.begin(), fw.end(), 0), '=') +
      string((fw.size() - 1) * spacer_width, '=');
  }
};

const char *
row_data::spacer = "  ";

const char *
row_data::field_names[n_fields] = {
  "Rank", "Accession", "Identifier", 
  "RC", "Divergence", "Alignment"
};

void write_alignment(vector<vector<row_data> > &rows, vector<Motif> &motifs) {
  ofstream table(alignmentfile.c_str());
  for (size_t i = 0; i < motifs.size(); ++i)
    if (!rows[i].empty()) {
      vector<size_t> fw = row_data::get_fw(rows[i]);
      string line = row_data::get_separator(fw);
      table << motifs[i].get_accession() << endl << line << endl
	    << row_data::get_headers(fw) << endl << line << endl;
      for (size_t j = 0; j < rows[i].size(); ++j) {
	table << rows[i][j].format(fw) << endl;
	if (j < rows[i].size() - 1 && j % 2) table << endl;
      }
      table << line << endl << endl;
    }
  table.close();
}

void get_alignment_rows(Motif &motif, vector<Motif> &motif_lib, 
			vector<match> &top_matches, vector<row_data> &rows) {
  for (vector<match>::iterator i = top_matches.begin(); 
       i != top_matches.end(); ++i) {
    vector<string> temp_row;
    Motif lm = motif_lib[i->index];
    Matrix libmat = lm.const_get_matrix().freqmat();
    Matrix matrix = motif.const_get_matrix().freqmat();
    if (i->rc) matrix = matrix.revcomp();
    string consensus = get_consensus(matrix);
    string lib_consensus = get_consensus(libmat);
    string filler(std::abs(i->offset), ' ');
    if ((i->offset >= 0) == (libmat.get_width() <= matrix.get_width()))
      lib_consensus = filler + lib_consensus;
    else consensus = filler + consensus;
    int right_fill = static_cast<int>(consensus.length()) - 
      static_cast<int>(lib_consensus.length());
    if (right_fill > 0)
      fill_n(back_inserter(lib_consensus), right_fill, ' ');
    else fill_n(back_inserter(consensus), -right_fill, ' ');
    temp_row.push_back(cread::toa(i - top_matches.begin() + 1));
    temp_row.push_back(lm.get_accession());
    temp_row.push_back(lm.get_identifier());
    temp_row.push_back(cread::toa(i->rc));
    temp_row.push_back(cread::toa(i->divergence));
    temp_row.push_back(lib_consensus);
    rows.push_back(row_data(consensus));
    rows.push_back(row_data(temp_row));
  }
}

string replace_bad_characters(string s) {
  string r = s;
  size_t length = strlen(bad_name_characters);
  for (size_t i = 0; i < length; ++i)
    replace(r.begin(), r.end(), bad_name_characters[i], replacement);
  return r;
}

void add_factors(Motif &motif, vector<Motif> &library,
		 vector<match> &top_matches) {
  string match_prefix = match_label;
  for (size_t i = 0; i < top_matches.size(); ++i) {
    Motif match = library[top_matches[i].index];
    string match_value = match.get_accession() +
      string(separator) + replace_bad_characters(match.get_factor_names()) +
      string(separator) + match.get_identifier() + 
      string(separator) + cread::toa(top_matches[i].divergence);
    motif.set_attribute(match_prefix + cread::toa(i + 1), match_value);
  }
}

void add_factors_module(Module &module, vector<Motif> &library,
			vector<vector<match> > &top_matches) {
  string match_prefix = match_label;
  for (size_t i = 0; i < module.size(); ++i) {
    string matrix_prefix = Module::add_index(match_prefix.c_str(), i);
    for (size_t j = 0; j < top_matches[i].size(); ++j) {
      Motif match = library[top_matches[i][j].index];
      string match_value = match.get_accession() + 
 	string(separator) + match.get_identifier() + 
 	string(separator) + replace_bad_characters(match.get_factor_names()) +
	string(separator) + cread::toa(top_matches[i][j].divergence);
      module.set_attribute(Module::add_index(matrix_prefix.c_str(), j + 1),
			   match_value);
    }
  }
}


void get_matches(Matrix matrix, vector<Motif> &lib, 
		 vector<match> &matches) {
  size_t min_width = matrix.get_width() - max_overhang;
  const Matrix matrix_rc = matrix.revcomp();
  
  typedef vector<Motif>::iterator motif_iter;
  
  if (use_fisher || use_chisquared){
    priority_queue<match, vector<match>, greater<match> > top_matches;
    size_t index = 0;
    for (motif_iter j = lib.begin(); j != lib.end(); ++j){
      if (j->get_width() >= min_width){
	Matrix b = j->get_matrix();
	if (use_fisher){
	  // do fisher test
	  float fisher_forward = 
	    MatCompMethods::sliding_fishertest(b, matrix, max_overhang);
	  float fisher_rc = 
	    MatCompMethods::sliding_fishertest(b, matrix_rc, max_overhang);
	  bool rc = false;
	  float fisher;
	  int offset = 0;
	  if(fisher_rc > fisher_forward){
	    fisher = fisher_rc;
	    rc = true;
	    offset = MatCompMethods::sliding_fishertest_offset(b, matrix_rc, 
							   max_overhang);
	  }
	  else{
	    fisher = fisher_forward;
	    offset = MatCompMethods::sliding_fishertest_offset(b, matrix, 
							 max_overhang);
	  }
	  if (fisher > threshold &&
	      (top_matches.size() < n_top ||
	       fisher > top_matches.top().divergence)){
	    top_matches.push(match(index, fisher, rc, offset));
	    if (top_matches.size() > n_top)
	      top_matches.pop();
	  }
	}
	
#ifdef HAVE_GSL
	else if (use_chisquared){
	  // do chi squared test
	  float chisquare_forward = 
	    MatCompMethods::sliding_chisquared(b, matrix, max_overhang);
	  float chisquare_rc = 
	    MatCompMethods::sliding_chisquared(b, matrix_rc, max_overhang);
	  bool rc = false;
	  float chisquare;
	  int offset = 0;
	  if (chisquare_rc > chisquare_forward){
	    chisquare = chisquare_rc;
	    rc = true;
	    offset = MatCompMethods::sliding_chisquared_offset(b, matrix_rc, 
							   max_overhang);
	  }
	  else{
	    chisquare = chisquare_forward;
	    offset = MatCompMethods::sliding_chisquared_offset(b, matrix, 
							   max_overhang);
	  }
	  if (chisquare > threshold &&
	      (top_matches.size() < n_top ||
	       chisquare > top_matches.top().divergence)){
	    top_matches.push(match(index, chisquare, rc, offset));
	    if (top_matches.size() > n_top)
	      top_matches.pop();
	  }
	}
#endif
      }
      index++;
    }
    
    
    matches.clear();
    while (top_matches.size() > 0) {
      matches.push_back(top_matches.top());
      top_matches.pop();
    }
    std::reverse(matches.begin(), matches.end());
  }

  else{
    // do divergence test
    priority_queue<match> top_matches;
    size_t index = 0;

    //TODO: changes to functions in MatCompMethods mean that most of
    //the code below is really mis-using those methods. In particular,
    //the work is actually done twice: once by the calls to
    //"sliding_divergence" and again by the calls to
    //"sliding_divergence_offset", which now produces the divergences
    //as a return value, and passes the offset back through a
    //parameter.

    for (motif_iter j = lib.begin(); j != lib.end(); ++j){
      if (j->get_width() >= min_width){
	Matrix b = j->get_matrix().freqmat();
	float divergence_forward = 
	  MatCompMethods::sliding_divergence(matrix, b.freqmat(), max_overhang);
	float divergence_rc = 
	  MatCompMethods::sliding_divergence(matrix_rc, b.freqmat(), max_overhang);
	bool rc = false;
	float divergence;
	int offset = 0;
	if (divergence_rc < divergence_forward) {
	  divergence = divergence_rc;
	  rc = true;
	  MatCompMethods::sliding_divergence_offset(b, matrix_rc,
						    offset, max_overhang);
	}
	else {
	  divergence = divergence_forward;
	  MatCompMethods::sliding_divergence_offset(b, matrix,
						    offset, max_overhang);
	}
	if (divergence < threshold && 
	    (top_matches.size() < n_top || 
	     divergence < top_matches.top().divergence)) {
	  top_matches.push(match(index, divergence, rc, offset));
	  if (top_matches.size() > n_top){
	    top_matches.pop();
	    
	  }
	}
      }
      index++;
    }
    matches.clear();
    while (top_matches.size() > 0) {
      matches.push_back(top_matches.top());
      top_matches.pop();
    }
    std::reverse(matches.begin(), matches.end());
  }
}


Matrix get_freqmat(Motif m) {
  return m.const_get_matrix().freqmat();
}



int main(int argc, const char **argv) {
  int processing_modules = false;
  int list_mode = false;
  int verbose = false;  
  string outfile;
  string motifsfile;
  string libfile;
  string base_comp_string;
  try {
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), " ",
                           " ");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                       false, outfile);
    opt_parse.add_opt("alignment-file", 'a',  "alignment output file"
                      " (default: None)", false, alignmentfile);
    opt_parse.add_opt("overhang", 'h', "greatest overhang when "
                      "comparing motifs", false, max_overhang);
    opt_parse.add_opt("lib", 'l', "motif lib file name",
                       false, libfile);
    opt_parse.add_opt("threshold", 't',"threshold divergence"
                      " for elimination", false, threshold);
    opt_parse.add_opt("modules", 'M',  "specifies that the file"
                      " contains modules", false, processing_modules);
    opt_parse.add_opt("ntop", 'n',  "number of top matches to report", 
                       false, n_top);
    opt_parse.add_opt("base-comp", 'b',  "comma separated base frequencies",
                       false, base_comp_string);
    opt_parse.add_opt("verbose", 'v', "print more match info in attributes",
                       false, verbose);
    opt_parse.add_opt("fisher", 'F', "default: K-L divergence", 
                       false, use_fisher);
    opt_parse.add_opt("chisquared", 'C', "default: K-L divergence",
                       false, use_chisquared);
    opt_parse.add_opt("list", 'L', "default: annotation mode",
                       false, list_mode);

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
    if (list_mode) {
      vector<Motif> library = Motif::ReadMotifVector(libfile.c_str());
      vector<Motif> motifs;
      vector<Module> modules;
      
      ostringstream outbuffer;
      if (!processing_modules){
	motifs = Motif::ReadMotifVector(motifsfile.c_str());

	typedef vector<Motif>::iterator motif_iter;
	for (motif_iter motif = motifs.begin(); motif != motifs.end(); ++motif) {
	  float best_score = (!use_chisquared && !use_fisher) ? 
	    numeric_limits<float>::max() : 0;
	  Matrix a(!use_chisquared && !use_fisher ? 
		   motif->const_get_matrix().freqmat() : 
		   motif->const_get_matrix());
	  string motif_name = motif->get_accession();
	  string best_name = null_name;
	  size_t min_width = a.get_width() - max_overhang;
	  for (motif_iter j = library.begin(); j != library.end(); ++j)
	    if (j->get_width() >= min_width) {
	      Matrix b(j->get_matrix());
	      if (use_fisher) {
		float fisher = MatCompMethods::sliding_fishertest(a, b, max_overhang);
		if (fisher > best_score) {
		  best_score = fisher;
		  best_name = j->get_accession();
		}
	      }
#ifdef HAVE_GSL
	      else if (use_chisquared) {
		float chi = MatCompMethods::sliding_chisquared(a, b, max_overhang);
		if (chi > best_score) {
		  best_score = chi;
		  best_name = j->get_accession();
		}
	      }
#endif	      
	      else {
		b = b.freqmat();
		float divergence = MatCompMethods::sliding_divergence(a, b, max_overhang);
		if (divergence < best_score) {
		  best_score = divergence;
		  best_name = j->get_accession();
		}
	      }
	    }
	  outbuffer << motif_name << " " << best_name << " " << best_score << endl;
	}
      }
      else{
	modules = Module::ReadModuleVector(motifsfile.c_str());
	for (size_t i = 0; i < modules.size(); ++i) {
	  for (size_t j = 0; j < modules[i].size(); ++j){
	    float best_score = (!use_chisquared && !use_fisher) ? 
	      numeric_limits<float>::max() : 0;
	    Matrix a = (!use_chisquared && !use_fisher) ? 
	      modules[i][j].freqmat() : modules[i][j];
	    string motif_name(modules[i].get_accession() + 
			      "." + cread::toa(j + 1));
	    string best_name = null_name;
	    size_t min_width = a.get_width() - max_overhang;
	    typedef vector<Motif>::iterator motif_iter;
	    for (motif_iter k = library.begin(); k != library.end(); ++k)
	      if (k->get_width() >= min_width){
		Matrix b(k->get_matrix());
		if (use_fisher) {
		  float fisher = MatCompMethods::sliding_fishertest(a, b, max_overhang);
		  if (fisher > best_score) {
		    best_score = fisher;
		    best_name = k->get_accession();
		  }
		}
#ifdef HAVE_GSL
		else if (use_chisquared) {
		  float chi = MatCompMethods::sliding_chisquared(a, b, max_overhang);
		  if (chi > best_score) {
		    best_score = chi;
		     best_name = k->get_accession();
		  }
		}
#endif
		else {
		  b = b.freqmat();
		  float divergence = MatCompMethods::sliding_divergence(a, b, max_overhang);
		  if (divergence < best_score) {
		    best_score = divergence;
		    best_name = k->get_accession();
		  }
		}
	      }
	    outbuffer << motif_name << " " << best_name << " " << best_score << endl;
	  }
	}
      }
      ostream* output = (outfile.c_str()) ? new ofstream(outfile.c_str()) : &cout;
      *output << outbuffer.str();
      if (output != &cout) delete output;
    }
    else {
      // use annotation mode
      float base_comp[alphabet_size];
      if (base_comp_string.c_str()) {
	vector<string> parts = cread::split(base_comp_string.c_str(), ",");
	if (parts.size() != alphabet_size) {
	  cerr << "ERROR: incorrect base composition format: " 
	       << base_comp_string.c_str() << endl;
	  return EXIT_FAILURE;
	}
	for (size_t i = 0; i < parts.size(); ++i)
	  base_comp[i] = std::atof(parts[i].c_str());
      }
      else std::fill_n(base_comp, alphabet_size, 1.0/alphabet_size);
    
      vector<Motif> motif_lib = Motif::ReadMotifVector(libfile.c_str());
      vector<Matrix> lib;
      transform(motif_lib.begin(), motif_lib.end(), back_inserter(lib),
		std::ptr_fun(&get_freqmat));
    
      vector<Motif> motifs;
      vector<Module> modules;
      if (!processing_modules) {
	vector<vector<row_data> > rows;
	motifs = Motif::ReadMotifVector(motifsfile.c_str());
	typedef vector<Motif>::iterator motif_iter;
	for (motif_iter motif = motifs.begin(); motif != motifs.end(); ++motif){
	  Matrix a = (!use_chisquared && !use_fisher) ? 
	    motif->const_get_matrix().freqmat() : motif->const_get_matrix();
	  vector<match> top_matches;
	  get_matches(a, motif_lib, top_matches);
	  add_factors(*motif, motif_lib, top_matches);
	  rows.push_back(vector<row_data>());
	  if (alignmentfile.c_str() && !top_matches.empty())
	    get_alignment_rows(*motif, motif_lib, top_matches, rows.back());
	}
	if (alignmentfile.c_str()) write_alignment(rows, motifs);
      }
      else {
	// doing modules
	modules = Module::ReadModuleVector(motifsfile.c_str());
	for (size_t i = 0; i < modules.size(); ++i) {
	  vector<vector<match> > top_matches;
	  for (size_t j = 0; j < modules[i].size(); ++j) {
	    Matrix a = (!use_chisquared && !use_fisher) ?
	      modules[i][j].freqmat() : modules[i][j];  
	    top_matches.push_back(vector<match>());
	    get_matches(a, motif_lib, top_matches.back());
	  }
	  add_factors_module(modules[i], motif_lib, top_matches);
	}
      }
      
      ostream* out = (outfile.c_str()) ? new ofstream(outfile.c_str()) : &cout;
      if (processing_modules)
      	copy(modules.begin(), modules.end(), ostream_iterator<Module>(*out, "\n"));
      else copy(motifs.begin(), motifs.end(), ostream_iterator<Motif>(*out, "\n"));
      if (out != &cout) delete out;
    }
    
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
