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
#include "Module.hpp"
#include "MatCompMethods.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"

using std::string;
using std::vector;
using std::list;
using std::ofstream;
using std::ostream;
using std::ostream_iterator;
using std::bind2nd;
using std::numeric_limits;
using std::ptr_fun;
using std::min;
using std::max;
using std::abs;
using std::cerr;
using std::cout;
using std::endl;

// FILE NAMES

size_t max_overhang = 2;
float threshold = 1.0;
const char *elim_label = "ELIMINATED_BY";

inline bool eliminated(Motif &m) {
  return m.get_attribute(elim_label) != "";
}

bool equiv_motifs(Motif m1, Motif m2) {
  if (eliminated(m1) || eliminated(m2)) return false;
  Matrix a = m1.const_get_matrix().freqmat();
  Matrix b = m2.const_get_matrix().freqmat();
  return (MatCompMethods::sliding_divergence(a, b, max_overhang) < threshold ||
	  MatCompMethods::sliding_divergence(a, b.revcomp(),
					 max_overhang) < threshold);
}

bool equiv_matrices(Matrix m1, Matrix m2) {
  return (MatCompMethods::sliding_divergence(m1, m2, max_overhang) < threshold ||
	  MatCompMethods::sliding_divergence(m1, m2.revcomp(),
					 max_overhang) < threshold);
}

bool equiv_modules(Module m1, Module m2) {
  Module *shorter, *longer;
  if (m1.size() < m2.size()) {
    shorter = &m1;
    longer = &m2;
  }
  else {
    shorter = &m2;
    longer = &m1;
  }
  for (vector<Matrix>::const_iterator i = shorter->begin(); 
       i != shorter->end(); ++i)
    if (find_if(longer->begin(), longer->end(), 
		bind2nd(ptr_fun(equiv_matrices), *i)) == longer->end())
      return false;
  return true;
}

int main(int argc, const char **argv) {

  string outfile;
  string motifsfile;  

  size_t ntop = numeric_limits<size_t>::max();
  int processing_modules = false;
  int keep_eliminated = false;
  string rank_label;
   
  // TODO: make sure the divergence and overhang are sane
  //       and that the program doesn't crash with bad input
  try {
  /***************** COMMAND LINE OPTIONS *******************/
    OptionParser opt_parse(strip_path(argv[0]), " ",
                           " ");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                       false, outfile);
    opt_parse.add_opt("overhang", 'h', "greatest overhang when"
                      " comparing motifs", false, max_overhang);
    opt_parse.add_opt("threshold", 't', "threshold divergence"
                      " for elimination", false, threshold);
    opt_parse.add_opt("modules", 'M', "processing modules", 
                       false, processing_modules);
    opt_parse.add_opt("ntop", 'n', "number of top motifs to outout", 
                       false, ntop);
    opt_parse.add_opt("keep", 'k', "keep motifs that were eliminated;"
                      " just re-order", false, keep_eliminated);
    opt_parse.add_opt("ranks", 'r', "update the ranks for this attribute",
                       false, rank_label);

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

    motifsfile = input_filenames[0];
    if (!processing_modules) {
      vector<Motif> motifs_v = Motif::ReadMotifVector(motifsfile.c_str());
      list<Motif> motifs;
      copy(motifs_v.begin(), motifs_v.end(), back_inserter(motifs));
      if (keep_eliminated)
	for (list<Motif>::iterator i = motifs.begin(); i != motifs.end();) {
	  list<Motif>::iterator found = find_if(motifs.begin(), i,
						bind2nd(ptr_fun(equiv_motifs), *i));
	  if (found == i) ++i;
	  else {
	    i->set_attribute(elim_label, found->get_accession());
	    ++found;
	    motifs.insert(found, *i);
	    i = motifs.erase(i);
	  }
	}
      else
	for (list<Motif>::iterator i = motifs.begin(); i != motifs.end();)
	  if (find_if(motifs.begin(), i,
		      bind2nd(ptr_fun(equiv_motifs), *i)) == i) ++i;
	  else i = motifs.erase(i);
      ostream* motifout = (outfile.c_str()) ? new ofstream(outfile.c_str()) : &cout;
      size_t rank = 0;
      for (list<Motif>::iterator i = motifs.begin(); i != motifs.end(); ++i) {
	if (i->get_attribute(elim_label) == "" && rank++ == ntop)
	  break;
	if (rank_label.c_str()){
	  if (i->get_attribute(elim_label) == "")
	    i->set_attribute(rank_label.c_str() + string("_RANK"), cread::toa(rank));
	  else i->set_attribute(rank_label.c_str() + string("_RANK"), "NA");
        }
	*motifout << *i << endl;
      }
      if (motifout != &cout) delete motifout;
    }
    else {
      vector<Module> modules_v = Module::ReadModuleVector(motifsfile.c_str());
      list<Module> modules;
      copy(modules_v.begin(), modules_v.end(), back_inserter(modules));
      modules_v.clear();
      if (keep_eliminated)
	for (list<Module>::iterator i = modules.begin(); i != modules.end();) {
	  list<Module>::iterator found = find_if(modules.begin(), i,
						 bind2nd(ptr_fun(equiv_modules), *i));
	  if (found == i) ++i;
	  else {
	    i->set_attribute(elim_label, found->get_accession());
	    ++found;
	    modules.insert(found, *i);
	    i = modules.erase(i);
	  }
	}
      else
	for (list<Module>::iterator i = modules.begin(); i != modules.end();)
	  if (find_if(modules.begin(), i, 
		      bind2nd(ptr_fun(equiv_modules), *i)) == i) ++i;
	  else i = modules.erase(i);
      ostream* moduleout = (outfile.c_str()) ? new ofstream(outfile.c_str()) : &cout;
      size_t rank = 0;
      for (list<Module>::iterator i = modules.begin(); i != modules.end(); ++i) {
	if (i->get_attribute(elim_label) == "" && rank++ == ntop)
	  break;
	if (rank_label.c_str()) {
	  if (i->get_attribute(elim_label) == "")
	    i->set_attribute(rank_label.c_str() + string("_RANK"), cread::toa(rank));
      	  else i->set_attribute(rank_label.c_str() + string("_RANK"), "NA");
        }
	*moduleout << *i << endl;
      }
      if (moduleout != &cout) delete moduleout;
    }
  }
  // TODO: there must be more exceptions to catch
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
