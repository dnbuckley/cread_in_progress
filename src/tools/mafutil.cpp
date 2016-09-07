/*
 * Copyright (C) 2007 Cold Spring Harbor Laboratory
 * Authors: Andrew D. Smith
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
#include "GenomicRegion.hpp"
#include "GenomeUtils.hpp"
#include "GenomeAlignment.hpp"
#include "BEDFile.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"

#include <errno.h>
#include <dirent.h>

#include <cassert>

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::cout;
using std::numeric_limits;
using std::max;
using std::min;
using std::map;
using std::set;
using std::mem_fun_ref;
using std::ostream_iterator;
using std::ofstream;
using std::ostream;

int VERBOSE = 0;


string 
path_join(const string& a, const string& b) {
  return a + "/" + b;
}


string
strip_suffix(string filename, string suffix) {
  return filename.substr(0, filename.length() - suffix.length() - 1);
}


GenomicRegion
filename2region(const string& filename) {
  string chrom;
  size_t start = 0, end = 0;
  parse_region_name(filename, chrom, start, end);
  if (start == 0 && start == end)
    end = numeric_limits<size_t>::max();
  return GenomicRegion(chrom, start, end, filename, 0.0, '+');
}


bool
is_valid_filename(const string name, const string& filename_suffix) {
  string chrom;
  size_t start = 0, end = 0;
  string basename(name.substr(0, name.find_first_of(".")));
  parse_region_name(basename, chrom, start, end);
  if (name.find(':') != string::npos)
    chrom += ":" + cread::toa(start) + "-" + cread::toa(end);
  return (name == chrom + '.' + filename_suffix);
}


vector<GenomicRegion> 
read_db_dir(const string& dirname, string filename_suffix) {
  DIR *dir;
  if (!(dir = opendir(dirname.c_str()))) {
    cerr << "ERROR: in read_db_dir(): could not open directory: " 
	 << dirname << endl << errno << endl;
    exit(EXIT_FAILURE);
  }
  vector<string> filenames;
  errno = 0;
  struct dirent *ent;
  while ((ent = readdir(dir))) {
    if (is_valid_filename(ent->d_name, filename_suffix))
      filenames.push_back(ent->d_name);
    errno = 0;
  }
  if (errno) {
    cerr << "ERROR: in read_db_dir(): " << errno << endl;
    exit(EXIT_FAILURE);
  }
  if (filenames.empty()) {
    cerr << "ERROR: no valid files found in dir: " << dirname << endl;
    exit(EXIT_FAILURE);
  }
  
  transform(filenames.begin(), filenames.end(),
	    filenames.begin(), std::bind2nd(std::ptr_fun(&strip_suffix),
					    filename_suffix));

  vector<GenomicRegion> file_regions;
  transform(filenames.begin(), filenames.end(),
	    back_inserter(file_regions), std::ptr_fun(&filename2region));
  sort(file_regions.begin(), file_regions.end());
  return file_regions;
}


void 
get_regions_by_filename(const vector<GenomicRegion>& file_regions,
			const vector<GenomicRegion>& regions,
			vector<vector<GenomicRegion> >& regions_by_file) {
  
  regions_by_file = vector<vector<GenomicRegion> >(file_regions.size());
  
  for (size_t i = 0; i < regions.size(); ++i) {
    const size_t file_id = 
      find_closest(file_regions, regions[i]) - file_regions.begin();
    regions_by_file[file_id].push_back(regions[i]);
  }
}

inline static bool
is_space_or_tab(const char &c) {
  return c == ' ' || c == '\t';
}

inline static void
parse_maf_block_line(const char *s,
		     const size_t len,
		     string &spec, 
		     string &chrom, 
		     size_t &start, 
		     size_t &end,
		     char &strand, 
		     string &seq) {
  size_t i = 2;
  
  // s name.chrom start length strand junk sequence
  
  // the species
  while (s[i] != '.' && i < len) ++i;
  if (s[i] != '.')
    throw GenomeAlignmentException("invalid line:\n" + string(s));
  spec = string(s + 2, i - 2);//s + i - i - j);
  ++i;
  // else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // the chromosome
  while (is_space_or_tab(s[i]) && i < len) ++i;
  size_t j = i;
  while (!is_space_or_tab(s[i]) && i < len) ++i;
  if (j < i) chrom = string(s + j, i - j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // start of the region (a positive integer)
  while (is_space_or_tab(s[i]) && i < len) ++i;
  j = i;
  while (!is_space_or_tab(s[i]) && i < len) ++i;
  if (j < i) start = atoi(s + j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // end of the region (a positive integer)
  while (is_space_or_tab(s[i]) && i < len) ++i;
  j = i;
  while (!is_space_or_tab(s[i]) && i < len) ++i;
  if (j < i) 
    end = start + atoi(s + j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // the strand, as single character ('+' or '-')
  while (is_space_or_tab(s[i]) && i < len) ++i;
  j = i;
  while (!is_space_or_tab(s[i]) && i < len) ++i;
  if (j < i) strand = s[j];
  else throw GenomeAlignmentException("invalid line:\n" + string(s));

  // the actual alignment sequence
  while (is_space_or_tab(s[i]) && i < len) ++i;
  j = i;
  while (!is_space_or_tab(s[i]) && i < len) ++i;
  if (j < i)
    string jnk = string(s + j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // the actual alignment sequence
  while (is_space_or_tab(s[i]) && i < len) ++i;
  j = i;
  while (!is_space_or_tab(s[i]) && i < len) ++i;
  if (j < i)
    seq = string(s + j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
}

void
generate_maf_file_index(const char *input_filename,
			const char *output_filename) {
  static const size_t buffer_size = 1000000; // Magic

  std::ostream* out = (output_filename) ? 
    new std::ofstream(output_filename) : &cout;
  
  // open and check the file
  std::ifstream in(input_filename, std::ios_base::binary);
  if (!in) 
    throw CREADException("cannot open input file " + 
			 cread::toa(input_filename));
  
  bool next_is_refspec = false;
  while (!in.eof()) {
    const size_t offset = in.tellg();
    
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw CREADException("Line too long in file: " + 
			   cread::toa(input_filename));
    
    if (buffer[0] == 'a')
      next_is_refspec = true;
    else {
      if (next_is_refspec) {
	if (buffer[0] != 's')
	  throw CREADException("Bad block in file: " + 
			       cread::toa(input_filename));
	else {
	  string spec, chrom, seq;
	  size_t start = 0, end = 0;
	  char strand = '+'; 
	  // assert(static_cast<size_t>(in.gcount()) == strlen(buffer) + 1);
	  parse_maf_block_line(buffer, in.gcount() - 1,
			       spec, chrom, start, end, strand, seq);
	  *out << SimpleGenomicRegion(chrom, start, end) << "\t"
	       << 'x' << "\t"
	       << offset << "\t"
	       << strand << endl;
	}
	next_is_refspec = false;
      }
    }
    in.peek();
  }
  if (out != &cout) delete out;
}


bool
test_index_file(const char *alnfile, const char *indexfile) {

  if (VERBOSE)
    cerr << "testing index file: " 
	 << indexfile 
	 << " on alignment file: "
	 << alnfile << endl;
  
  static const size_t buffer_size = 1000000; // Magic

  vector<GenomicRegion> index;
  ReadBEDFile(indexfile, index);

  vector<size_t> offsets;
  GenomeAlignmentOnDisk::ReadIndexFile(indexfile, offsets);

  assert(index.size() == offsets.size());
  
  // open and check the file
  std::ifstream in(alnfile, std::ios_base::binary);
  if (!in) 
    throw CREADException("cannot open input file " + 
			 cread::toa(alnfile));
  
  for (size_t i = 0; i < index.size(); ++i) {
    in.seekg(offsets[i]);
    
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw CREADException("Line too long in file: " + cread::toa(alnfile));
    
    string spec, chrom, seq;
    size_t start = 0, end = 0;
    char strand = '+'; 
    // assert(static_cast<size_t>(in.gcount()) == strlen(buffer) + 1);
    parse_maf_block_line(buffer, in.gcount() - 1,
			 spec, chrom, start, end, strand, seq);
    if (GenomicRegion(chrom, start, end, "x", offsets[i], strand) != index[i]) {
      if (VERBOSE) {
	cerr << "test failed" << endl
	     << "invalid index at line: " << i << endl;
      }
      return false;
    }
  }
  if (VERBOSE)
    cerr << "test passed" << endl;
  return true;
}


void
read_species_names(const string filename, 
		   set<string> &species_names) {
  static const size_t buffer_size = 1000; // Magic

  // open and check the file
  std::ifstream in(filename.c_str());
  if (!in) 
    throw "cannot open input file " + filename;
  
  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw "Line too long in file: " + filename;
    species_names.insert(species_names.end(), buffer);
    in.peek();
  }
  in.close();

  if (VERBOSE) {
    cerr << "species to keep:" << endl;
    for (set<string>::iterator i(species_names.begin()); 
	 i != species_names.end(); ++i)
      cerr << *i << endl;
  }
}


void
get_species_alignment(const set<string> &species_names,
		      vector<GenomicRegion> &regions,
		      vector<string> &segments) {
  
  vector<GenomicRegion> good_regions;
  vector<string> good_segments;
  
  for (size_t i = 0; i < regions.size(); ++i)
    if (species_names.find(regions[i].get_name()) != species_names.end()) {
      good_regions.push_back(regions[i]);
      good_segments.push_back(segments[i]);
    }

  segments = vector<string>(good_segments.size());
  size_t prev_index = 0;
  for (size_t i = 0; i < good_segments.front().length(); ++i) {
    bool has_non_gap = false;
    for (size_t j = 0; j < good_segments.size() && !has_non_gap; ++j)
      if (good_segments[j][i] != GenomeAlignment::gap_symbol)
	has_non_gap = true;
    if (!has_non_gap) {
      for (size_t j = 0; j < good_segments.size(); ++j)
	segments[j].append(good_segments[j].begin() + prev_index,
			   good_segments[j].begin() + i);
      prev_index = i + 1;
    }
  }
  for (size_t j = 0; j < good_segments.size(); ++j)
    segments[j].append(good_segments[j].begin() + prev_index,
		       good_segments[j].end());
  
  regions.swap(good_regions);
}



void
extract_regions_dir(const string alndir,
		    const string alnsuffix,
		    const vector<GenomicRegion> &query_regions,
		    const set<string> &species,
		    const bool get_depth,
		    const char *outfile) {
  
  vector<GenomicRegion> file_regions = read_db_dir(alndir, alnsuffix);
  if (VERBOSE) {
    cerr << "valid regions found:" << endl;
    std::copy(file_regions.begin(), file_regions.end(), 
	      std::ostream_iterator<GenomicRegion>(cerr, "\n"));
  }

  vector<vector<GenomicRegion> > query_regions_by_file;
  get_regions_by_filename(file_regions, query_regions, query_regions_by_file);

  ostream* out = (outfile) ? new ofstream(outfile) : &cout;

  for (size_t i = 0; i < query_regions_by_file.size(); ++i) {

    if (!query_regions_by_file[i].empty()) {
      if (VERBOSE)
	cerr << "loading: " << file_regions[i].get_name();
      GenomeAlignment aln(path_join(alndir, file_regions[i].get_name()) + '.' +
			  alnsuffix);
      if (VERBOSE) 
	cerr << ". processing " << query_regions_by_file[i].size() 
	     << " regions." << endl;
      
      for (size_t j = 0; j < query_regions_by_file[i].size(); ++j)
	if (aln.contains(query_regions_by_file[i][j])) {
	  
	  vector<GenomicRegion> regions;
	  vector<string> segments;
	  aln.get_slice(query_regions_by_file[i][j], regions, segments);
	  
	  if (!species.empty())
	    get_species_alignment(species, regions, segments);
	  
	  if (get_depth) {
	    query_regions_by_file[i][j].set_score(regions.size());
	    *out << query_regions_by_file[i][j] << endl;
	  }
	  else
	    *out << AlignmentBlock(regions, segments).maf_block_string() 
		 << endl;
	}
	else if (VERBOSE)
	  cerr << "Region not found:\t" << query_regions_by_file[i][j] << endl;
    }
  }
  if (out != &cout) delete out;
}

void
extract_regions_dir_idx(const string alndir,
			const string alnsuffix,
			const string indexsuffix,
			const vector<GenomicRegion> &query_regions,
			const set<string> &species,
			const bool get_depth,
			const char *outfile) {
  
  vector<GenomicRegion> file_regions = read_db_dir(alndir, alnsuffix);
  if (VERBOSE) {
    cerr << "valid regions found:" << endl;
    std::copy(file_regions.begin(), file_regions.end(), 
	      std::ostream_iterator<GenomicRegion>(cerr, "\n"));
  }

  vector<vector<GenomicRegion> > query_regions_by_file;
  get_regions_by_filename(file_regions, query_regions, query_regions_by_file);

  ostream* out = (outfile) ? new ofstream(outfile) : &cout;

  for (size_t i = 0; i < query_regions_by_file.size(); ++i) {

    if (!query_regions_by_file[i].empty()) {
      if (VERBOSE)
	cerr << "loading: " << file_regions[i].get_name();
      GenomeAlignmentOnDisk aln(path_join(alndir, file_regions[i].get_name()) + '.' +
				indexsuffix,
				path_join(alndir, file_regions[i].get_name()) + '.' +
				alnsuffix);
      if (VERBOSE) 
	cerr << ". processing " << query_regions_by_file[i].size() 
	     << " regions." << endl;
      
      for (size_t j = 0; j < query_regions_by_file[i].size(); ++j)
	if (aln.contains(query_regions_by_file[i][j])) {
	  
	  vector<GenomicRegion> regions;
	  vector<string> segments;
	  aln.get_slice(query_regions_by_file[i][j], regions, segments);
	  
	  if (!species.empty())
	    get_species_alignment(species, regions, segments);
	  
	  if (get_depth) {
	    query_regions_by_file[i][j].set_score(regions.size());
	    *out << query_regions_by_file[i][j] << endl;
	  }
	  else
	    *out << AlignmentBlock(regions, segments).maf_block_string() 
		 << endl;
	}
	else if (VERBOSE)
	  cerr << "Region not found:\t" << query_regions_by_file[i][j] << endl;
    }
  }
  
  if (out != &cout) delete out;
}


void
extract_regions(const string alnfile,
		const vector<GenomicRegion> &query_regions,
		const set<string> &species,
		const bool get_depth,
		const char *outfile) {
  
  GenomeAlignment aln(alnfile);

  ostream* out = (outfile) ? new ofstream(outfile) : &cout;
  
  for (size_t i = 0; i < query_regions.size(); ++i)
    if (aln.contains(query_regions[i])) {

      vector<GenomicRegion> regions;
      vector<string> segments;
      aln.get_slice(query_regions[i], regions, segments);

      if (!species.empty())
	get_species_alignment(species, regions, segments);
      
      if (get_depth) {
	GenomicRegion with_depth(query_regions[i]);
	with_depth.set_score(regions.size());
	*out << with_depth << endl;
      }
      else
	*out << AlignmentBlock(regions, segments).maf_block_string() 
	     << endl;
    }
    else if (VERBOSE)
      cerr << "Region not found:\t" << query_regions[i] << endl;

  if (out != &cout) delete out;
}

void
extract_regions_idx(const string alnfile,
		    const string indexfile,
		    const vector<GenomicRegion> &query_regions,
		    const set<string> &species,
		    const bool get_depth,
		    const char *outfile) {
  
  GenomeAlignmentOnDisk aln(indexfile, alnfile);
  
  ostream* out = (outfile) ? new ofstream(outfile) : &cout;
  
  for (size_t i = 0; i < query_regions.size(); ++i)
    if (aln.contains(query_regions[i])) {
      
      vector<GenomicRegion> regions;
      vector<string> segments;
      aln.get_slice(query_regions[i], regions, segments);
      
      if (!species.empty())
	get_species_alignment(species, regions, segments);
      
      if (get_depth) {
	GenomicRegion with_depth(query_regions[i]);
	with_depth.set_score(regions.size());
	*out << with_depth << endl;
      }
      else
	*out << AlignmentBlock(regions, segments).maf_block_string() 
	     << endl;
    }
    else if (VERBOSE)
      cerr << "Region not found:\t" << query_regions[i] << endl;
  
  if (out != &cout) delete out;
}

// command line parameters
string alnfile;
string alnsuffix;
string outfile;
string species_file;
string regions_file;

string indexfile;
string indexsuffix;
string index_to_generate;
string index_to_test;

int get_depth = false;

int main(int argc, const char **argv) {
  
  try {
    /***************** GET COMMAND LINE ARGUMENTS *******************/
    OptionParser opt_parse(strip_path(argv[0]), " ",
                           " "); 
    opt_parse.add_opt("generate", 'g', "generate index in specified file",
                       false, index_to_generate);
    opt_parse.add_opt("test", 't', "test index in specified file",
                       false, index_to_test);
    opt_parse.add_opt("index", 'i', "index file name", false, indexfile);
    opt_parse.add_opt("index-suff", 'I', "index file name suffix",
                       false, indexsuffix);
    opt_parse.add_opt("aln-suff", 'A', "alignment file name suffix",
                       false, alnsuffix);
    opt_parse.add_opt("output", 'o', "output file (default: stdout;"
                      " no use with -g)", false, outfile);
    opt_parse.add_opt("species", 's', "file of species names",
                       false, species_file);
    opt_parse.add_opt("verbose", 'v',"print more run information",
                       false, VERBOSE);
    opt_parse.add_opt("regions", 'r',"regions to use in testing",
                       false, regions_file);
    opt_parse.add_opt("depth", 'd',"report depth of alignment at regions",
                       false, get_depth);
    /****************************************************************/
    if (index_to_test.c_str())
      return test_index_file(alnfile.c_str(), index_to_test.c_str());
    else if (index_to_generate.c_str()) {
      generate_maf_file_index(alnfile.c_str(), index_to_generate.c_str());
    }
    else {
      set<string> species;
      if (species_file.c_str())
	read_species_names(species_file.c_str(), species);
      
      // Specified regions file means that an alignment file
      // containing the alignments at those regions is to be obtained.
      if (regions_file.c_str()) {
	
	vector<GenomicRegion> regions;
	ReadBEDFile(regions_file.c_str(), regions);
	sort(regions.begin(), regions.end());

	if (alnsuffix.c_str()) {
	  if (indexsuffix.c_str())
	    extract_regions_dir_idx(alnfile.c_str(), alnsuffix.c_str(), indexsuffix.c_str(), regions,
				    species, get_depth, outfile.c_str());
	  else
	    extract_regions_dir(alnfile.c_str(), alnsuffix.c_str(), regions,
				species, get_depth, outfile.c_str());
	}
	else if (indexfile.c_str())
	  extract_regions_idx(alnfile.c_str(), indexfile.c_str(), regions,
			      species, get_depth, outfile.c_str());
	else
	  extract_regions(alnfile.c_str(), regions,
			  species, get_depth, outfile.c_str());
      }
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
