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
#include "Motif.hpp"
#include "ScoringMatrix.hpp"
#include "GenomicRegion.hpp"
#include "smithlab_utils.hpp"
#include "GenomeUtils.hpp"
#include "GenomeAlignment.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"


#include <errno.h>
#include <dirent.h>

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

int VERBOSE = 0;

bool
 inline valid_base_id(int c) {
 return (c < static_cast<int>(smithlab::alphabet_size) && c >= 0);
  }

string
motif_and_site_progress_string(const size_t current_motif,
			       const size_t total_motifs,
			       const size_t current_site,
			       const size_t total_sites) {
  const char *motif_progress_part = "processing motifs:";
  const char *site_progress_part = "sites:";
  std::ostringstream s;
  s << "\r" 
    << motif_progress_part << "\t"
    << min((current_motif + 1), total_motifs) << "/" 
    << total_motifs << "\t"
    << site_progress_part << "\t"
    << min((current_site + 1), total_sites) << "/"
    << total_sites;
  return s.str();
}



string
sites_set_and_site_progress_string(const size_t current_set,
				   const size_t total_sets,
				   const size_t current_site,
				   const size_t total_sites) {
  const char *set_progress_part = "processing alignment file";
  std::ostringstream s;
  s << "\r" 
    << set_progress_part << "\t"
    << min((current_set + 1), total_sets) << "/" 
    << total_sets << "\t("
    << static_cast<size_t>((100.0*min((current_site + 1), 
				      total_sites))/total_sites)
    << "%)";
  return s.str();
}



void
read_species_names(const string filename, set<string> &species_names) {
  static const size_t buffer_size = 1000; // Magic
  
  // open and check the file
  std::ifstream in(filename.c_str());
  if (!in) 
    throw CREADException("cannot open species file: " + filename);
  
  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw "Line too long in file: " + filename;
    species_names.insert(species_names.end(), buffer);
    in.peek();
  }
  in.close();
}

GenomicRegion
site2region(const MotifSite& site) {
  string chrom;
  size_t start = 0, end = 0;
  parse_region_name(site.get_seq_name(), chrom, start, end);
  start += site.get_start();
  return GenomicRegion(chrom, start, start + site.get_length());
}


// string 
// path_join(const string& a, const string& b) {
//   return a + "/" + b;
// }

string 
append_suffix(const string& a, const string& b) {
  return a + "." + b;
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



// bool
// is_valid_filename(const string name, const string& filename_suffix) {
//   string chrom;
//   size_t start = 0, end = 0;
//   string basename(name.substr(0, name.find_first_of(".")));
//   parse_region_name(basename, chrom, start, end);
//   if (name.find(':') != string::npos)
//     chrom += ":" + cread::toa(start) + "-" + cread::toa(end);
//   return (name == chrom + '.' + filename_suffix);
// }



vector<GenomicRegion> 
read_db_dir(const string& dirname, string filename_suffix) {
  DIR *dir;
  if (!(dir = opendir(dirname.c_str())))
    throw CREADException("could not open directory: " + dirname);
  
  vector<string> filenames;
  errno = 0;
  struct dirent *ent;
  while ((ent = readdir(dir))) {
    if (is_valid_filename(ent->d_name, filename_suffix))
      filenames.push_back(ent->d_name);
    errno = 0;
  }

  // check for some errors
  if (errno)
    throw CREADException("error reading db directory: " + dirname);
  if (filenames.empty())
    throw CREADException("no valid files found in dir: " + dirname);
  
  // strip the suffix from the filenames
  transform(filenames.begin(), filenames.end(),
	    filenames.begin(), std::bind2nd(std::ptr_fun(&strip_suffix),
					    filename_suffix));
  
  // turn the filenames into the regions they represent
  vector<GenomicRegion> file_regions;
  transform(filenames.begin(), filenames.end(),
	    back_inserter(file_regions), std::ptr_fun(&filename2region));
  sort(file_regions.begin(), file_regions.end());
  return file_regions;
}



void 
get_sites_by_filename(const vector<GenomicRegion>& file_regions,
		      const vector<Motif>& motifs,
		      vector<vector<MotifSite> >& sites_by_file,
		      vector<vector<size_t> >& motif_id) {
  
  // extract the sites from the input motifs
  transform(motifs.begin(), motifs.end(), back_inserter(sites_by_file),
	    std::mem_fun_ref(&Motif::get_sites));
  
  sites_by_file = vector<vector<MotifSite> >(file_regions.size());
  motif_id = vector<vector<size_t> >(file_regions.size());
  
  for (size_t i = 0; i < motifs.size(); ++i) {
    vector<MotifSite> sites(motifs[i].get_sites());
    for (size_t j = 0; j < sites.size(); ++j) {
      GenomicRegion region(site2region(sites[j]));
      const size_t file_id = 
	find_closest(file_regions, region) - file_regions.begin();
      sites_by_file[file_id].push_back(sites[j]);
      motif_id[file_id].push_back(i);
    }
  }
}



void
cut_out_non_ref(vector<string> &seqs) {
  vector<bool> keep(seqs.front().length());
  for (size_t i = 0; i < seqs.front().length(); ++i)
    keep[i] = (seqs.front()[i] != '-');
  vector<string> work_seqs(seqs.size());
  for (size_t i = 0; i < seqs.size(); ++i)
    for (size_t j = 0; j < seqs[i].length(); ++j)
      if (keep[j])
	work_seqs[i] += seqs[i][j];
  seqs.swap(work_seqs);
}



void
prepare_alignment_data(const set<string> &species,
		       vector<GenomicRegion> &regions,
		       vector<string> &seqs) {

  vector<GenomicRegion> good_regions;
  vector<string> good_seqs;
  for (size_t i = 0; i < seqs.size(); ++i) {
    const string species_name = regions[i].get_name();
    if ((species.empty() || species.find(species_name) != species.end()) &&
	seqs[i].find_first_of("-") == string::npos) {
      good_regions.push_back(regions[i]);
      good_seqs.push_back(seqs[i]);
    }
  }
  regions.swap(good_regions);
  seqs.swap(good_seqs);
}



float 
match(const ScoringMatrix &sm, const vector<string>& seqs) {
  const size_t matwidth = sm.get_width();
  float score = 0;
  for (size_t s = 0; s < seqs.size(); ++s) {
    size_t j = 0;
    const size_t seqlen = seqs[s].length();
    for (size_t i = 0; i < seqlen && j < matwidth; ++i) {
      const size_t base_id = base2int(seqs[s][i]);
      if (valid_base_id(base_id)) {
	score += sm[j][base_id];
	++j;
      }
    }
  }
  return score/seqs.size();
}



float 
count_matches(const ScoringMatrix &sm, const vector<string>& seqs,
	      const float cutoff) {
  float counts = 0;
  for (size_t s = 0; s < seqs.size(); ++s) {
    size_t j = 0;
    float score = 0;
    for (size_t i = 0; i < seqs[s].length() && j < sm.get_width(); ++i) {
      const size_t base_id = base2int(seqs[s][i]);
      if (valid_base_id(base_id)) {
	score += sm[j][base_id];
	++j;
      }
    }
    counts += (j == sm.get_width() && score > cutoff);
  }
  return counts;
}



void
verify_site_aln_consistency(const MotifSite &site, string ref_site) {
  transform(ref_site.begin(), ref_site.end(), 
	    ref_site.begin(), std::ptr_fun(&::toupper));
  if (!site.posstrand())
    ref_site = revcomp(ref_site);
  if (ref_site != site.get_site()) {
    const string message = "inconsistent alignment data for site:\n" + 
      site.tostring() + " (" + ref_site + ")\n"
      "likely causes: inconsistent genome versions\n"
      "               incorrect site locations";
    throw CREADException(message);
  }
}



float
score_site(const set<string> &species, 
	   const StadenPValue &spv, 
	   const MotifSite &site,
	   const ScoringMatrix &sm, 
	   const ScoringMatrix &smrc,
	   const float cutoff,
	   vector<GenomicRegion> &regions, 
	   vector<string> &slice) {
  
  cut_out_non_ref(slice);
  prepare_alignment_data(species, regions, slice);

  verify_site_aln_consistency(site, slice.front());
  
  if (site.posstrand())
    return (cutoff != -numeric_limits<float>::max()) ?
      count_matches(sm, slice, cutoff) :
      spv.get_pvalue(sm, match(sm, slice));
  else
    return (cutoff != -numeric_limits<float>::max()) ?
      count_matches(smrc, slice, cutoff) :
      spv.get_pvalue(smrc, match(smrc, slice));
}



float
score_site_aln(const set<string> &species, 
	       const StadenPValue &spv,
	       const MotifSite &site, 
	       const GenomeAlignment &aln,
	       const ScoringMatrix &sm, 
	       const ScoringMatrix &smrc,
	       const float cutoff) {
  GenomicRegion region(site2region(site));
  if (aln.contains(region)) {
    vector<string> slice;
    vector<GenomicRegion> regions;
    aln.get_slice(region, regions, slice);
    return score_site(species, spv, site, sm, smrc, cutoff, regions, slice);
  }
  if (VERBOSE)
    cerr << "\n" << region << "\t(not found)" << endl;
  return (cutoff == -numeric_limits<float>::max()) ? 1.0 : 0.0;
}

float
score_site_aln_on_disk(const set<string> &species, 
		       const StadenPValue &spv,
		       const MotifSite &site, 
		       const GenomeAlignmentOnDisk &aln,
		       const ScoringMatrix &sm, 
		       const ScoringMatrix &smrc,
		       const float cutoff) {
  GenomicRegion region(site2region(site));
  if (aln.contains(region)) {
    vector<string> slice;
    vector<GenomicRegion> regions;
    aln.get_slice(region, regions, slice);
    return score_site(species, spv, site, sm, smrc, cutoff, regions, slice);
  }
  if (VERBOSE)
    cerr << "\n" << region << "\t(not found)" << endl;
  return (cutoff == -numeric_limits<float>::max()) ? 1.0 : 0.0;
}

  

void
process_alnfile(const string alnfile,
		const vector<float> &base_comp,
		const float given_cutoff,
		const set<string> &species,
		vector<Motif>& motifs) {
  /*
    This constructor reads the data from the file
  */
  if (VERBOSE)
    cerr << "loading: " << alnfile;
  GenomeAlignment aln(alnfile);
  if (VERBOSE) 
    cerr << ". processing" << endl;
  
  StadenPValue spv(base_comp, 1000);
  
  /*
    Process each motif by first copying its sites, then erasing
    them. Then iterate over all the sites and obtain their new scores
    (if one was able to be obtained), adding each back to the motif
    with the updated score.
  */
  
  for (size_t i = 0; i < motifs.size(); ++i) {
    
    vector<MotifSite> sites(motifs[i].get_sites());
    motifs[i].clear_sites();
    
    const ScoringMatrix sm =
      ScoringMatrix::StormoScoringMatrix(motifs[i].get_matrix(), base_comp);
    const ScoringMatrix smrc(sm.revcomp());
    
    const float cutoff = (given_cutoff != -numeric_limits<float>::max()) ?
      spv.get_score(sm, given_cutoff) : given_cutoff;
    
    for (size_t j = 0; j < sites.size(); ++j) {
      if (VERBOSE) 
	cerr << motif_and_site_progress_string(i, motifs.size(),
					       j, sites.size());
      sites[j].set_score(score_site_aln(species, spv, sites[j], aln,
					sm, smrc, cutoff));
      motifs[i].add_site(sites[j]);
    }
  }
  if (VERBOSE) 
    cerr << endl;
}



void
process_alnfile(const string alnfile,
		const string idfile,
		const vector<float> &base_comp,
		const float given_cutoff,
		const set<string> &species,
		vector<Motif>& motifs) {
  /*
    This constructor reads the data from the file
  */
  if (VERBOSE)
    cerr << "loading: " << alnfile;
  GenomeAlignmentOnDisk aln(idfile, alnfile);
  if (VERBOSE) 
    cerr << ". processing" << endl;
  
  StadenPValue spv(base_comp, 1000);
  
  /*
    Process each motif by first copying its sites, then erasing
    them. Then iterate over all the sites and obtain their new scores
    (if one was able to be obtained), adding each back to the motif
    with the updated score.
  */
  
  for (size_t i = 0; i < motifs.size(); ++i) {
    
    vector<MotifSite> sites(motifs[i].get_sites());
    motifs[i].clear_sites();
    
    const ScoringMatrix sm =
      ScoringMatrix::StormoScoringMatrix(motifs[i].get_matrix(), base_comp);
    const ScoringMatrix smrc(sm.revcomp());
    
    const float cutoff = (given_cutoff != -numeric_limits<float>::max()) ?
      spv.get_score(sm, given_cutoff) : given_cutoff;
    
    for (size_t j = 0; j < sites.size(); ++j) {
      if (VERBOSE) 
	cerr << motif_and_site_progress_string(i, motifs.size(), j, sites.size());
      sites[j].set_score(score_site_aln_on_disk(species, spv, sites[j], aln,
						sm, smrc, cutoff));
      motifs[i].add_site(sites[j]);
    }
  }
  if (VERBOSE) 
    cerr << endl;
}



void
process_alndir(const string alndir, 
	       const string alnsuff,
	       const vector<float> & base_comp, 
	       const float cutoff,
	       const set<string> &species, 
	       vector<Motif>& motifs) {
  
  vector<GenomicRegion> file_regions = read_db_dir(alndir, alnsuff);
  if (VERBOSE) {
    cerr << "valid regions found:" << endl;
    std::copy(file_regions.begin(), file_regions.end(), 
	      std::ostream_iterator<GenomicRegion>(cerr, "\n"));
  }
  
  StadenPValue spv(base_comp, 1000);
  
  vector<vector<size_t> > motif_id;
  vector<vector<MotifSite> > sites;
  get_sites_by_filename(file_regions, motifs, sites, motif_id);
  
  vector<ScoringMatrix> sm;
  vector<ScoringMatrix> smrc;
  vector<float> cutoffs;
  for (size_t i = 0; i < motifs.size(); ++i) {
    motifs[i].clear_sites();
    sm.push_back(ScoringMatrix::StormoScoringMatrix(motifs[i].get_matrix(),
						    base_comp));
    smrc.push_back(sm.back().revcomp());
    cutoffs.push_back((cutoff != -numeric_limits<float>::max())
		      ? spv.get_score(sm[i], cutoff) : cutoff);
    if (VERBOSE && cutoff != -numeric_limits<float>::max())
      cerr << motifs[i].get_accession() << " cutoff " << cutoffs.back() << endl;
  }
  
  // iterate over the sites from the input
  for (size_t i = 0; i < sites.size(); ++i)
    if (!sites[i].empty()) {
      GenomeAlignment aln(path_join(alndir, 
				    append_suffix(file_regions[i].get_name(),
						  alnsuff)));
      if (VERBOSE)
	cerr << "loaded " << file_regions[i].get_name()
	     << " (" << sites[i].size() << " sites)" << endl;
      for (size_t j = 0; j < sites[i].size(); ++j) {
	if (VERBOSE)
	  cerr << sites_set_and_site_progress_string(i, sites.size(),
						     j, sites[i].size());
	
	sites[i][j].set_score(score_site_aln(species, spv, 
					     sites[i][j], aln,
					     sm[motif_id[i][j]], 
					     smrc[motif_id[i][j]],
					     cutoffs[motif_id[i][j]]));
	motifs[motif_id[i][j]].add_site(sites[i][j]);
      }
    }
  if (VERBOSE)
    cerr << endl;
}



void
process_alndir(const string alndir, 
	       const string alnsuff, 
	       const string iddir,
	       const string idsuff,
	       const vector<float> &base_comp, 
	       const float cutoff,
	       const set<string> &species, 
	       vector<Motif>& motifs) {
  
  // Read the relevant files from the specified directory.  The
  // directory essentially acts as the database, and each file
  // contains the part associated with a particular chromosome or
  // chunk of a chromosome.
  vector<GenomicRegion> file_regions = read_db_dir(alndir, alnsuff);
  if (VERBOSE) {
    cerr << "valid regions found:" << endl;
    std::copy(file_regions.begin(), file_regions.end(), 
	      std::ostream_iterator<GenomicRegion>(cerr, "\n"));
  }

  // Create the StadenPValue object to convert between motif match
  // scores and their associated p-values.
  StadenPValue spv(base_comp, 1000);
  
  // Take the sites associated with each motif, and partition them
  // according to the alignment files in which they are located.
  vector<vector<size_t> > motif_id;
  vector<vector<MotifSite> > sites;
  get_sites_by_filename(file_regions, motifs, sites, motif_id);
  
  vector<ScoringMatrix> sm;
  vector<ScoringMatrix> smrc;
  vector<float> cutoffs;
  for (size_t i = 0; i < motifs.size(); ++i) {
    motifs[i].clear_sites();
    sm.push_back(ScoringMatrix::StormoScoringMatrix(motifs[i].get_matrix(),
						    base_comp));
    smrc.push_back(sm.back().revcomp());
    cutoffs.push_back((cutoff != -numeric_limits<float>::max())
		      ? spv.get_score(sm[i], cutoff) : cutoff);
    if (VERBOSE && cutoff != -numeric_limits<float>::max())
      cerr << motifs[i].get_accession() << " cutoff " << cutoffs.back() << endl;
  }
  
  // Top level iteration is over sets of sites that correspond to each
  // alignment file.
  for (size_t i = 0; i < sites.size(); ++i)
    if (!sites[i].empty()) {
      if (VERBOSE)
	cerr << "loading: " << file_regions[i].get_name();
      GenomeAlignmentOnDisk aln(path_join(iddir,
					  append_suffix(file_regions[i].get_name(),
							idsuff)),
				path_join(alndir,
					  append_suffix(file_regions[i].get_name(),
							alnsuff)));
      if (VERBOSE)
	cerr << ". (" << sites[i].size() << " sites)" << endl;
      for (size_t j = 0; j < sites[i].size(); ++j) {
	if (VERBOSE)
	  cerr << sites_set_and_site_progress_string(i, sites.size(),
						     j, sites[i].size());
	sites[i][j].set_score(score_site_aln_on_disk(species, spv, 
						     sites[i][j], aln,
						     sm[motif_id[i][j]],
						     smrc[motif_id[i][j]],
						     cutoffs[motif_id[i][j]]));
	motifs[motif_id[i][j]].add_site(sites[i][j]);
      }
    }
  if (VERBOSE)
    cerr << endl;
}



vector<float>
parse_base_comp_str(string bc_str) {
  const float base_comp_tolerance = 1e-05;
  
  vector<float> base_comp(alphabet_size);
  
  vector<string> parts(cread::split(bc_str, ","));
  if (parts.size() != alphabet_size)
    throw CREADException("incorrect base composition: " + 
			 cread::toa(bc_str));
  float total = 0;
  for (size_t i = 0; i < alphabet_size; ++i) {
    base_comp[i] = atof(parts[i].c_str());
    total += base_comp[i];
  }
  if (std::abs(static_cast<float>(1) - total) > base_comp_tolerance)
    throw CREADException("base composition does not sum to 1.0");

  return base_comp;
}



int
main(int argc, const char **argv) {

  std::cout << "in main" << std::endl;

  try {
    
    // filenames
    string outfile;
    string motifsfile; //this is not being assigned
    string alnfile;
    string alnsuff;
    string idfile;
    string idsuff;
    string species_file;
    


    // other command line parameters
    string base_comp_str;     // string of base composition
    float cutoff = -std::numeric_limits<float>::max();
    
 //    static struct poptOption optionsTable[] = {
 //      { "output", 'o', POPT_ARG_STRING, &outfile, 0, 
	// "output file (default: stdout)" },

 //      { "alignment", 'a', POPT_ARG_STRING, &alnfile, 0, 
	// "alignment file or directory name" },

 //      { "alnsuff", 'A', POPT_ARG_STRING, &alnsuff, 0, 
	// "alignment file suffix (assumes alignment directory specified)" },

 //      { "index", 'i', POPT_ARG_STRING, &idfile, 0,
	// "index file or directory name" },

 //      { "idsuff", 'I', POPT_ARG_STRING, &idsuff, 0,
	// "suffix of index files (assumes index directory specified)" },

 //      { "species", 's', POPT_ARG_STRING, &species_file, 0,
	// "file containing list of species" },

 //      { "base-comp", 'C', POPT_ARG_STRING, &base_comp_str, 0,
	// "comma separated base freqs (order:A,C,G,T)" },

 //      { "cutoff", 'c', POPT_ARG_FLOAT, &cutoff, 0,
	// "score cutoff for sites" },

 //      { "verbose", 'v', POPT_ARG_NONE, &VERBOSE, 0,
	// "print more run information" },
 //      POPT_TABLEEND
 //    };
    
    /***************** GET COMMAND LINE ARGUMENTS *******************/

    OptionParser opt_parse(strip_path(argv[0]),
                "evaluates storm output");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                       false, outfile);
    opt_parse.add_opt("alignment", 'a', "alignment file or directory name",
                       false, alnfile);
    opt_parse.add_opt("alnsuff", 'A', "alignment file suffix (assumes alignment directory specified)",
                       false, alnsuff);
    opt_parse.add_opt("index", 'i', "index file or directory name",
                       false, idfile);
    opt_parse.add_opt("idsuff", 'I', "suffix of index files (assumes index directory specified)",
                       false, idsuff);
    opt_parse.add_opt("species", 's', "file containing list of species",
                       false, species_file);
    opt_parse.add_opt("base-comp", 'C', "comma separated base freqs (order:A,C,G,T)",
                       false, base_comp_str);
    opt_parse.add_opt("cutoff", 'c', "score cutoff for sites",
                       false, cutoff);
    opt_parse.add_opt("verbose", 'v', "print more run information",
                       false, VERBOSE);

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

    motifsfile = argv[argc - 1];   

    std::cout << "outfile = " << outfile << "  motifsfile = " << motifsfile
    << "  alnfile = " << alnfile << "  argc = " << argc << std::endl;




    // poptContext opCo = poptGetContext("multistorm", argc, argv, optionsTable, 0);
    // poptSetOtherOptionHelp(opCo, "[OPTIONS] motifs-file");
    // if (argc < 2) {
    //   poptPrintHelp(opCo, stderr, 0);
    //   return EXIT_SUCCESS;
    // }
    // char c;
    // if ((c = poptGetNextOpt(opCo)) < -1) {
    //   cerr << "multistorm: bad argument " 
	   // << poptBadOption(opCo, POPT_BADOPTION_NOALIAS) << ": "
	   // << poptStrerror(c) << endl;
    //   return EXIT_FAILURE;
    // }
    // if (!poptPeekArg(opCo)) {
    //   poptPrintHelp(opCo, stderr, 0);
    //   return EXIT_FAILURE;
    // }
    // else motifsfile = poptGetArg(opCo);
    // if (poptPeekArg(opCo)) {
    //   cerr << "multistorm: leftover argument " << poptGetArg(opCo) << endl;
    //   return EXIT_FAILURE;
    // }
    // if (!alnfile) {
    //   cerr << "ERROR: must specify alignment file or directory" << endl;
    //   poptPrintHelp(opCo, stderr, 0);
    //   return EXIT_FAILURE;
    // }
    // if (idfile && alnsuff && !idsuff) {
    //   cerr << "ERROR: must specify index file suffix if using "
	   // << "index with alignment directory" << endl;
    //   poptPrintHelp(opCo, stderr, 0);
    //   return EXIT_FAILURE;
    // }
    // if (idsuff && !idfile && !alnsuff) {
    //   cerr << "ERROR: only specify index file suffix if using "
	   // << "index with alignment directory" << endl;
    //   poptPrintHelp(opCo, stderr, 0);
    //   return EXIT_FAILURE;
    // }
    // poptFreeContext(opCo);
    /****************************************************************/
    
    vector<Motif> motifs(Motif::ReadMotifVector(motifsfile));
    
    const vector<float> base_comp = (!base_comp_str.empty()) ? 
      parse_base_comp_str(base_comp_str) : 
      vector<float>(alphabet_size, 1.0/alphabet_size);

    // Tree stuff
    set<string> species;
    if (!species_file.empty())
      read_species_names(species_file, species);
    
    // processing alignment directory
    if (!alnsuff.empty()) {
      
      if (!idsuff.empty()) { 
	if (!idfile.empty())
	  process_alndir(alnfile, alnsuff, idfile, idsuff, 
			 base_comp, cutoff, species, motifs);
	else
	  process_alndir(alnfile, alnsuff, alnfile, idsuff,
			 base_comp, cutoff, species, motifs);
      }
      else
	process_alndir(alnfile, alnsuff, base_comp, cutoff, species, motifs);
    }
    
    // processing SINGLE alignment file
    else {
      
      if (!idfile.empty())
	process_alnfile(alnfile, idfile, base_comp, cutoff, species, motifs);
      else 
	process_alnfile(alnfile, base_comp, cutoff, species, motifs);
    
    }
    
    std::ostream* out = (!outfile.empty()) ? new std::ofstream(outfile) : &cout;
    copy(motifs.begin(), motifs.end(),
	 std::ostream_iterator<Motif>(*out, "\n"));
    if (out != &cout) delete out;
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
