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
#include "Motif.hpp"
#include "SuffixTree.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "PatternFactory.hpp"
#include "Pattern.hpp"

// For getpid() in the random number seed
#include <sys/types.h>
#include <unistd.h>

using std::string;
using std::vector;
using std::greater_equal;
using std::ofstream;
using std::numeric_limits;
using std::greater;
using std::max;
using std::min;
using std::ostream;
using std::cout;
using std::cerr;
using std::endl;
using std::ostream_iterator;


static bool VERBOSE = false;
bool calculate_error_rate = false;
bool calculate_relative_error = false;
float threshold_param = numeric_limits<float>::max();

float spec_param = 0.0;
string label;

const char *err_label          = "ERRORRATE";
const char *errsens_label      = "ERRORRATE_SENSITIVITY";
const char *errspec_label      = "ERRORRATE_SPECIFICITY";
const char *errthresh_label    = "ERRORRATE_THRESHOLD";
const char *errfuncdepth_label = "ERRORRATE_FUNCTIONAL_DEPTH";
const char *err_pval_label     = "ERRORRATE_PVALUE";

const char *relerr_label          = "RELATIVE_ERRORRATE";
const char *relerr_pval_label     = "RELATIVE_ERRORRATE_PVALUE";
const char *relerrsens_label      = "RELATIVE_ERRORRATE_SENSITIVITY";
const char *relerrspec_label      = "RELATIVE_ERRORRATE_SPECIFICITY";
const char *relerrthresh_label    = "RELATIVE_ERRORRATE_THRESHOLD";
const char *relerrfuncdepth_label = "RELATIVE_ERRORRATE_FUNCTIONAL_DEPTH";

const char *binom_label          = "BINOMIAL_PVALUE";
const char *binomsens_label      = "BINOMIAL_SENSITIVITY";
const char *binomspec_label      = "BINOMIAL_SPECIFICITY";
const char *binomthresh_label    = "BINOMIAL_THRESHOLD";
const char *binomfuncdepth_label = "BINOMIAL_FUNCTIONAL_DEPTH";

const char *rank_suffix = "_RANK";

size_t
count_valid_bases(const std::string& s) {
  return count_if(s.begin(), s.end(), &valid_base);
}

size_t
count_valid_bases(const std::vector<std::string>& s) {
  size_t n_valid = 0;
  for (std::vector<std::string>::const_iterator i = s.begin(); i != s.end(); ++i)
    n_valid += count_valid_bases(*i);
  return n_valid;
}

void
get_base_comp(const std::vector<std::string>& sequences, float *base_comp) {
  std::fill(base_comp, base_comp + smithlab::alphabet_size, 0.0);
  float total = 0;
  for (std::vector<std::string>::const_iterator i = sequences.begin();
       i != sequences.end(); ++i)
    for (std::string::const_iterator j = i->begin(); j != i->end(); ++j)
      if (valid_base(*j)) {
        base_comp[base2int(*j)]++;
        total++;
      }
  transform(base_comp, base_comp + smithlab::alphabet_size, base_comp,
            std::bind(std::divides<float>(), std::placeholders::_1, total));
}

double lnchoose(const size_t n, const size_t k) {
  double sum = 0;
  for (size_t i = n; i > k; --i) sum += log(i);
  for (size_t i = n; i > k; --i) sum -= log(n - i + 1);
  // log factorials: it would be good if we had fast functions for
  // these:
  // const size_t n_lnfact = lnfact(n);
  // const size_t m_lnfact = lnfact(m);
  // const size_t nm_lnfact = lnfact(n - m);
  // const size_t val = n_lnfact - m_lnfact - nm_lnfact;
  return sum;
}


#ifdef HAVE_GSL
#include <gsl/gsl_sf_gamma.h>
#endif

double binomial_pdf(const size_t k, const double p,
                        const size_t n) {
  if (k > n) return 0;
#ifdef HAVE_GSL
  const double ln_nCk = gsl_sf_lnchoose(n, k);
#else
  const double ln_nCk = lnchoose(n, k);
#endif
  return exp(ln_nCk + k*log(p) + (n - k)*log(1 - p));
}

// return the sum from 0 to 'k' for the binomial distribution with
// parameters 'n' and 'p'
inline double
binomial(size_t N, size_t k, double p) {
  double z = 0.0;
  if (k < N/2) {
    for (size_t i = 0; i < k; i++)
      z += binomial_pdf(i, p, N);
    return 1 - z;
  }
  else {
    for (size_t i = k; i <= N; i++)
      z += binomial_pdf(i, p, N);
    return z;
  }
}

inline double binomial(float N, float k, double p) {
  return binomial(static_cast<size_t>(N), static_cast<size_t>(k), p);
}

float optimize_binomial(vector<float> &fgmax, vector<float> &bgmax, float p) {
  sort(fgmax.begin(), fgmax.end(), greater<float>());
  sort(bgmax.begin(), bgmax.end(), greater<float>());
  size_t successes = 0, trials = 0;
  float best_threshold = numeric_limits<float>::max(), best_p_value = 1.0;

  for (vector<float>::iterator i = fgmax.begin(), j = bgmax.begin();
       i != fgmax.end(); ++i) {
    ++successes;
    ++trials;

    for (; *j >= *i && j != bgmax.end(); ++j) ++trials;

    float temp_p_value = binomial(trials, successes, p);
    if (temp_p_value <= best_p_value) {
      best_p_value = temp_p_value;
      best_threshold = *i;
    }
  }
  return best_threshold;
}

float optimize_relative_error(vector<float> &fgmax,
                              vector<float> &bgmax, float r) {
  sort(fgmax.begin(), fgmax.end(), greater<float>());
  sort(bgmax.begin(), bgmax.end(), greater<float>());
  float successes = r*bgmax.size(), most_successes = r*bgmax.size();
  float best_threshold = numeric_limits<float>::max();
  float spec = r*bgmax.size();
  float min_spec = spec_param*r*bgmax.size();
  for (vector<float>::iterator i = fgmax.begin(), j = bgmax.begin();
       i != fgmax.end(); ++i) {
    ++successes;

    for (; *j >= *i && j != bgmax.end(); ++j) {
      successes -= r;
      spec -= r;
    }
    if (successes >= most_successes && spec >= min_spec) {
      most_successes = successes;
      best_threshold = *i;
    }
  }
  return best_threshold;
}

float
optimize_relative_error_get_score(vector<float> &fgmax,
                                  vector<float> &bgmax, float r) {
  sort(fgmax.begin(), fgmax.end(), greater<float>());
  sort(bgmax.begin(), bgmax.end(), greater<float>());
  float successes = r*bgmax.size(), most_successes = r*bgmax.size();
  // float best_threshold = numeric_limits<float>::max();
  float spec = r*bgmax.size();
  float min_spec = spec_param*r*bgmax.size();
  for (vector<float>::iterator i = fgmax.begin(), j = bgmax.begin();
       i != fgmax.end(); ++i) {
    ++successes;

    for (; *j >= *i && j != bgmax.end(); ++j) {
      successes -= r;
      spec -= r;
    }
    if (successes >= most_successes && spec > min_spec) {
      most_successes = successes;
      // best_threshold = *i;
    }
  }
  return 1 - most_successes/(r*bgmax.size()+fgmax.size());
}

void
get_max_scores(const vector<ScoringMatrix> &sm,
               const vector<ScoringMatrix> &smrc,
               const vector<string> &sequences,
               const string progress_prefix,
               size_t max_width,
               vector<vector<float> > &maxs) {
  const float denom = sequences.size();
  // get the max scores in each sequence
  for (size_t i = 0; i < sequences.size(); ++i) {
    if (VERBOSE)
      cerr << "\r" << progress_prefix << "\t"
           << static_cast<size_t>(100*(i/denom)) << "%";
    SuffixTree t(sequences[i], max_width);
    for (size_t j = 0; j < sm.size(); ++j)
      maxs[j].push_back(max(t.top_score(sm[j]), t.top_score(smrc[j])));
  }
  if (VERBOSE)
    cerr << "\r" << progress_prefix << "\t100%" << endl;
}



float functional_depth(Motif& motif, float threshold, float bc[]) {
  Matrix m = motif.get_matrix();
  return ScoringMatrix::StormoScoringMatrix(m, bc).functional_depth(threshold);
}


//TODO: need better error handling here
float get_threshold(Motif& motif) {
  if (threshold_param != numeric_limits<float>::max())
    return threshold_param;
  else if (!label.empty()) {
    if (motif.get_attribute(label) == "") {
      cerr << "Motif " << motif.get_accession()
           << " has no attribute with label: " << label << endl;
      exit(EXIT_FAILURE);
    }
    // TODO: test if this is a reasonable threshold (e.g. a number)
    return atof(motif.get_attribute(label).c_str());
  }
  // TODO: check properly for errors
  return numeric_limits<float>::max();
}

void get_counts_above_threshold(vector<float>& fgmax, vector<float>& bgmax,
                                float threshold, float& fgcount,
                                float &bgcount) {
  fgcount = count_if(fgmax.begin(), fgmax.end(),
                     bind2nd(greater_equal<float>(), threshold));
  bgcount = count_if(bgmax.begin(), bgmax.end(),
                     bind2nd(greater_equal<float>(), threshold));
}

void
set_relative_error_attributes(Motif& motif, float fgsize, float bgsize,
                              float fgcount, float bgcount, float fg_bg_ratio,
                              float base_comp[],
                              float thresh = numeric_limits<float>::max()) {
  motif.set_attribute(relerr_label, (fgsize - fgcount + bgcount*fg_bg_ratio)/
                      (fgsize + bgsize*fg_bg_ratio));
  motif.set_attribute(relerrsens_label, fgcount/fgsize);
  motif.set_attribute(relerrspec_label, 1 - bgcount/bgsize);
  motif.set_attribute(relerrthresh_label, thresh);
  float fd = functional_depth(motif, thresh, base_comp);
  motif.set_attribute(relerrfuncdepth_label, fd);
}

void
set_error_rate_attributes(Motif& motif, float fgsize, float bgsize,
                          float fgcount, float bgcount, float base_comp[],
                          float thresh = numeric_limits<float>::max()) {
  motif.set_attribute(err_label,
                      (fgsize - fgcount + bgcount)/(fgsize + bgsize));
  motif.set_attribute(errsens_label, fgcount/fgsize);
  motif.set_attribute(errspec_label, 1 - bgcount/bgsize);
  motif.set_attribute(errthresh_label, thresh);
  float fd = functional_depth(motif, thresh, base_comp);
  motif.set_attribute(errfuncdepth_label, fd);
}

void
set_binomial_attributes(Motif& motif, float fgsize, float bgsize,
                        float fgcount, float bgcount, float p_value, float base_comp[],
                        float thresh = numeric_limits<float>::max()) {
  motif.set_attribute(binom_label, p_value);
  motif.set_attribute(binomsens_label, fgcount/fgsize);
  motif.set_attribute(binomspec_label, 1 - bgcount/bgsize);
  motif.set_attribute(binomthresh_label, thresh);
  float fd = functional_depth(motif, thresh, base_comp);
  motif.set_attribute(binomfuncdepth_label, fd);
}

void
get_relative_errorrate_quantiles_corrected(const vector<vector<float> > &fgmax,
                                           const vector<vector<float> > &bgmax,
                                           const string progress_prefix,
                                           const float fg_bg_ratio,
                                           const size_t n_samples,
                                           vector<float> &quantiles) {
  const size_t n_motifs = fgmax.size();
  const size_t fgsize = fgmax.front().size();
  const size_t bgsize = bgmax.front().size();

  quantiles.clear();

  vector<vector<float> > fgshuff(n_motifs, vector<float>(fgsize));
  vector<vector<float> > bgshuff(n_motifs, vector<float>(bgsize));

  for (size_t i = 0; i < n_samples; ++i) {
    if (VERBOSE)
      cerr << "\r" << progress_prefix << "\t"
           << static_cast<size_t>((100.0*i)/n_samples) << "%";

    vector<vector<float> > tmp_max(fgsize + bgsize, vector<float>(n_motifs));

    for (size_t j = 0; j < n_motifs; ++j) {
      for (size_t k = 0; k < fgsize; ++k)
        tmp_max[k][j] = fgmax[j][k];
      for (size_t k = 0; k < bgsize; ++k)
        tmp_max[fgsize + k][j] = bgmax[j][k];
    }
    random_shuffle(tmp_max.begin(), tmp_max.end());

    for (size_t j = 0; j < n_motifs; ++j) {
      for (size_t k = 0; k < fgsize; ++k)
        fgshuff[j][k] = tmp_max[k][j];
      for (size_t k = 0; k < bgsize; ++k)
        bgshuff[j][k] = tmp_max[fgsize + k][j];
    }

    float best = numeric_limits<float>::max();
    if(calculate_relative_error)
      for (size_t j = 0; j < n_motifs; ++j)
        best = min(best,
                   optimize_relative_error_get_score(fgshuff[j],
                                                     bgshuff[j], fg_bg_ratio));
    if(calculate_error_rate)
      for (size_t j = 0; j < n_motifs; ++j)
        best = min(best,
                   optimize_relative_error_get_score(fgshuff[j],
                                                     bgshuff[j], 1.0));
    quantiles.push_back(best);
  }
  if (VERBOSE)
    cerr << "\r" << progress_prefix << "\t100%" << endl;
  sort(quantiles.begin(), quantiles.end());
}

void
get_relative_errorrate_quantiles(const vector<vector<float> > &fgmax,
                                 const vector<vector<float> > &bgmax,
                                 const string progress_prefix,
                                 const float fg_bg_ratio,
                                 const size_t n_samples,
                                 vector<vector<float> > &quantiles) {
  const size_t n_motifs = fgmax.size();
  const size_t fgsize = fgmax.front().size();
  const size_t bgsize = bgmax.front().size();

  quantiles = vector<vector<float> >(n_motifs);

  vector<vector<float> > fgshuff(n_motifs, vector<float>(fgsize));
  vector<vector<float> > bgshuff(n_motifs, vector<float>(bgsize));

  for (size_t i = 0; i < n_samples; ++i) {
    if (VERBOSE)
      cerr << "\r" << progress_prefix << "\t"
           << static_cast<size_t>((100.0*i)/n_samples) << "%";

    vector<vector<float> > tmp_max(fgsize + bgsize, vector<float>(n_motifs));

    for (size_t j = 0; j < n_motifs; ++j) {
      for (size_t k = 0; k < fgsize; ++k)
        tmp_max[k][j] = fgmax[j][k];
      for (size_t k = 0; k < bgsize; ++k)
        tmp_max[fgsize + k][j] = bgmax[j][k];
    }
    random_shuffle(tmp_max.begin(), tmp_max.end());

    for (size_t j = 0; j < n_motifs; ++j) {
      for (size_t k = 0; k < fgsize; ++k)
        fgshuff[j][k] = tmp_max[k][j];
      for (size_t k = 0; k < bgsize; ++k)
        bgshuff[j][k] = tmp_max[fgsize + k][j];
    }

    if(calculate_relative_error)
      for (size_t j = 0; j < n_motifs; ++j)
        quantiles[j].push_back(optimize_relative_error_get_score(fgshuff[j],
                                                    bgshuff[j], fg_bg_ratio));
    if(calculate_error_rate)
      for (size_t j = 0; j < n_motifs; ++j)
        quantiles[j].push_back(optimize_relative_error_get_score(fgshuff[j],
                                                    bgshuff[j], 1.0));
  }
  if (VERBOSE)
    cerr << "\r" << progress_prefix << "\t100%" << endl;
  for (size_t i = 0; i < quantiles.size(); ++i)
    sort(quantiles[i].begin(), quantiles[i].end());
}

void
set_relerr_pvalue(const vector<float> &quantiles, Motif &motif) {
  const float relerr = atof(motif.get_attribute(relerr_label).c_str());
  vector<float>::const_iterator i = upper_bound(quantiles.begin(),
                                                quantiles.end(), relerr);
  motif.set_attribute(relerr_pval_label,
                      static_cast<double>(i - quantiles.begin())/
                                                  quantiles.size());
}

void
set_err_pvalue(const vector<float> &quantiles, Motif &motif) {
  const float err = atof(motif.get_attribute(err_label).c_str());
  vector<float>::const_iterator i = upper_bound(quantiles.begin(),
                                                quantiles.end(), err);
  motif.set_attribute(err_pval_label,
                      static_cast<double>(i - quantiles.begin())/
                                                  quantiles.size());
}

/*void AssignRanks(vector<Pattern*>& patterns, string key) {
  string key_rank = key + string(rank_suffix);
  for (size_t i = 0, rank = 1; i < patterns.size(); ++i) {
    patterns[i]->set_attribute(key_rank, rank);
    if (i == patterns.size() - 1 ||
        patterns[i]->get_attribute(key) != patterns[i + 1]->get_attribute(key))
      rank = i + 2;
  }
}*/

int main(int argc, const char *argv[]) {
  /* INPUT PARAMETERS */
  string fgfile;
  string bgfile;
  string outfile;
  string motifs_file_name;
  string key;
  string cutoff;

  bool optimize = false;
  bool calculate_binomial = false;
  bool input_threshold_is_functional_depth = false;

  bool reverse_order = false;
  bool numeric_sort_order = false;
  //bool assign_rank = false;

  int relerr_pval_samples = 0;
  int correct_for_multiple_testing = 0;
  static int max_patterns_to_output_param = -1;

  srand(time(0) + getpid());

  try {

  /***************** GET COMMAND LINE ARGUMENTS *******************/
    OptionParser opt_parse(strip_path(argv[0]),
               "uses motifs to classify two sets of sequences ",
               "-f foreground -b background motif-file ");
    opt_parse.add_opt("threshold", 't', "threshold; default is find best",
                       false, threshold_param);
    opt_parse.add_opt("threshold-label", 'l', "label of attribute"
                      " specifying threshold", false, label);
    opt_parse.add_opt("optimize", 'O',"optimize thresholds w.r.t. score"
                      " (e.g. error rate)", false, optimize);
    opt_parse.add_opt("error", 'e', "calculate error rate",
                       false, calculate_error_rate );
    opt_parse.add_opt("func-depth", 'd', "Thresholds are specified as"
                      "functional depth", false,
                       input_threshold_is_functional_depth);
    opt_parse.add_opt("specificity", 'p',"minimum specificity",
                       false, spec_param);
    opt_parse.add_opt("relative-error", 'r',"optimize relative error rate",
                       false, calculate_relative_error);
    opt_parse.add_opt("binomial", 'B', "optimize binomial p-value",
                       false, calculate_binomial);
    opt_parse.add_opt("foreground", 'f', "Foreground", false, fgfile);
    opt_parse.add_opt("background", 'b',"Background", false, bgfile);
    opt_parse.add_opt("output", 'o', "name of output file",
                       false, outfile);
    opt_parse.add_opt("verbose", 'v', "print progress information",
                       false, VERBOSE);
    opt_parse.add_opt("calcpval", 'P', "samples to calculate pval"
                      " (needs -r or -e)", false, relerr_pval_samples);
    opt_parse.add_opt("correct", 'c', "correct for multiple testing"
                      " (needs -P)", false, correct_for_multiple_testing);
    opt_parse.add_opt("key", 'k', "name of attribute on which to sort",
                         false, key);
    opt_parse.add_opt("reverse", 'R', "sort in reverse order",
                         false, reverse_order);
    opt_parse.add_opt("numeric", 'n',"sort as though key has numeric value",
                         false, numeric_sort_order);
    opt_parse.add_opt("cutoff", 'C', "cutoff value for patterns to return",
                       false, cutoff);
    opt_parse.add_opt("ntop", 'm', "only print this many of top outputs",
                       false, max_patterns_to_output_param);

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
    // READ IN THE FOREGROUND SEQUENCES
    FastaFile fg = FastaFile(const_cast<char*>(fgfile.c_str()));
    vector<string> foreground = fg.get_sequences();
    size_t fgsize = foreground.size();
    size_t fglen = count_valid_bases(foreground);

    // READ IN THE BACKGROUND SEQUENCES
    FastaFile bg = FastaFile(const_cast<char*>(bgfile.c_str()));
    vector<string> background = bg.get_sequences();
    size_t bgsize = background.size();
    size_t bglen = count_valid_bases(background);

    // GET THE BASE COMPOSITION
    float base_comp[alphabet_size], base_comp_fg[alphabet_size], base_comp_bg[alphabet_size];
    get_base_comp(foreground, base_comp_fg);
    get_base_comp(background, base_comp_bg);

    for (size_t i = 0; i < alphabet_size; ++i)
      base_comp[i] =
        (fglen*base_comp_fg[i] + bglen*base_comp_bg[i])/(fglen + bglen);

    // CALCULATE THE PROBABILITY TO USE IN THE BINOMIAL (EQUAL TO THE
    // PROPORTION OF FOREGROUND SEQUENCES)
    float null_prob = static_cast<float>(fgsize)/(fgsize + bgsize);
    float fg_bg_ratio = static_cast<float>(fgsize)/bgsize;

    // READ IN THE MOTIFS
    motifs_file_name = input_filenames[0];
    vector<Motif> motifs(Motif::ReadMotifVector(motifs_file_name.c_str()));

    size_t max_motif_width = 0;
    for (vector<Motif>::iterator i = motifs.begin(); i != motifs.end(); ++i)
      max_motif_width = std::max(max_motif_width, i->get_width());

    // MAKE A VECTOR OF SCORING MATRICES (AND REVERSE COMPLEMENT
    // MATRICES IF NECESSARY)
    vector<ScoringMatrix> sm, smrc;
    for (vector<Motif>::iterator i = motifs.begin(); i != motifs.end(); ++i) {
      if (!i->get_matrix().is_count_mat())
        sm.push_back(ScoringMatrix(i->get_matrix().freqmat(), base_comp));
      else
        sm.push_back(ScoringMatrix::StormoScoringMatrix(i->get_matrix(),
                                                        base_comp));
      smrc.push_back(sm.back().revcomp());
    }

    // GET THE ARRAYS OF MAX SCORES
    vector<vector<float> > fgmax(motifs.size()), bgmax(motifs.size());
    get_max_scores(sm, smrc, foreground, "getting foreground scores",
                   max_motif_width, fgmax);
    get_max_scores(sm, smrc, background, "getting background scores",
                   max_motif_width, bgmax);


    // Get the ordered foreground and background scores to keep their
    // order preserved if the relative error p-value is to be
    // calucated
    vector<vector<float> > fgmax_ordered, bgmax_ordered;
    if (relerr_pval_samples > 0 &&
          (calculate_relative_error || calculate_error_rate) &&
          optimize) {
      fgmax_ordered = fgmax;
      bgmax_ordered = bgmax;
    }

    // ITERATE OVER EACH MATRIX AND GET THE REQUESTED QUALITY MEASURES
    for (size_t i = 0; i < motifs.size(); ++i) {
      if (VERBOSE)
        cerr << "\r" << "evaluating motifs\t"
             << static_cast<size_t>((100.0*i)/motifs.size()) << "%";

      float threshold, fgcount, bgcount;

      // IF THE USER DOESN'T WANT OPTIMIZATION:
      if (!optimize) {
        threshold = get_threshold(motifs[i]);
        if (input_threshold_is_functional_depth)
          threshold = sm[i].functional_depth_to_score(threshold);
        get_counts_above_threshold(fgmax[i], bgmax[i], threshold,
                                   fgcount, bgcount);
        if (calculate_error_rate)
          set_error_rate_attributes(motifs[i], fgsize, bgsize,
                                    fgcount, bgcount, base_comp);
        if (calculate_relative_error)
          set_relative_error_attributes(motifs[i], fgsize, bgsize, fgcount,
                                        bgcount, fg_bg_ratio, base_comp,
                                        threshold);
        float p_value;
        if (calculate_binomial) {
          p_value = binomial(fgcount + bgcount, fgcount, null_prob);
          set_binomial_attributes(motifs[i], fgsize, bgsize,
                                  fgcount, bgcount, p_value, base_comp);
        }
      }
      // IF USER SPECIFIED TO OPTIMIZE THE THRESHOLD(S):
      else {
        if (calculate_error_rate) {
           // 1.0 makes it absolute error:
          threshold = optimize_relative_error(fgmax[i], bgmax[i], 1.0);
          get_counts_above_threshold(fgmax[i], bgmax[i],
                                     threshold, fgcount, bgcount);
          set_error_rate_attributes(motifs[i], fgsize, bgsize,
                                    fgcount, bgcount, base_comp, threshold);


        }
        if (calculate_relative_error) {
          threshold = optimize_relative_error(fgmax[i], bgmax[i], fg_bg_ratio);
          get_counts_above_threshold(fgmax[i], bgmax[i], threshold,
                                     fgcount, bgcount);
          set_relative_error_attributes(motifs[i], fgsize, bgsize, fgcount,
                                        bgcount, fg_bg_ratio, base_comp, threshold);
        }
        float p_value;
        if (calculate_binomial) {
          threshold = optimize_binomial(fgmax[i], bgmax[i], null_prob);
          get_counts_above_threshold(fgmax[i], bgmax[i],
                                     threshold, fgcount, bgcount);
          p_value = binomial(fgcount + bgcount, fgcount, null_prob);
          set_binomial_attributes(motifs[i], fgsize, bgsize, fgcount,
                                  bgcount, p_value, base_comp, threshold);
        }
      }
    }
    if (VERBOSE)
      cerr << "\r" << "evaluating motifs\t100%" << endl;

    if (relerr_pval_samples > 0 &&
        (calculate_relative_error || calculate_error_rate)
        && optimize) {
      if (correct_for_multiple_testing) {
        vector<float> quantiles;
        get_relative_errorrate_quantiles_corrected(fgmax_ordered,
                                                   bgmax_ordered,
                                                   "calculating p-values",
                                                   fg_bg_ratio,
                                                   relerr_pval_samples,
                                                   quantiles);
        if (calculate_relative_error)
          for (size_t i = 0; i < motifs.size(); ++i)
            set_relerr_pvalue(quantiles, motifs[i]);
        if(calculate_error_rate)
          for (size_t i = 0; i < motifs.size(); ++i)
            set_err_pvalue(quantiles, motifs[i]);
      }
      else {
        vector<vector<float> > quantiles;
        get_relative_errorrate_quantiles(fgmax_ordered, bgmax_ordered,
                                         "calculating p-values",
                                         fg_bg_ratio, relerr_pval_samples,
                                         quantiles);
        if (calculate_relative_error)
          for (size_t i = 0; i < motifs.size(); ++i)
            set_relerr_pvalue(quantiles[i], motifs[i]);
        if(calculate_error_rate)
          for (size_t i = 0; i < motifs.size(); ++i)
            set_err_pvalue(quantiles[i], motifs[i]);
      }
    }
    if (!key.empty()){
      PatternOrder po(key.c_str(), reverse_order, numeric_sort_order);
      std::stable_sort(motifs.begin(), motifs.end(), po);

      if (!cutoff.empty()) {
        PatternCutoff pc(key.c_str(), cutoff.c_str(), !reverse_order, numeric_sort_order);
        motifs.erase(std::find_if(motifs.begin(), motifs.end(), pc),
                       motifs.end());
      }
      //if (assign_rank) AssignRanks(motifs, key);
    }
    const size_t max_patterns_to_output =
      (max_patterns_to_output_param < 1) ?
      motifs.size() : static_cast<size_t>(max_patterns_to_output_param);
    // OUTPUT THE UPDATED MOTIFS
    ostream* output = ((!outfile.empty()) ? new ofstream(outfile.c_str()) : &cout);
    //copy(motifs.begin(), motifs.end(), ostream_iterator<Motif>(*output, "\n"));
    for (size_t i = 0; i < max_patterns_to_output; ++i)
      *output << (motifs[i]) << endl;
    if (output != &cout) delete output;
  }
  catch (SMITHLABOptionException &e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
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
