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
#include "Observation.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::pair;
using std::vector;
using std::ostream;
using std::ostream_iterator;
using std::ofstream;
using std::ostringstream;
using std::cout;
using std::cerr;
using std::endl;

// FILE NAMES


typedef Observation<float, bool> FBObs;

struct feature_info {
  feature_info(string &n, bool fas, float sv, float se, float sp, float e) :
    name(n), fg_abv_split(fas), split_val(sv), sens(se), spec(sp), err(e) {}
  bool predict(FBObs &, size_t) const;
  void addstats(vector<FBObs> &, size_t);
  string name;
  bool fg_abv_split;
  float split_val, sens, spec, err;
};

void
feature_info::addstats(vector<FBObs> &obs, size_t index) {
  float pweight = 0, nweight = 0;
  sens = spec = err = 0;
  for (size_t i = 0; i < obs.size(); ++i) {
    float weight = obs[i].get_weight();
    if (obs[i].get_outcome()) {
      pweight += weight;
      sens += (predict(obs[i], index)) ? weight : 0;
    }
    else {
      nweight += weight;
      spec += (!predict(obs[i], index)) ? weight : 0;
    }
  }
  err = (pweight - sens + nweight - spec)/(pweight + nweight);
  sens = (pweight == 0) ? 0 : sens/pweight;
  spec = (nweight == 0) ? 0 : spec/nweight;
}

bool
feature_info::predict(FBObs &obs, size_t index) const {
  return (fg_abv_split && obs[index] > split_val) ||
    (!fg_abv_split && obs[index] < split_val); 
}

std::ostream& operator<<(std::ostream &s, const feature_info &fi) {
  return s << fi.name << "\t" << fi.fg_abv_split << "\t"<< fi.split_val 
	   << "\t" << fi.sens << "\t" << fi.spec << "\t" << fi.err;
}

class InfoOrder : public std::binary_function<const feature_info&, 
		  const feature_info&, bool> {
public:
  explicit InfoOrder() {}
  bool operator()(const feature_info &i1, const feature_info &i2) const {
    return i1.err < i2.err;
  }
};

template <class U, class T> T get2nd(pair<U, T> &p) {
  return p.second;
}

feature_info feateval(vector<FBObs> &obs, 
		      string &name, size_t index) {
  float tot_weight = accumulate(obs.begin(), obs.end(), 0.0, FBObs::accum_weight);
  float pos_weight = accumulate(obs.begin(), obs.end(), 0.0, FBObs::accum_pos_weight);
  
  vector<pair<float, FBObs*> > helper(obs.size());
  for (size_t i = 0; i < obs.size(); ++i) {
    helper[i].first = obs[i][index];
    helper[i].second = &obs[i];
  }
  sort(helper.begin(), helper.end());
  vector<FBObs*> obs_ptrs;
  transform(helper.begin(), helper.end(), back_inserter(obs_ptrs),
	    std::ptr_fun(get2nd<float, FBObs*>));
  
  float best_err = std::min(pos_weight, tot_weight - pos_weight);
  bool abv_splt = true;
  float best_split = 0;
  
  float temp_pos = 0.0, temp_tot = 0.0;
  bool temp_above = true;
  for (size_t i = 0; i < obs.size() - 1; ++i) {
    temp_tot += obs_ptrs[i]->get_weight();
    if (obs_ptrs[i]->get_outcome()) 
      temp_pos += obs_ptrs[i]->get_weight();
    if ((*obs_ptrs[i])[index] != (*obs_ptrs[i + 1])[index]) {
      temp_above = (temp_pos < (temp_tot - temp_pos));
      float err = std::min(temp_pos + (tot_weight - pos_weight) - 
			   (temp_tot - temp_pos),
			   (temp_tot - temp_pos) + (pos_weight - temp_pos));
      if (err < best_err) {
	best_err = err;
	best_split = ((*obs_ptrs[i])[index] + (*obs_ptrs[i + 1])[index])/2;
	abv_splt = temp_above;
      }
    }
  }
  feature_info fi(name, abv_splt, best_split, 0, 0, 0);
  fi.addstats(obs, index);
  return fi;
}

void get_train_test(vector<FBObs> &obs, size_t fold, size_t part,
		    vector<FBObs> &training, vector<FBObs> &testing) {
  training.clear();
  testing.clear();
  float subset_size = obs.size()/static_cast<float>(fold);
  size_t s = static_cast<size_t>(part*subset_size);
  size_t e = static_cast<size_t>((part + 1)*subset_size);
  copy(obs.begin(), obs.begin() + s, back_inserter(training));
  copy(obs.begin() + s, obs.begin() + e, back_inserter(testing));
  copy(obs.begin() + e, obs.end(), back_inserter(training));
}

int main (int argc, const char **argv) {
  
  string datafile;
  string outfile;

  // int print_estimates = false;
  // size_t print_performance = false;
  bool cross_validate = false;
  bool balance = false;

  try {
  /**************** GET COMMAND LINE ARGUMENTS **********************/

  OptionParser opt_parse(strip_path(argv[0])," ", " ");
  opt_parse.add_opt("balance", 'b', 
         "balance contribution of positive and negative data", false, balance);
  // TODO: implement these
  //opt_parse.add_opt("performance", 'p',  "output training error", 
                    // false, &print_performance);
  //opt_parse.add_opt("estimates", 'e',
       //"print estimated outcome for each observation", false, &print_estimates);

  opt_parse.add_opt("cross-validate", 'c', "perform cross validation", 
                     false, cross_validate);
  opt_parse.add_opt("output", 'o', "name of output file", 
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
  const vector<string> files(leftover_args);//need to change

  /**********************************************************************/
    // READ IN THE DATA
    datafile = files[0];
    vector<FBObs> observations;
    vector<string> feature_names;
    FBObs::read_observations(datafile.c_str(), feature_names, observations);
    if (balance) {
      float total = accumulate(observations.begin(), observations.end(), 0.0,
			       FBObs::accum_weight);
      float positive = accumulate(observations.begin(), observations.end(), 0.0,
				  FBObs::accum_pos_weight);
      for (size_t i = 0; i < observations.size(); ++i)
	if (observations[i].get_outcome()) 
	  observations[i].set_weight(1);
	else observations[i].set_weight(positive/(total - positive));
    }
    
    vector<feature_info> fi;
    if (cross_validate) {
      Observation<float, bool>::shuffle(observations);
      vector<FBObs> training, testing;
      get_train_test(observations, cross_validate, 0, training, testing);
      float test_weight = accumulate(testing.begin(), testing.end(), 0.0,
				     FBObs::accum_weight);
      float test_pos = accumulate(testing.begin(), testing.end(), 0.0,
					 FBObs::accum_pos_weight);
      float test_neg = test_weight - test_pos;
      float total_weight = test_weight;
      float total_pos = test_pos;
      float total_neg = test_neg;
      for (size_t j = 0; j < feature_names.size(); ++j) {
	fi.push_back(feateval(training, feature_names[j], j));
	fi.back().addstats(testing, j);
      }
      for (size_t i = 1; i < cross_validate; ++i) {
	get_train_test(observations, cross_validate, i, training, testing);
	test_weight = accumulate(testing.begin(), testing.end(), 0.0,
				 FBObs::accum_weight);
	test_pos = accumulate(testing.begin(), testing.end(), 0.0,
				     FBObs::accum_pos_weight);
	test_neg = test_weight - test_pos;
	for (size_t j = 0; j < feature_names.size(); ++j) {
	  feature_info temp_fi = feateval(training, feature_names[j], j);
	  temp_fi.addstats(testing, j);
	  fi[j].err = (total_weight*fi[j].err + test_weight*temp_fi.err)/
	    (total_weight + test_weight);
	  fi[j].sens = (total_pos*fi[j].sens + test_pos*temp_fi.sens)/
	    ((total_pos + test_pos) ? (total_pos + test_pos) : 1);
	  fi[j].spec = (total_neg*fi[j].spec + test_neg*temp_fi.spec)/
	    ((total_neg + test_neg) ? (total_neg + test_neg) : 1);
	}
	total_weight += test_weight;
	total_pos += test_pos;
	total_neg += test_neg;
      }
    }
    else 
      for (size_t i = 0; i < feature_names.size(); ++i)
	fi.push_back(feateval(observations, feature_names[i], i));
    
    sort(fi.begin(), fi.end(), InfoOrder());
    ostream* output = (!outfile.empty()) ? new ofstream(outfile.c_str()) : &cout;
    copy(fi.begin(), fi.end(), ostream_iterator<feature_info>(*output, "\n"));
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
