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

#ifndef OBSERVATION_HPP
#define OBSERVATION_HPP

#include "cread.hpp"
#include "Feature.hpp"

#include <sys/types.h>
#include <unistd.h>
#include <cassert>

class InconsistentObservations : public CREADException {
public:
  InconsistentObservations() {
    message = "Inconsistent sizes for observations provided as input";
  }
};

class ObservationFileException : public CREADException{
public:
  ObservationFileException(std::string m = std::string(), int l = -1) : 
    CREADException(m), line_number(l) {}
  int GetLineNumber() const {return line_number;}
private:
  int line_number;
};


/**
 * Observation class exists to represent the values of all predictors
 * for an observation, the outcome of the observation, the name of the
 * observation, and the weight of the observation.
 */
template <class T = float, class U = bool> class Observation {

  const static size_t obs_line_size = 10000000;
  
  std::string name;
  std::vector<T> predictors;
  U outcome;
  float weight;

public:
  /** Constructor with parameters to initialize each instance variable (no defaults) */
  Observation(std::string nm, std::vector<T>& p, U o, float w = 1) : 
    name(nm), predictors(p), outcome(o), weight(w) {
    std::replace(name.begin(), name.end(), '\t', ' ');
  }
  Observation(std::string nm, U o, float w = 1) : 
    name(nm), predictors(std::vector<T>()), outcome(o), weight(w) {
    std::replace(name.begin(), name.end(), '\t', ' ');
  }
  // ACCESSORS
  std::string get_name() const {return name;}
  std::vector<T>& get_predictors() {return predictors;}
  size_t size() const {return predictors.size();}
  const T& operator[](int i) const {return predictors[i];}
  U get_outcome() const {return outcome;}
  float get_weight() const {return weight;}
  
  // MUATATORS
  void set_weight(T t) {weight = t;}
  void set_outcome(U u) {outcome = u;}
  void add_feature(T t) {predictors.push_back(t);}

  static float accum_weight(float w, const Observation<T, U>& o) {
    return w + o.get_weight();
  }
  static float accum_outcome(float v, const Observation<T, U>& o) {
    return v + o.get_outcome();
  }
  static float accum_pos_weight(float w, const Observation<T, U>& o) {
    return w + ((o.get_outcome() > 0) ? o.get_weight() : 0);
  }
  
  // TODO: what is required here for in-class definition vs out-of-class
  // definition for friend functions w.r.t. template parameters?
  template <class TT, class UU> 
  friend std::ostream& operator<<(std::ostream &s, const Observation<TT, UU> &o) {
    s << o.name << "\t";
    copy(o.predictors.begin(), o.predictors.end(), 
	 std::ostream_iterator<float>(s, "\t"));
    return s << "\t" << o.outcome << "\t" << o.weight;
  }
  static void
  test_consistency(std::vector<Observation<T, U> >& obs) {
    for (typename std::vector<Observation<T, U> >::iterator i = obs.begin(); 
	 i != obs.end(); ++i)
      if (i->size() != obs.front().size()) 
	throw InconsistentObservations();
  }
  static void read_observations(std::string, std::vector<std::string>&, 
				std::vector<Observation<T, U> >&);

  static void shuffle(std::vector<Observation<T, U> >& obs) {
    srand(time(0) + getpid());
    random_shuffle(obs.begin(), obs.end());
  }
};

template <class T> T froma(std::string s) {
  std::istringstream ss(s);
  T t;
  ss >> t;
  return t;
}

template<class T, class U> void
Observation<T, U>::read_observations(std::string filename, std::vector<std::string> &features,
				     std::vector<Observation<T, U> >& obs) {
  std::ifstream in(filename.c_str());
  if (!in) { 
    std::string message = std::string("could not open observation file: ") + filename;
    throw ObservationFileException(message);
  }

  // allow comments at the beginning of the file
  size_t line_count = 0;
  char buffer[obs_line_size];
  while (in.peek() == '#') {
    in.getline(buffer, obs_line_size);
    ++line_count;
  }  

  // first non-comment line must be feature names
  features.clear();
  in.getline(buffer, obs_line_size);
  features = cread::split(buffer, " \t");
  
  size_t n_features = features.size();
  
  obs.clear();
  std::vector<std::string> v;
  while (!in.eof()) {
    in.getline(buffer, obs_line_size);
    ++line_count;
    std::vector<std::string> v = cread::split(buffer, "\t");
    if (v.size() < n_features + 2)
      throw ObservationFileException("observation file: bad observation line", 
				     line_count);
    U outcome = froma<U>(v[n_features + 1]);
    std::vector<T> pred;
    float weight = 1;
    if (v.size() == n_features + 3) {
      transform(v.begin() + 1, v.end() - 2,
		back_inserter(pred), std::ptr_fun(&froma<T>));
      weight = froma<float>(v.back());
    }
    else
      transform(v.begin() + 1, v.end() - 1,
		back_inserter(pred), std::ptr_fun(&froma<T>));
    obs.push_back(Observation<T, U>(v.front(), pred, outcome, weight));
    // TODO: figure out why I was using this peek() function:
    in.peek();
  }
}

/**
 * Function object used when sorting a sequence of observations based
 * on values of a specific predictor.
 */
template <class T, class U> 
class ObservationOrder : public std::binary_function<const Observation<T, U> &, 
						const Observation<T, U> &, bool> {
  size_t key;
public:
  explicit ObservationOrder(size_t k) : key(k) {}
  bool operator()(const Observation<T, U> &o1, const Observation<T, U> &o2) const {
    return o1[key] < o2[key];
  }
};

/**
 * Function object used to help identify a subset of predictors, from
 * a given sequence of predictors, with the value for a specific predictor
 * less than or equal to a specific value.
 */
template <class T, class U>
class ObservationLessEqual : public std::unary_function<const Observation<T, U>&, bool> {
  size_t key;
  T value;
public:
  explicit ObservationLessEqual(size_t k, T v) : key(k), value(v) {}
  bool operator()(const Observation<T, U>& o) const {
    return o[key] <= value;
  }
};

/**
 * Same as ObservationLessEqual, but provides the greater than or
 * equal to relation.
 */
template <class T, class U> 
class ObservationGreater : public std::unary_function<const Observation<T, U>&, bool> {
  size_t key;
  T value;
public:
  explicit ObservationGreater(size_t k, T v) : key(k), value(v) {}
  bool operator()(const Observation<T, U>& o) const {
    return o[key] > value;
  }
};

#endif
