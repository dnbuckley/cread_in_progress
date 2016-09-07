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

#ifndef FEATURE_HPP
#define FEATURE_HPP

#include "Pattern.hpp"
#include "SuffixTree.hpp"
#include "cread.hpp"

class feat_val {
public:
  virtual std::string srep() const = 0;
  virtual ~feat_val() {}
};

class Feature {
public:
  virtual feat_val* evaluate(const SuffixTree&) const = 0;
  virtual feat_val* evaluate(const std::string&) const = 0;
  virtual std::string name() const = 0;
  virtual ~Feature() {}
};

struct feat_val_double : public feat_val {
  double val;
  feat_val_double(double v = 0.0) : val(v) {}
  std::string srep() const {
    std::ostringstream ss;
    ss << val;
    return ss.str();
  }
};

struct feat_val_float : public feat_val {
  float val;
  feat_val_float(float v = 0.0) : val(v) {}
  std::string srep() const {
    std::ostringstream ss;
    ss << val;
    return ss.str();
  }
};

struct feat_val_bool : public feat_val {
  bool val;
  feat_val_bool(bool v = false) : val(v) {}
  std::string srep() const {
    std::ostringstream ss;
    ss << val;
    return ss.str();
  }
};

struct feat_val_int : public feat_val {
  int val;
  feat_val_int(int v = 0) : val(v) {}
  std::string srep() const {
    std::ostringstream ss;
    ss << val;
    return ss.str();
  }
};

#endif
