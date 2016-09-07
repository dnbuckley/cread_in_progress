/*
 * Copyright (C) 2006 Cold Spring Harbor Laboratory
 * Authors: Andrew D. Smith, Dustin E. Schones, Pavel Sumazin and Michael Q. Zhang
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

#include "ModuleSite.hpp"

using std::string;
using std::ostringstream;
using std::vector;
using std::copy;
using std::sort;


ModuleSite::ModuleSite(const std::string& sn, int st, size_t l,
	     const std::vector<patternSite>& psv):
	     seq_name(sn), start(st), length(l), module_score(0.0){
  copy(psv.begin(), psv.end(), back_inserter(sites));
  sort(sites.begin(), sites.end());
}


ModuleSite::ModuleSite(const string& sn, int st, size_t l, float ms,
			const vector<patternSite>& psv):
			seq_name(sn), start(st), length(l), module_score(ms) {
  copy(psv.begin(), psv.end(), back_inserter(sites));
  sort(sites.begin(), sites.end());
}


ModuleSite::ModuleSite(const string &sn, const string &st, const string &l, 
		       const string &ms, const vector<string> &s,
		       const vector<string> &sts, const vector<string> &ls,
		       const vector<string> &o, const vector<string> &scs) :
		       seq_name(sn) {
  if (s.size() != sts.size() || s.size() != ls.size() || s.size() != o.size())
    throw InvalidModuleSiteException("Bad initialization: Sizes do not match");
  
  start = atoi(st.c_str());
  int temp_length = atoi(l.c_str());
  if (temp_length > 0)
    length = static_cast<size_t>(temp_length);
  else throw InvalidModuleSiteException("Bad initialization: negative length");
  module_score = atof(ms.c_str());

  for (size_t i = 0; i < sts.size(); ++i){
    if(o[i].size() != 1)
      throw InvalidModuleSiteException("Bad initialization: missing start-offset information");
    sites.push_back(patternSite(s[i], atoi(sts[i].c_str()), i,
			o[i].c_str()[0], atof(scs[i].c_str())));
  }
  sort(sites.begin(), sites.end());
}


ModuleSite::ModuleSite(const string &sn, const float &ms, const vector<string> &s,
		       const vector<int> &sts, const vector<float> &scs,
		       const vector<size_t> &ls, const vector<string> &o) :
  seq_name(sn), module_score(ms) {

  if (s.size() != sts.size() || s.size() != scs.size() || s.size() != ls.size() ||
      s.size() != o.size())
    throw InvalidModuleSiteException("Bad initialization: Sizes do not match");
  
  for (size_t i = 0; i < sts.size(); ++i) {
    if(o[i].size() != 1)
      throw InvalidModuleSiteException("Bad initialization: invalid orientation");
    sites.push_back(patternSite(s[i], sts[i], i, o[i].c_str()[0], scs[i]));
  }
  
  //identify ModuleSite extreme points
  int maxst = -std::numeric_limits<int>::max();
  start = std::numeric_limits<int>::max();
  for (size_t i = 0; i < s.size(); ++i) {
    maxst = std::max(maxst, static_cast<int>(sites[i].start + sites[i].site.size()));
    start = std::min(start, sites[i].start);
  }
  length = static_cast<size_t>(maxst - start);
  sort(sites.begin(), sites.end());
}


ModuleSite::ModuleSite(const string &sn, const vector<string> &s,
		       const vector<int> &sts, const vector<size_t> &ls,
		       const vector<string> &o, const vector<float> &scs) :
  seq_name(sn) {
  if (s.size() != sts.size() || s.size() != ls.size() || s.size() != o.size() ||
      s.size() != scs.size())
    throw InvalidModuleSiteException("Bad initialization: negative length");
  
  for (size_t i = 0; i < sts.size(); ++i) {
    if(o[i].size() != 1)
      throw InvalidModuleSiteException("Bad initialization: invalid orientation");
    sites.push_back(patternSite(s[i], sts[i], i, o[i].c_str()[0], scs[i]));
  }
  
  //identify ModuleSite extreme points
  int maxst = -std::numeric_limits<int>::max();
  start = std::numeric_limits<int>::max();
  for (size_t i = 0; i < s.size(); ++i) {
    maxst = std::max(maxst, static_cast<int>(sites[i].start + sites[i].site.size()));
    start = std::min(start, sites[i].start);
  }
  length = static_cast<size_t>(maxst - start);
  sort(sites.begin(), sites.end());
}


ModuleSite::ModuleSite(const string &sn, const int st, const size_t l, 
		       const float ms, const vector<string> &s, const vector<int> &sts,
		       const vector<size_t> &ls, const vector<string> &o,
		       const vector<float> &scs) :
  seq_name(sn), start(st), length(l), module_score(ms) {
  if (s.size() != sts.size() || s.size() != ls.size() || s.size() != o.size())
    throw InvalidModuleSiteException("Bad initialization: negative length");
  for (size_t i = 0; i < sts.size(); ++i) {
    if(o[i].size() != 1)
      throw InvalidModuleSiteException("Bad initialization: invalid orientation");
    sites.push_back(patternSite(s[i], sts[i], i, o[i].c_str()[0], scs[i]));
  }
  sort(sites.begin(), sites.end());
}


ModuleSite::ModuleSite(string s) {
  vector<string> parts(cread::split(s, ";", false));
  if (parts.size() < 5)
    throw InvalidModuleSiteException("Bad initialization: too few site components");
  seq_name = parts[1];
  start = atoi(parts[2].c_str());
  int temp_length = atoi(parts[3].c_str());
  if (temp_length > 0)
    length = static_cast<size_t>(temp_length);
  else throw InvalidModuleSiteException("Bad initialization: missing information");
  module_score = atof(parts[4].c_str());
  parts = cread::split(parts.front(), "+", true);

  for (size_t i = 0; i < parts.size(); ++i) {
    const vector<string> subparts(cread::split(
	      parts[i].substr(1, parts[i].length() - 2), ":", false));
    sites.push_back(patternSite(subparts[1], atoi(subparts.front().c_str()),
			   atoi(subparts[2].c_str()), subparts[3].c_str()[0],
					      atof(subparts.back().c_str())));
  }
}


string
ModuleSite::tostring() const {
  ostringstream os;
  os << "(" << sites[0].start << ":" << sites[0].site << ":" 
     << sites[0].pattern << ":" << sites[0].orientation << ":"
				      << sites[0].score << ")";
  for (size_t i = 1; i < sites.size(); ++i)
    os << "+(" << sites[i].start << ":" << sites[i].site << ":" 
       << sites[i].pattern << ":" << sites[i].orientation << ":"
					<< sites[i].score <<  ")";
  os << "; " << seq_name << "; " << start << "; " 
     << length << "; " << module_score <<";";
  return os.str();
}


vector<string> ModuleSite::get_sites() const {
  vector<string> s(sites.size());
  for (size_t i = 1; i < sites.size(); ++i)
    s[i] = sites[i].site;
  return s;
}


vector<size_t> ModuleSite::get_lengths() const {
  vector<size_t> s(sites.size());
  for (size_t i = 1; i < sites.size(); ++i)
    s[i] = sites[i].site.size();
  return s;
}


vector<int> ModuleSite::get_starts() const {
  vector<int> s(sites.size());
  for (size_t i = 1; i < sites.size(); ++i)
    s[i] = sites[i].start;
  return s;
}


vector<char> ModuleSite::get_orientations() const{
  vector<char> s(sites.size());
  for (size_t i = 1; i < sites.size(); ++i)
    s[i] = sites[i].orientation;
  return s;
}


vector<float> ModuleSite::get_scores() const {
  vector<float> s(sites.size());
  for (size_t i = 1; i < sites.size(); ++i)
    s[i] = sites[i].score;
  return s;
}
