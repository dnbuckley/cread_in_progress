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

#include "WordSite.hpp"

using std::string;
using std::ostringstream;

WordSite::WordSite(std::string s, std::string sn, std::string st, 
	   std::string l, std::string o) :
  site(s), seq_name(sn), orientation(o) {
  start = atoi(st.c_str());
  int temp_length = atoi(l.c_str());
  if (temp_length > 0)
    length = static_cast<size_t>(temp_length);
  else throw InvalidWordSiteException();
}

WordSite::WordSite(string s) {
  size_t seq_name_offset = s.find_first_of(";");
  site = s.substr(0, seq_name_offset);
  seq_name_offset = s.find_first_not_of(" ", seq_name_offset + 1);
  
  size_t start_offset = s.find_first_of(";", seq_name_offset);
  seq_name = s.substr(seq_name_offset, start_offset - seq_name_offset);
  start_offset = s.find_first_not_of(" ", start_offset + 1);
  
  size_t length_offset = s.find_first_of(";", start_offset);
  start = atoi(s.substr(start_offset, length_offset - start_offset).c_str());
  length_offset = s.find_first_not_of(" ", length_offset + 1);
  
  //size_t gaps_offset = s.find_first_of(";", length_offset + 1);
  size_t orientation_offset = s.find_first_of(";", length_offset + 1);
  int temp_length = atoi(s.substr(length_offset, 
				  orientation_offset - length_offset).c_str());
  //   int temp_length = atoi(s.substr(length_offset, 
  // 				  gaps_offset - length_offset).c_str());
  //gaps_offset = s.find_first_not_of(" ", gaps_offset + 1);
  
  //size_t orientation_offset = s.find_first_of(";", gaps_offset);
  //gaps = s.substr(gaps_offset, orientation_offset - gaps_offset);
  orientation_offset = s.find_first_not_of(" ", orientation_offset + 1);
  
  size_t end_offset = s.find_first_of(";.", orientation_offset);
  orientation = s.substr(orientation_offset, end_offset - orientation_offset);
  
  if (temp_length >= 0)
    length = static_cast<size_t>(temp_length);
  else throw InvalidWordSiteException();
}

string
WordSite::tostring() const {
  ostringstream os;
  os << site << "; " << seq_name << "; " << start << "; " 
     << length << "; " << orientation;
  return os.str();
}
