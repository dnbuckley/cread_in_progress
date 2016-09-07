/*
 * Copyright (C) 2006 Cold Spring Harbor Laboratory
 * Authors: Andrew D. Smith, Julie M Granka and Michael Q. Zhang
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

#include "BEDFile.hpp"
#include "GenomicRegion.hpp"

#include <cassert>

using std::string;
using std::vector;
using std::ostringstream;
using std::ostream_iterator;
using std::endl;
using std::map;
using std::numeric_limits;
using std::cerr;
using std::ofstream;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////////////// UCSCGenomeBrowserHeader functions ////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

const char *UCSCGenomeBrowserHeader::browser_line_id = "browser";
const size_t UCSCGenomeBrowserHeader::browser_line_id_len = 7;

UCSCGenomeBrowserHeader::UCSCGenomeBrowserHeader(const vector<string>& lines) {
  for (size_t i = 0; i < lines.size(); ++i) {
    vector<string> parts;
    cread::split_whitespace(lines[i], parts);
    if (parts.front() != browser_line_id)
      throw BEDFileException("bad browser line: " + lines[i]);
    header_lines.push_back(lines[i]);
  }
}


string 
UCSCGenomeBrowserHeader::tostring() const {
  ostringstream s;
  if (!header_lines.empty()) {
    copy(header_lines.begin(), header_lines.end() - 1,
	 std::ostream_iterator<string>(s, "\n"));
    s << header_lines.back();
  }
  return s.str();
}


std::ostream& operator<<(std::ostream& s, 
			 const UCSCGenomeBrowserHeader& the_header) {
  return s << the_header.tostring();
}


void
UCSCGenomeBrowserHeader::set_position(const GenomicRegion& region) {
  static const char *position = "position";
  ostringstream ss;
  ss << position << "\t" 
     << region.get_chrom() << "\t" 
     << region.get_start() << "\t" << region.get_end();
  header_lines.push_back(ss.str());
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////// UCSCGenomeBrowserTrack functions ////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

const char *
UCSCGenomeBrowserTrack::valid_attributes[] = {
  "description",
  "visibility",
  "color",
  "itemRgb",
  "useScore",
  "group",
  "priority",
  "offset",
  "url",
  "htmlUrl"
};


string
UCSCGenomeBrowserTrack::get_attribute(string label) const {
  if (attributes.find(label) != attributes.end()) 
    return attributes.find(label)->second;
  else return "";
}


bool
UCSCGenomeBrowserTrack::is_valid_attribute_label(string attr) {
  for (size_t i = 0; i < n_valid_attributes; ++i)
    if (attr == valid_attributes[i]) return true;
  return false;
};


UCSCGenomeBrowserTrack::UCSCGenomeBrowserTrack(const string& track_line) {
  vector<string> parts = cread::split_whitespace_quoted(track_line);
  
  assert(!parts.empty() && parts[0] == "track" && 
	 "\"track\" not at start of track line");
  
  for (size_t i = 1; i < parts.size(); ++i) {
    const size_t label_end = parts[i].find_first_of("=");
    const string label(parts[i].substr(0, label_end));
    if (label != "name" && !is_valid_attribute_label(label))
      throw BEDFileException("invalid track attribute label: " + label);
    size_t valstart = label_end + 1;
    size_t valend = parts[i].length();
    if (parts[i][valstart] == '\"') {
      valstart += 1;
      valend -= 1;
    }
    const string value(parts[i].substr(valstart, valend - valstart));
    
    if (label == "name")
      name = value;
    else attributes[label] = value;
  }
}


string 
UCSCGenomeBrowserTrack::tostring() const {
  ostringstream s;
  if (!name.empty() || !attributes.empty())
    s << "track";
  if (!name.empty())
    s << " name=" << name;
  typedef map<string, string>::const_iterator attr_itr;
  for (attr_itr i = attributes.begin(); i != attributes.end(); ++i) {
    s << " " << i->first << "=";
    const bool print_quotes = (i->second.find(' ') != string::npos ||
			       i->second.find('\t') != string::npos);
    if (print_quotes) s << "\"";
    s << i->second;
    if (print_quotes) s << "\"";
  }
  return s.str();
}


std::ostream& 
operator<<(std::ostream& s, 
	   const UCSCGenomeBrowserTrack& the_track) {
  return s << the_track.tostring();
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////////////////////// Non-member functions /////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


// utilities to check line types
static bool 
is_header_line(const string& line) {
  return line.substr(0, 7) == string("browser");
}


static bool 
is_track_line(const string& line) {
  return line.substr(0, 5) == string("track");
}


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////                                                ////
/////                MULTIPLE TRACKS                 ////
/////                                                ////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


void
ReadBEDTracks(string filename,
	      vector<UCSCGenomeBrowserTrack> &the_tracks,
	      vector<vector<GenomicRegion> > &the_regions,
	      UCSCGenomeBrowserHeader &the_header) {
  static const size_t buffer_size = 1000; // Magic
  
  // open and check the file
  std::ifstream in(filename.c_str());
  if (!in) 
    throw BEDFileException("cannot open input file " + filename);

  // separate according to tracks and strip off any header lines
  vector<string> header_lines;
  vector<string> track_lines;
  
  bool no_tracks = false;
  
  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw BEDFileException("Line too long in file: " + filename);
    const string line(buffer);
    if (is_header_line(line))
      header_lines.push_back(line);
    else if (is_track_line(line)) {
      if (no_tracks)
	throw BEDFileException("track line missing for first track in " + 
			       filename);
      track_lines.push_back(line);
      the_regions.push_back(vector<GenomicRegion>());
    }
    else {
      if (track_lines.empty() && !no_tracks) {
	no_tracks = true;
	the_regions.push_back(vector<GenomicRegion>());
      }
      the_regions.back().push_back(GenomicRegion(line));
    }
    in.peek();
  }
  in.close();
  
  for (size_t i = 0; i < track_lines.size(); ++i) {
    the_tracks.push_back(UCSCGenomeBrowserTrack(track_lines[i]));
  }
  
  assert(((the_regions.size() == 1 && the_tracks.empty()) ||
	  the_regions.size() == the_tracks.size()) &&
	 "the_regions.size() != the_tracks.size()");
  
  the_header = UCSCGenomeBrowserHeader(header_lines);
}


void
ReadBEDTracks(string filename,
	      vector<UCSCGenomeBrowserTrack> &the_tracks,
	      vector<vector<GenomicRegion> > &the_regions) {
  UCSCGenomeBrowserHeader the_header;
  ReadBEDTracks(filename, the_tracks, the_regions, the_header);
}


void
ReadBEDTracks(string filename,
	      vector<UCSCGenomeBrowserTrack> &the_tracks,
	      vector<vector<SimpleGenomicRegion> > &the_regions,
	      UCSCGenomeBrowserHeader &the_header) {
  static const size_t buffer_size = 1000; // Magic
  
  // open and check the file
  std::ifstream in(filename.c_str());
  if (!in) 
    throw BEDFileException("cannot open input file " + filename);

  // separate according to tracks and strip off any header lines
  vector<string> header_lines;
  vector<string> track_lines;

  bool no_tracks = false;

  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw BEDFileException("Line too long in file: " + filename);
    const string line(buffer);
    if (is_header_line(line))
      header_lines.push_back(line);
    else if (is_track_line(line)) {
      if (no_tracks)
	throw BEDFileException("track line missing for first track in " + 
			       filename);
      track_lines.push_back(line);
      the_regions.push_back(vector<SimpleGenomicRegion>());
    }
    else {
      if (the_regions.empty()) {
	no_tracks = true;
	the_regions.push_back(vector<SimpleGenomicRegion>());
      }
      the_regions.back().push_back(SimpleGenomicRegion(line));
    }
    in.peek();
  }
  in.close();
  
  for (size_t i = 0; i < track_lines.size(); ++i) {
    the_tracks.push_back(UCSCGenomeBrowserTrack(track_lines[i]));
  }
  
  assert((the_regions.size() == 1 && the_tracks.empty()) ||
	 the_regions.size() == the_tracks.size() &&
	 "the_regions.size() != the_tracks.size()");
  
  the_header = UCSCGenomeBrowserHeader(header_lines);
}


void
ReadBEDTracks(string filename,
	      vector<UCSCGenomeBrowserTrack> &the_tracks,
	      vector<vector<SimpleGenomicRegion> > &the_regions) {
  UCSCGenomeBrowserHeader the_header;
  ReadBEDTracks(filename, the_tracks, the_regions, the_header);
}


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//// SINGLE TRACK
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////


void
ReadBEDTrack(string filename,
	     UCSCGenomeBrowserTrack &the_track,
	     vector<GenomicRegion> &the_regions,
	     UCSCGenomeBrowserHeader &the_header) {
  
  vector<UCSCGenomeBrowserTrack> the_tracks;
  vector<vector<GenomicRegion> > the_region_sets;
  ReadBEDTracks(filename, the_tracks, the_region_sets, the_header);
  
  if (!the_tracks.empty())
    the_track = the_tracks.front();
  
  vector<vector<GenomicRegion> >::const_iterator i;
  for (i = the_region_sets.begin(); i != the_region_sets.end(); ++i)
    the_regions.insert(the_regions.end(), i->begin(), i->end());
}


void
ReadBEDTrack(string filename,
	     UCSCGenomeBrowserTrack &the_tracks,
	     vector<GenomicRegion> &the_regions) {
  UCSCGenomeBrowserHeader the_header;
  ReadBEDTrack(filename, the_tracks, the_regions, the_header);
}


void
ReadBEDTrack(string filename,
	     UCSCGenomeBrowserTrack &the_track,
	     vector<SimpleGenomicRegion> &the_regions,
	     UCSCGenomeBrowserHeader &the_header) {
  vector<UCSCGenomeBrowserTrack> the_tracks;
  vector<vector<SimpleGenomicRegion> > the_region_sets;
  ReadBEDTracks(filename, the_tracks, the_region_sets, the_header);
  
  if (!the_tracks.empty())
    the_track = the_tracks.front();
  
  vector<vector<SimpleGenomicRegion> >::const_iterator i;
  for (i = the_region_sets.begin(); i != the_region_sets.end(); ++i)
    the_regions.insert(the_regions.end(), i->begin(), i->end());
}


void
ReadBEDTrack(string filename,
	     UCSCGenomeBrowserTrack &the_track,
	     vector<SimpleGenomicRegion> &the_regions) {
  UCSCGenomeBrowserHeader the_header;
  ReadBEDTrack(filename, the_track, the_regions, the_header);
}


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//// JUST THE REGIONS
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////


void
ReadBEDFile(string filename, vector<GenomicRegion> &the_regions) {
  static const size_t buffer_size = 1000; // Magic
  
  // open and check the file
  std::ifstream in(filename.c_str());
  if (!in) 
    throw BEDFileException("cannot open input file " + filename);
  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw BEDFileException("Line too long in file: " + filename);
    const string line(buffer);
    if (!is_header_line(line) && !is_track_line(line)) {
      the_regions.push_back(GenomicRegion(line));
    }
    in.peek();
  }
  in.close();
}


void
ReadBEDFile(string filename,
	    vector<SimpleGenomicRegion> &the_regions) {
  static const size_t buffer_size = 1000; // Magic
  
  // open and check the file
  std::ifstream in(filename.c_str());
  if (!in) 
    throw BEDFileException("cannot open input file " + filename);
  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw BEDFileException("Line too long in file: " + filename);
    const string line(buffer);
    if (!is_header_line(line) && !is_track_line(line)) {
      the_regions.push_back(SimpleGenomicRegion(line));
    }
    in.peek();
  }
  in.close();
}

/***********************************************************************
 *
 *  Functions for writing bed files
 *
 */

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////                                                ////
/////                Multiple tracks                 ////
/////                                                ////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

void
WriteBEDTracks(std::string filename,
	       std::vector<UCSCGenomeBrowserTrack> &the_tracks,
	       std::vector<std::vector<GenomicRegion> > &the_regions,
	       UCSCGenomeBrowserHeader &the_header) {
  assert(the_tracks.size() == the_regions.size() &&
	 "the_tracks.size() is not equal to the_regions.size()");
  ofstream out(filename.c_str());
  if (!out)
    throw BEDFileException("could not open file for writing: " + filename);
  out << the_header << std::endl;
  for (size_t i = 0; i < the_tracks.size(); ++i)
    out << the_tracks[i] << std::endl
	<< the_regions[i] << std::endl;
  out.close();
}

void
WriteBEDTracks(std::string filename,
	       std::vector<UCSCGenomeBrowserTrack> &the_tracks,
	       std::vector<std::vector<GenomicRegion> > &the_regions) {
  assert(the_tracks.size() == the_regions.size() &&
	 "the_tracks.size() is not equal to the_regions.size()");
  ofstream out(filename.c_str());
  if (!out)
    throw BEDFileException("could not open file for writing: " + filename);
  for (size_t i = 0; i < the_tracks.size(); ++i)
    out << the_tracks[i] << std::endl
	<< the_regions[i] << std::endl;
  out.close();
}

void
WriteBEDTracks(std::string filename,
	       std::vector<UCSCGenomeBrowserTrack> &the_tracks,
	       std::vector<std::vector<SimpleGenomicRegion> > &the_regions) {
  assert(the_tracks.size() == the_regions.size() &&
	 "the_tracks.size() is not equal to the_regions.size()");
  ofstream out(filename.c_str());
  if (!out)
    throw BEDFileException("could not open file for writing: " + filename);
  for (size_t i = 0; i < the_tracks.size(); ++i)
    out << the_tracks[i] << std::endl
	<< the_regions[i] << std::endl;
  out.close();
}


void
WriteBEDTracks(std::string filename,
	       std::vector<UCSCGenomeBrowserTrack> &the_tracks,
	       std::vector<std::vector<SimpleGenomicRegion> > &the_regions,
	       UCSCGenomeBrowserHeader &the_header) {
  assert(the_tracks.size() == the_regions.size() &&
	 "the_tracks.size() is not equal to the_regions.size()");
  ofstream out(filename.c_str());
  if (!out)
    throw BEDFileException("could not open file for writing: " + filename);
  out << the_header << std::endl;
  for (size_t i = 0; i < the_tracks.size(); ++i)
    out << the_tracks[i] << std::endl
	<< the_regions[i] << std::endl;
  out.close();
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//// SINGLE TRACK
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

void
WriteBEDTrack(std::string filename,
	      UCSCGenomeBrowserTrack &the_track,
	      std::vector<GenomicRegion> &the_regions) {
  ofstream out(filename.c_str());
  if (!out)
    throw BEDFileException("could not open file for writing: " + filename);
  out << the_track << std::endl
      << the_regions << std::endl;
  out.close();
}

void
WriteBEDTrack(std::string filename,
	      UCSCGenomeBrowserTrack &the_track,
	      std::vector<GenomicRegion> &the_regions,
	      UCSCGenomeBrowserHeader &the_header) {
  ofstream out(filename.c_str());
  if (!out)
    throw BEDFileException("could not open file for writing: " + filename);
  out << the_header << std::endl
      << the_track << std::endl
      << the_regions << std::endl;
  out.close();
}


void
WriteBEDTrack(std::string filename,
	      UCSCGenomeBrowserTrack &the_track,
	      std::vector<SimpleGenomicRegion> &the_regions) {
  ofstream out(filename.c_str());
  if (!out)
    throw BEDFileException("could not open file for writing: " + filename);
  out << the_track << std::endl
      << the_regions << std::endl;
  out.close();
}


void
WriteBEDTrack(std::string filename,
	      UCSCGenomeBrowserTrack &the_track,
	      std::vector<SimpleGenomicRegion> &the_regions,
	      UCSCGenomeBrowserHeader &the_header) {
  ofstream out(filename.c_str());
  if (!out)
    throw BEDFileException("could not open file for writing: " + filename);
  out << the_header << std::endl
      << the_track << std::endl
      << the_regions << std::endl;
  out.close();
}


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//// JUST THE REGIONS
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

void
WriteBEDFile(std::string filename,
	     std::vector<GenomicRegion> &the_regions) {
  ofstream out(filename.c_str());
  if (!out)
    throw BEDFileException("could not open file for writing: " + filename);
  out << the_regions << std::endl;
  out.close();
}


void
WriteBEDFile(std::string filename,
	     std::vector<SimpleGenomicRegion> &the_regions) {
  ofstream out(filename.c_str());
  if (!out)
    throw BEDFileException("could not open file for writing: " + filename);
  out << the_regions << std::endl;
  out.close();
}

