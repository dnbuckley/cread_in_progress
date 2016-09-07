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

#ifndef BEDFILE_HPP
#define BEDFILE_HPP

/*!
  \file BEDFile.hpp

  \brief Declarations of functions related to reading and writing
  files in BED format, elements of which are represented as GenomicRegion
  or SimpleGenomicRegion objects.

  Information on this format can be found in the Help pages for the
  UCSC Genome Browser.
*/

#include "cread.hpp"

class GenomicRegion;
class SimpleGenomicRegion;
class UCSCGenomeBrowserHeader;
class UCSCGenomeBrowserTrack;

/***********************************************************************
 *
 *  Functions for reading bed files
 *
 */

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////                                                ////
/////                Multiple tracks                 ////
/////                                                ////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

/*!
  \brief read a BED file consisting of multiple tracks

  \param filename
  Name of the BED format file to read.
  \param the_tracks Reference to vector of genome browser tracks in
  which to pass back that information.
  \param the_regions Reference to vector of genomic regions in
  which to pass back that info.
  
*/
void
ReadBEDTracks(std::string filename,
	      std::vector<UCSCGenomeBrowserTrack> &the_tracks,
	      std::vector<std::vector<GenomicRegion> > &the_regions
	      );

/*!
  \brief read a BED file consisting of multiple tracks, and
  keep the file header.
  
  \param filename
  Name of the BED format file to read.
  \param the_tracks Reference to vector of genome browser tracks in
  which to pass back that information.
  \param the_regions Reference to vector of genomic regions in
  which to pass back that info.
  \param the_header Reference to a genome browser header in which
  that information will be passed back.
*/
void
ReadBEDTracks(std::string filename,
	      std::vector<UCSCGenomeBrowserTrack> &the_tracks,
	      std::vector<std::vector<GenomicRegion> > &the_regions,
	      UCSCGenomeBrowserHeader &the_header
	      );

/*!
  \brief read a BED file, consisting of multiple tracks, as simple
  genomic regions (chromosome, start, end).

  \param filename
  Name of the BED format file to read
  \param the_tracks reference to vector of genome browser tracks in which
  to pass back that information.
  \param the_regions Reference to vector of simple genomic regions in
  which to pass back that info.
*/
void
ReadBEDTracks(std::string filename,
	      std::vector<UCSCGenomeBrowserTrack> &the_tracks,
	      std::vector<std::vector<SimpleGenomicRegion> > &the_regions
	      );

/*!
  \brief read a BED file, consisting of multiple tracks, as
  simple genomic regions, and also keep the file header.

  \param filename
  Name of the BED format file to read
  \param the_tracks reference to vector of genome browser tracks in which
  to pass back that information.
  \param the_regions Reference to vector of simple genomic regions in
  which to pass back that info.
  \param the_header Reference to a genome browser header in which
  that information will be passed back.
  
*/
void
ReadBEDTracks(std::string filename,
	      std::vector<UCSCGenomeBrowserTrack> &the_tracks,
	      std::vector<std::vector<SimpleGenomicRegion> > &the_regions,
	      UCSCGenomeBrowserHeader &the_header
	      );

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//// SINGLE TRACK
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/*!
  \brief read a BED file assuming all lines are for a single track.

  \param filename
  Name of the BED format file to read
  \param the_track Reference to a genome browser track object in which
  to pass back the track information.
  \param the_regions Reference to vector of genomic regions in
  which to pass back that info.
*/
void
ReadBEDTrack(std::string filename,
	     UCSCGenomeBrowserTrack &the_track,
	     std::vector<GenomicRegion> &the_regions
	     );

/*!
  \brief read a BED file assuming all lines are for a single track,
  and keep the file header.

  \param filename
  Name of the BED format file to read
  \param the_track Reference to a genome browser track object in which
  to pass back the track information.
  \param the_regions Reference to vector of genomic regions in
  which to pass back that info.
  \param the_header Reference to a genome browser header in which
  that information will be passed back.
*/
void
ReadBEDTrack(std::string filename,
	     UCSCGenomeBrowserTrack &the_track,
	     std::vector<GenomicRegion> &the_regions,
	     UCSCGenomeBrowserHeader &header
	     );

/*!
  \brief read a BED file, assuming all lines are for a single track,
  as simple genomic regions.

  \param filename
  Name of the BED format file to read
  \param the_track Reference to a genome browser track object in which
  to pass back the track information.
  \param the_regions Reference to vector of simple genomic regions in
  which to pass back that info.
*/
void
ReadBEDTrack(std::string filename,
	     UCSCGenomeBrowserTrack &the_track,
	     std::vector<SimpleGenomicRegion> &the_regions
	     );

/*!
  \brief read a BED file, assuming all lines are for a single track,
  as simple genomic regions, and keep the file header.

  \param filename
  Name of the BED format file to read
  \param the_track Reference to a genome browser track object in which
  to pass back the track information.
  \param the_regions Reference to vector of simple genomic regions in
  which to pass back that info.
  \param the_header Reference to a genome browser header in which
  that information will be passed back.
*/
void
ReadBEDTrack(std::string filename,
	     UCSCGenomeBrowserTrack &the_track,
	     std::vector<SimpleGenomicRegion> &the_regions,
	     UCSCGenomeBrowserHeader &header
	     );

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//// JUST THE REGIONS
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/*!
  \brief read a BED file, but do not read track information.
  
  \param filename
  Name of the BED format file to read
  \param the_regions Reference to vector of genomic regions in
  which to pass back that info.
  
*/
void
ReadBEDFile(std::string filename,
	    std::vector<GenomicRegion> &the_regions
	    );

/*!
  \brief read a BED file as simple genomic regions, but do not keep
  the track information.

  \param filename
  Name of the BED format file to read
  \param the_regions Reference to vector of simple genomic regions in
  which to pass back that info.
  
*/
void
ReadBEDFile(std::string filename,
	    std::vector<SimpleGenomicRegion> &the_regions
	    );

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
	       std::vector<std::vector<GenomicRegion> > &the_regions
	       );

void
WriteBEDTracks(std::string filename,
	       std::vector<UCSCGenomeBrowserTrack> &the_tracks,
	       std::vector<std::vector<GenomicRegion> > &the_regions,
	       UCSCGenomeBrowserHeader &the_header
	       );

void
WriteBEDTracks(std::string filename,
	       std::vector<UCSCGenomeBrowserTrack> &the_tracks,
	       std::vector<std::vector<SimpleGenomicRegion> > &the_regions
	       );

void
WriteBEDTracks(std::string filename,
	       std::vector<UCSCGenomeBrowserTrack> &the_tracks,
	       std::vector<std::vector<SimpleGenomicRegion> > &the_regions,
	       UCSCGenomeBrowserHeader &the_header
	       );

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//// SINGLE TRACK
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

void
WriteBEDTrack(std::string filename,
	      UCSCGenomeBrowserTrack &the_track,
	      std::vector<GenomicRegion> &the_regions
	      );

void
WriteBEDTrack(std::string filename,
	      UCSCGenomeBrowserTrack &the_track,
	      std::vector<GenomicRegion> &the_regions,
	      UCSCGenomeBrowserHeader &header
	      );

void
WriteBEDTrack(std::string filename,
	      UCSCGenomeBrowserTrack &the_track,
	      std::vector<SimpleGenomicRegion> &the_regions
	      );

void
WriteBEDTrack(std::string filename,
	      UCSCGenomeBrowserTrack &the_track,
	      std::vector<SimpleGenomicRegion> &the_regions,
	      UCSCGenomeBrowserHeader &header
	      );

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//// JUST THE REGIONS
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

void
WriteBEDFile(std::string filename,
	     std::vector<GenomicRegion> &the_regions
	     );

void
WriteBEDFile(std::string filename,
	     std::vector<SimpleGenomicRegion> &the_regions
	     );

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////// UCSCGenomeBrowserHeader class ///////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

/**
 * Representation of UCSC Genome Browser custom annotation track file
 * headers. These headers must appear at the beginning of the file
 * (before any "tracks") and each line must start with the word
 * "browser". The format of these browser lines is as follows:
 *
 * browser attribute_name attribute_value(s)
 *
 * The valid attribute names and values are (from UCSC Genome Browser Help):
 *
 * position <position> - Determines the part of the genome that the
 * Genome Browser will initially open to, in chromosome:start-end
 * format.
 *
 * pix <width> - Sets the Genome Browser window to the specified width
 * in pixels.
 * 
 * hide all - Hides all annotation tracks except for custom ones.
 * hide <track_name(s)> - Hides the listed tracks. Multiple track
 * names should be space-separated.
 * 
 * dense all - Displays all tracks in dense mode.
 * dense <track_name(s)> - Displays the specified tracks in dense
 * mode. Symbolic names must be used. Multiple track names should be
 * space-separated.
 *  
 * pack all - Displays all tracks in pack mode.
 * pack <track_name(s)> - Displays the specified tracks in pack
 * mode. Symbolic names must be used. Multiple track names should be
 * space-separated.
 *
 * squish all - Displays all tracks in squish mode. NOTE: If the
 * browser display includes a large number of tracks or a large
 * position range, this option may overload your browser's resources
 * and cause an error or timeout.
 * squish <track_name(s)> - Displays the specified tracks in squish
 * mode. Symbolic names must be used. Multiple track names should be
 * space-separated.
 *
 * full all - Displays all tracks in full mode. See NOTE for squish
 * all.
 *
 * full <track_name(s)> - Displays the specified tracks in full
 * mode. Symbolic names must be used. Multiple track names should be
 * space-separated.
 */
class UCSCGenomeBrowserHeader {
private:
  static const char* browser_line_id; 
  static const size_t browser_line_id_len;
  
  std::vector<std::string> header_lines;
  size_t pix;
public:
  /**
   * Construct an empty header to be filled later
   */
  UCSCGenomeBrowserHeader() {}
  /**
   \brief Construct a header with lines that are valid browser lines.
   
   \param lines 
   a vector of valid browser lines
   */
  UCSCGenomeBrowserHeader(const std::vector<std::string>& lines);
  
  /**
     \brief Get a string representation (for printing) of this header.
     
     \return
     the string representation
   */
  std::string tostring() const;
  
  // mutators
  
  /**
     \brief Set the initial position when the annotation file is loaded.
     
     \param 
     region a GenomicRegion whose values specify the initial
     position of the browser window
  */
  void set_position(const GenomicRegion& region);

  /**
     \brief Set the initial position when the annotation file is loaded.
     
     \param 
     region a SimpleGenomicRegion whose values specify the
     initial position of the browser window
  */
  void set_position(const SimpleGenomicRegion& region);
};

std::ostream& 
operator<<(std::ostream& s, 
	   const UCSCGenomeBrowserHeader& the_header);

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////////////////// UCSCGenomeBrowserTrack class ///////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
 	
/**
 * Class representing the "track" lines that give information about
 * custom tracks in the UCSC Genome Browser. These immediately precede
 * the lines in a custom annotation track file that specify the actual
 * annotations. All must be printed on a line, and that line begins
 * with the word "track". Each attribute-value pair must be separated
 * by whitespace, and if whitespace is within an attribute value, that
 * value must be quoted. The attributes and their values must appear as:
 *
 * attribute_name=<attribute_value>
 *
 * The valid names and the values they can take are as follows (from
 * the UCSC Genome Browser Help):
 * 
 * name=<track_label> - up to 15 characters, and must be
 * enclosed in quotes if contains spaces. Default is "User Track".
 * 
 * description=<center_label> - Center label of track in GB window. Up
 * to 60 characters, quotes required if contains spaces. Default="User
 * Supplied Track".
 * 
 * visibility=<display_mode> - display mode of track. Values: 0 - hide,
 * 1 - dense, 2 - full, 3 - pack, 4 - squish.  The numerical values or
 * words can be used. Default="1".
 * 
 * color=<RRR,GGG,BBB> - Defines the main color for the annotation
 * track. The track color consists of three comma-separated RGB values
 * from 0-255. The default value is 0,0,0 (black).
 * 
 * itemRgb=On - If present and set to "On", Browser will use RGB value
 * shown in the itemRgb field in each data line of the associated BED
 * track to determine the display color of the data on that
 * line. (works for 9 or less fields)
 * 
 * useScore=<use_score> - If present and set to 1, score field in
 * each line will determine the level of shading. Color attribute
 * is 100,50,0 (shades of brown) or 0,60,120 (shades of blue).
 * Default setting for useScore is "0".
 * 
 * group=<group> - Defines the annotation track group in which the
 * custom track will display in the Genome Browser window. By default,
 * group is set to "user", which causes custom tracks to display at the
 * top of the window.
 * 
 * priority=<priority> - When the group attribute is set, defines the
 * display position of the track relative to other tracks within the
 * same group in the Genome Browser window. If group is not set, the
 * priority attribute defines the track's order relative to other
 * custom tracks displayed in the default group, "user".
 * 
 * offset=<offset> - Defines a number to be added to all coordinates in
 * the annotation track. The default is "0".
 * 
 * url=<external_url> - Defines a URL for an external link associated
 * with this track. used in details page for track. 
 * '$$' in will be substituted with item name. no default
 * 
 * htmlUrl=<external_url> - Defines a URL for an HTML description page
 * to be displayed with this track. no default
 */
class UCSCGenomeBrowserTrack {
private:
  static const size_t n_valid_attributes = 11;
  static const char *valid_attributes[];
  //   = {
  //     "description",
  //     "visibility",
  //     "color",
  //     "itemRgb",
  //     "useScore",
  //     "group",
  //     "priority",
  //     "offset",
  //     "url",
  //     "htmlUrl"
  //   };
  
  // instance variables
  std::map<std::string, std::string> attributes;
  std::string name;
  
  static bool is_valid_attribute_label(std::string attr);
  
public:
  UCSCGenomeBrowserTrack() {};
  UCSCGenomeBrowserTrack(const std::string& track_line);
  
  // mutators
  template<class LABEL_VALUE_TYPE> void set_attribute(std::string name, 
						      LABEL_VALUE_TYPE val);
  void set_name(std::string n) {name = n;}
  
  // accessors
  std::string get_name() const {return name;}
  std::string get_attribute(std::string name) const;
  std::string tostring() const;
};


class BEDFileException : public CREADException {
public:
  BEDFileException(std::string s = std::string()) throw() : CREADException(s) {}
};


template<class LABEL_VALUE_TYPE>
void
UCSCGenomeBrowserTrack::set_attribute(std::string label, LABEL_VALUE_TYPE val) {
  if (is_valid_attribute_label(label))
    attributes[label] = cread::toa(val);
  else throw BEDFileException("invalid track attribute label: " + label);
}

std::ostream& operator<<(std::ostream& s, 
			 const UCSCGenomeBrowserTrack& the_track);

#endif
