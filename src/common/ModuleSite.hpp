/*!
 * Copyright (C) 2006 Cold Spring Harbor Laboratory
 * \author Andrew D. Smith, Pavel Sumazin and Michael Q. Zhang
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

#ifndef MODULESITE_HPP
#define MODULESITE_HPP

#include "cread.hpp"
#include "Pattern.hpp"

/*!
 * \class ModuleSite
 * \brief Provides data structures and operators to read, write and store
 *	  module sites (cis-regulatory modules).
 */
class ModuleSite {
public:
  /*!
   * \struct patternSite
   * \brief Provides storage for sites associated with individual
   *	    pattern components of the module.
   *
   *  \param site
   *	The sequence of the pattern site
   *  \param start
   *	The start offset relative to the left most position in the sequence
   *  \param pattern
   *	The index of the pattern corresponding to the site
   *  \param orientation
   *	The orientation is 'p' or 'n'
   *  \param score
   *	The score of a site is optional with default value 0
   */
  struct patternSite{
    //! The site sequence
    std::string site;
    //! The start offset
    int start;
    //! The index of the associated pattern in the module
    size_t pattern;
    //! The orientation ('n' or 'p')
    char orientation;
    //! The score
    float score;
    //! A constructor
    patternSite(const std::string site, int start, size_t pattern,
			      char orientation, float score=0):
		site(site), start(start), pattern(pattern),
		orientation(orientation), score(score) {}
    //! A copy constructor
    patternSite(const patternSite& ps):
		site(ps.site), start(ps.start), pattern(ps.pattern),
		orientation(ps.orientation), score(ps.score) {}
    //! The less operator
    bool operator<(const patternSite& ps) const { return start < ps.start; }
    //! The equality operator
    bool operator==(const patternSite& ps)
       { return (ps.site == site && ps.start == start && ps.pattern == pattern
				 && ps.orientation == orientation) ? true :
								     false; }
    //! The inequality operator uses the equality operator
    bool operator!=(const patternSite& ps) { return (*this == ps) ? false :
								    true; }
  };

  //! A constructor that expects full initialization
  ModuleSite(const std::string& sn,
	     int st,
	     size_t l,
	     float ms,
	     const std::vector<patternSite>& psv);
  //! A constructor with score initialized to 0
  ModuleSite(const std::string& sn,
	     int st,
	     size_t l,
	     const std::vector<patternSite>& psv);
  /*! \brief A constructor that takes vectors of strings, positions, start sites
   * and lengths.
   *
   * This constructor will be deprecated once dependencies are resolved.
   * \warning To be deprecated
   */
  ModuleSite(const std::string& sn,
	     int st,
	     size_t l,
	     float ms,
	     const std::vector<std::string>& s,
	     const std::vector<int>& sts,
	     const std::vector<size_t>& ls,
	     const std::vector<std::string>& o,
	     const std::vector<float>& scs);
  /*! \brief A constructor that takes vectors of strings, positions, start sites
   * and lengths, and computes length and and start sites from the vectors.
   * The module score is not set.
   *
   * This constructor will be deprecated once dependencies are resolved.
   * \warning To be deprecated
   */
  ModuleSite(const std::string &sn,                 // sequence id
	     const std::vector<std::string>& s,    // sites
	     const std::vector<int>& sts,           // starts
	     const std::vector<size_t>& ls,         // lengths
	     const std::vector<std::string>& o,     // orientations
	     const std::vector<float>& scs);        // scores
  /*! \brief A constructor that takes vectors of strings, positions, start sites
   * and lengths, and computes length and and start sites from the vectors.
   *
   * This constructor will be deprecated once dependencies are resolved.
   * \warning To be deprecated
   */
  ModuleSite(const std::string &sn,  // seq id  
	     const float &ms,        // module score
	     const std::vector<std::string> &s, // sites
	     const std::vector<int> &sts,       // starts
	     const std::vector<float> &scs,     // scores
	     const std::vector<size_t> &ls,     // lengths
	     const std::vector<std::string> &o); // orientations
  /*! \brief A constructor that takes vectors of strings, positions, start sites
   * and lengths.  All information is taken as strings.
   *
   * This constructor will be deprecated once dependencies are resolved.
   * \warning To be deprecated
   */
  ModuleSite(const std::string &sn,
	     const std::string &st,
	     const std::string &l,
	     const std::string &ms,                // module score
	     const std::vector<std::string>& s,
	     const std::vector<std::string>& sts,
	     const std::vector<std::string>& ls,
	     const std::vector<std::string>& o,
	     const std::vector<std::string>& scs);
  /*! \brief A constructor that parses a string to initialize the module site.
   *
   * The input string must have at least 5 components, giving
   * \li List of pattern sites enclosed in parenthesis and separated
   *	 by '+'.  Each has the form (offset:sequence:patternIndex:Orientation:score)
   * \li Sequence name
   * \li Module-site start offset
   * \li Module-site length
   * \li Module-site score
   *
   * \param [in] s input string
   *
   */
  ModuleSite(std::string s);

  //! Represent module data in string form
  std::string tostring() const;
  //! Prints modules
  friend std::ostream& operator<<(std::ostream& s, const ModuleSite &ms) {
    return s << ms.tostring();
  }
  
  //! Returns the number of sites
  size_t size() const {return sites.size();}

  //! Returns the sequence name
  std::string get_seq_name() const {return seq_name;}
  //! Returns the module start site
  int get_start() const {return start;}
  /*! \brief
   * Returns the module-site length which is computed as the distance
   * from the beginning of the left most site to the end of the right
   * most site in the module.
   */
  size_t get_length() const {return length;}
  //! Returns the module score
  float get_module_score() const {return module_score;}
  
  //! Sets the module score
  void set_module_score(float s){module_score = s;};


  // individual sites
  //! Returns a vector of the pattern-site sequences
  std::vector<std::string> get_sites() const;
  //! Returns a vector of the pattern-site lengths
  std::vector<size_t> get_lengths() const;
  //! Returns a vector of the pattern-site start offsets
  std::vector<int> get_starts() const;
  //! Returns a vector of the pattern-site orientations
  std::vector<char> get_orientations() const;
  //! Returns a vector of the pattern-site scores
  std::vector<float> get_scores() const;
  
  // parts of the site for a particular module
  //! Returns the sequence of the pattern-site indexed by index
  std::string get_site(size_t index) const {return sites[index].site;}
  //! Returns the start offset of the pattern-site indexed by index
  size_t get_start(size_t index) const {return sites[index].start;}
  //! Returns the size of the pattern-site indexed by index
  size_t get_length(size_t index) const {return sites[index].site.size();}
  //! Returns the orientation of the pattern-site indexed by index
  char get_orientation(size_t index) const { return sites[index].orientation; }
  //! Returns the score of the pattern-site indexed by index
  float get_score(size_t index) const {return sites[index].score;}
  
  //! Sets the score of the pattern-site indexed by index
  void set_score(size_t index, float s) {sites[index].score = s;}

private:
  //! Sequence name
  std::string seq_name;
  //! Start offset for the module site
  int start;
  //! Total length of the module site from left most site to end of right most
  size_t length;
  //! Module score
  float module_score;
  //! Pattern sites
  std::vector<patternSite> sites;
};

/*!
   \class InvalidModuleSiteException
   \brief Exception class for handling invalid module site exceptions
*/
class InvalidModuleSiteException : public PatternSiteException {
  public:
  //! Constructor that initializes the message
  InvalidModuleSiteException(const std::string m) : PatternSiteException(m){}
};


#endif
