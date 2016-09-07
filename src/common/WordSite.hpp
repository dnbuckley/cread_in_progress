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

#ifndef WORDSITE_HPP
#define WORDSITE_HPP

#include "cread.hpp"
#include "Pattern.hpp"

class WordSite {
public:
  WordSite(std::string s, std::string sn, int st, 
	   size_t l, std::string o) : 
    site(s), seq_name(sn), orientation(o), start(st), length(l) {}
  WordSite(std::string s, std::string sn, std::string st, 
	   std::string l, std::string o);
  WordSite(std::string);
  std::string tostring() const;
  friend std::ostream& operator<<(std::ostream &s, const WordSite &bs) {
    return s << bs.tostring();
  }
  
  std::string get_site() const {return site;}
  std::string get_seq_name() const {return seq_name;}
  std::string get_orientation() const {return orientation;}
  int get_start() const {return start;}
  size_t get_length() const {return length;}

private:
  // TODO: make the start and the length into numbers
  std::string site;
  std::string seq_name;
  //std::string gaps;
  std::string orientation;
  int start;
  size_t length;
};

/*!
  \class InvalidWordSiteException
  \brief Exception class for handling invalid word site exceptions
 */
class InvalidWordSiteException : public PatternSiteException {
public:
  //! Constructor that initializes the message
  InvalidWordSiteException(const std::string m = "") : PatternSiteException(m){}
};

#endif
