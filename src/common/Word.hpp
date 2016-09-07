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

#ifndef WORD_HPP
#define WORD_HPP
/*!
  \file Word.hpp
  \brief Contains the Word pattern and associated exception structures.
 */

#include "cread.hpp"
#include "WordSite.hpp"
#include "Pattern.hpp"
//#include "Alphabet.hpp"
#include "smithlab_utils.hpp"

/*!
  \class Word
  \brief The word pattern includes patterns consisting
	 of the 4 nucleotides and 'N'.
 */
class Word : public Pattern {
public:
  /*! \brief This constructor populates the pattern from a string vector
      and calls the base class Pattern to retrieve standard attributes.
   */
  Word(std::vector<std::string>&);
  //! This minimal constructor sets the word, accession and type
  Word(const std::string, const std::string);
  //! The copy constructor
  Word(const Word &w);
  //! The assignment operator
  Word& operator=(const Word&);

  // Accessors
  /*! \brief Returns a minimally initialized Word
	     with word set to the reverse complement of this word
  */
  Word revcomp_w() const;
  //! Returns the width of this word
  size_t get_width() const {return word.length();}
  //! Returns this word
  std::string get_word() const {return word;}
  //! Returns true if the string begins with the type Word label
  static bool is_type(std::string s) {
    return !s.compare(0, type_id_size, type_id);
  }

  //! Returns a vector of words retrieved from file_name
  static std::vector<Word> ReadWordVector(std::string file_name);

protected:
  //! The Word type id
  static const char *type_id;
  //! The length of the Word type label
  static const size_t type_id_size;
  //! The prefix of the Word type label
  static const char *word_start;
  //! Print the Word representation
  void format_representation(std::ostream& os) const;
  //! Print the sites
  void format_sites(std::ostream& os) const;
  //! The word pattern
  std::string word;
  //! The site vector
  std::vector<WordSite> sites;
};

/*!
  \exception WordFormatterException
  \brief Used for reporting exceptions while writing the pattern.
*/
class WordFormatterException : public PatternFormatterException {};
/*!
  \exception WordFormatterException
  \brief Used for reporting exceptions while reading the pattern.
*/
class WordFormatException : public PatternFormatException {};

#endif
