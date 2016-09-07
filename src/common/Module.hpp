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

#ifndef MODULE_HPP
#define MODULE_HPP
/*!
    \file Module.hpp
    \brief Contains the Module pattern and associated exception classes.
*/


/*!
  \file Module.hpp

  \brief This header file contains class definitions for the Module
  class, a subclass of Pattern.
*/

#include "cread.hpp"
#include "ModuleSite.hpp"
#include "Matrix.hpp"
#include "Pattern.hpp"

/*!
  \brief Represents cis-regulatory modules, which are sets of motifs
  whose sites may be under organizational constraints.

  \class Module
*/
class Module : public Pattern {
public:
  /// This constructor parses a vector of strings to populate the Module pattern
  Module(std::vector<std::string>&);
  /// The copy constructor
  Module(const Module&);
  /// Assignment operator
  Module& operator=(const Module&);
  /*!
      \brief This constructor creates a minimal motif pattern from a matrix;
	     the accession is set to (none).
  */
  Module(std::vector<Matrix>&, std::string);
  
  /// Returns the number of matrices
  size_t size() const {return matrices.size();}
  
  typedef Matrix& reference;
  typedef const Matrix& const_reference;
  typedef std::vector<Matrix>::iterator iterator;
  typedef std::vector<Matrix>::const_iterator const_iterator;

  /// Returns a const iterator pointing to the front of the matrix vector
  const_iterator begin() const {return matrices.begin();}
  /// Returns a const iterator pointing to the back of the matrix vector
  const_iterator end() const {return matrices.end();}
  
  /// Returns an iterator pointing to the front of the matrix vector
  iterator begin() {return matrices.begin();}
  /// Returns an iterator pointing to the back of the matrix vector
  iterator end() {return matrices.end();}
  
  /// Get the occurrences associated with this module
  std::vector<ModuleSite> get_sites() const {return sites;}
  
  /// Returns true if the type is a module type
  static bool is_type(std::string s) {
    return !s.compare(0, type_id_size, type_id);
  }
  
  /// The reference operator returns the nth matrix
  reference operator[](int n) {return matrices[n];}
  /// The constant reference operator returns the nth matrix
  const_reference operator[](int n) const {return matrices[n];}
  
  /// Returns the first matrix
  reference front() {return matrices.front();}
  /// Returns the last matrix
  reference back() {return matrices.back();} 
  
  /// Clear sites
  void clear_sites() {sites.clear();}
  /// Add a module sites
  void add_site(ModuleSite& bs) {sites.push_back(bs);}
  
  /// Read modules from a file
  static std::vector<Module> ReadModuleVector(std::string file_name);
  /*!
     \brief Returns a string consisting of a concatenation of l and index,
	    separated by '_'
     \param[in] l A character array.
     \param[in] index An unsigned integer.
     \return "%s_%s", l, index
  */
  static std::string add_index(const char *l, size_t index);
  
private:
  
  static const char *type_id;
  static const size_t type_id_size;
  
  // Allowing both 'PO' and 'P0' for now
  //static const char *matrix_binding_site_start;
  static const char *matrix_start;
  
  void format_representation(std::ostream& os) const;
  void format_sites(std::ostream& os) const;
  std::vector<Matrix> matrices;
  std::vector<ModuleSite> sites;
};

inline std::string
Module::add_index(const char *l, size_t index) {
  return std::string(l) + std::string("_") + cread::toa(index);
}

/*!
  \class ModuleFormatterException
  \brief Reports writing exceptions.
 */
class ModuleFormatterException : public PatternFormatterException {};
/*!
  \class ModuleFormatException
  \brief Reports reading exceptions.
 */
class ModuleFormatException : public PatternFormatException {};

#endif
