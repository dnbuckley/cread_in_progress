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

#include "Module.hpp"

using std::string;
using std::vector;
using std::pair;
using std::map;
using std::ostringstream;
using std::ostream_iterator;
using std::endl;
using std::ostream;
using std::copy;

const char *
Module::type_id = "Module";

const size_t
Module::type_id_size = 6;

const char *
Module::matrix_start = "P0";


Module::Module(const Module& m) : Pattern(m) {
  matrices = m.matrices;
  type = type_id;
  sites = m.sites;
}


Module&
Module::operator=(const Module& m) {
  if (this != &m) {
    Pattern::operator=(m);
    matrices = m.matrices;
  }
  return *this;
}

Module::Module(vector<Matrix>& m, string a) : matrices(m) {
  accession = a;
  type = type_id;
}

void
Module::format_representation(ostream& os) const {
  for (size_t i = 0; i < matrices.size(); ++i)
    os << matrices[i] << endl << PatternID::BLANK_PATTERN_LINE << endl;
}

void
Module::format_sites(ostream& os) const {
  string sep("\n");
  sep += PatternID::BINDING_SITE_START + string("  ");
  if (sites.size() > 0) {
    os << PatternID::BINDING_SITE_START << "  ";
    copy(sites.begin(), sites.end() - 1,
         ostream_iterator<ModuleSite>(os, sep.c_str()));
    os << sites.back() << endl
       << PatternID::BLANK_PATTERN_LINE << endl;
  }
}

Module::Module(vector<string>& lines) : Pattern(lines) {
  /* ADS: the read_matrix variable was not used. Probably this
     function should be reconsidered. */
  // bool read_matrix = false,
  bool read_sites = false;
  for (size_t i = 0; i < lines.size(); ++i) {
    string line(lines[i]);
    if (line_type(line, matrix_start)) {
      vector<string> matrix_lines;
      matrix_lines.push_back(line);
      for (i += 1; i < lines.size() && isdigit(lines[i][0]); ++i)
        matrix_lines.push_back(lines[i]);
      matrices.push_back(Matrix(matrix_lines));
      // read_matrix = true;
    }
    else if (!read_sites &&
             line_type(line, PatternID::BINDING_SITE_START)) {
      while (line_type(lines[i], PatternID::BINDING_SITE_START))
        sites.push_back(ModuleSite(remove_line_id(lines[i++])));
      read_sites = true;
    }
    // TODO: need to verify format better in here
  }
  if (matrices.empty())
    // TODO: reset state (deallocate matrices)
    throw ModuleFormatException();
  type = type_id;
}

vector<Module>
Module::ReadModuleVector(string file_name) {
  vector<vector<string> > module_lines;
  ReadPatternLines(file_name, module_lines);
  vector<Module> modules;
  for (size_t i = 0; i < module_lines.size(); ++i)
    modules.push_back(Module(module_lines[i]));
  return modules;
}
