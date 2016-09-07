/*
 * Copyright (C) 2007 Cold Spring Harbor Laboratory
 * Authors: Andrew D. Smith
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

#include "GenomeUtils.hpp"

using std::string;
using std::vector;

static const char *digits = "0987654321";
static const char *whitespace = " \t";


// define what to do if parts are missing: if no ':' or '-', is whole
// thing chrom?

void
parse_region_name(string region_name, 
		  string& chrom, size_t &start, size_t &end) {
  
  const size_t colon_offset = region_name.find(":");
  
  // get the chromosome
  size_t chr_offset = region_name.find_last_of(whitespace, colon_offset);
  if (chr_offset == string::npos)
    chr_offset = 0;
  else
    chr_offset += 1;
  chrom = region_name.substr(chr_offset, colon_offset - chr_offset);
  
  // get the start
  const size_t start_end = region_name.find("-", colon_offset + 1);
  const string start_string = region_name.substr(colon_offset + 1,
						 start_end - colon_offset + 1);
  start = static_cast<size_t>(atoi(start_string.c_str()));
  
  // get the end
  const size_t end_end =
    region_name.find_first_not_of(digits, start_end + 1);
  const string end_string = region_name.substr(start_end + 1,
					       end_end - start_end - 1);
  end = static_cast<size_t>(atoi(end_string.c_str()));
}
