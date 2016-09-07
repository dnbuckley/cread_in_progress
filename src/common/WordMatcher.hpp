/*
 * Copyright (C) 2006 Cold Spring Harbor Laboratory
 * Authors: Andrew D. Smith and Michael Q. Zhang
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

#ifndef WORDMATCHER_H
#define WORDMATCHER_H

#include "cread.hpp"
#include "Word.hpp"

void
get_matches(const Word& word, const std::vector<std::string>& seqs,
	    const size_t min_matching_positions,
	    std::vector<std::string>& matches,
	    size_t max_matches = std::numeric_limits<size_t>::max());

size_t
count_sequences_with_match(const Word& word, 
			   const std::vector<std::string>& seqs,
			   const size_t min_matching_positions);

#endif
