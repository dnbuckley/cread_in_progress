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

#include "Word.hpp"

using std::string;
using std::vector;
using std::ostringstream;
using std::ostream;
using std::endl;

const char *
Word::type_id = "Word";

const size_t
Word::type_id_size = 4;

const char *
Word::word_start = "WO";

Word::Word(const Word& w) : Pattern(w), word(w.word) {}

Word& 
Word::operator=(const Word& w) {
  if (this != &w) {
    Pattern::operator=(w);
    word = w.word;
  }
  return *this;
}

Word::Word(const string w, const string a) : word(w) {
  accession = a;
  type = type_id;
}

Word::Word(vector<string>& lines) : Pattern(lines) {
  bool read_word = false;
  for (size_t i = 0; i < lines.size(); ++i) {
    string line = lines[i];
    if (Pattern::line_type(line, word_start)) {
      word = Pattern::remove_line_id(line);
      read_word = true;
    }
    // TODO: need to verify format better in here
  }
  type = type_id;
  if (!read_word)
    throw WordFormatException();
}


void
Word::format_representation(ostream& os) const {
  os << word_start << "  " << word << endl
     << PatternID::BLANK_PATTERN_LINE << endl;
}

void
Word::format_sites(ostream& os) const {
  string sep("\n");
  sep += PatternID::BINDING_SITE_START + string("  ");
  if (sites.size() > 0) {
    os << PatternID::BINDING_SITE_START << "  ";
    std::copy(sites.begin(), sites.end() - 1,
	      std::ostream_iterator<WordSite>(os, sep.c_str()));
    os << sites.back() << endl;
    os << PatternID::BLANK_PATTERN_LINE << endl;
  }
}

vector<Word> 
Word::ReadWordVector(string file_name) {
  vector<vector<string> > word_lines;
  ReadPatternLines(file_name, word_lines);
  vector<Word> words;
  for (size_t i = 0; i < word_lines.size(); ++i)
    words.push_back(Word(word_lines[i]));
  return words;
}

Word
Word::revcomp_w() const {
  return Word(revcomp(this->word), accession + string("_rc"));
}
