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

#include "GenomeAlignment.hpp"
//#include "Alphabet.hpp"
#include "BEDFile.hpp"
#include "smithlab_utils.hpp"

#include <cassert>

using std::string;
using std::vector;
using std::ostringstream;
using std::ostream;
using std::numeric_limits;
using std::endl;


using std::endl;
using std::cerr;

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////        AlignmentBlock class         ///////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

string
AlignmentBlock::tostring() const {
  assert(src.size() == aln.size() && 
	 "src.size() != aln.size() in AlignmentBlock::tostring()");
  
  vector<size_t> info_lengths;
  size_t max_length = 0;
  for (vector<GenomicRegion>::const_iterator i = src.begin();
       i != src.end(); ++i) {
    ostringstream s;
    s << i->get_name() << " " << i->get_chrom() << " " 
      << i->get_start() << " " << i->get_end() << " " 
      << i->get_strand();
    info_lengths.push_back(s.str().length());
    max_length = std::max(info_lengths.back(), max_length);
  }
  
  ostringstream s;
  vector<GenomicRegion>::const_iterator i = src.begin();
  vector<string>::const_iterator j = aln.begin();
  vector<size_t>::const_iterator k = info_lengths.begin();
  while (i != src.end() - 1) {
    s << i->get_name() << " " << i->get_chrom() << " " 
      << i->get_start() << " " << i->get_end() << " " 
      << i->get_strand() << string(max_length - *k + 1, ' ') << *j << endl;
    ++i; ++j; ++k;
  }
  s << i->get_name() << " " << i->get_chrom() << " " 
    << i->get_start() << " " << i->get_end() << " " 
    << i->get_strand() << string(max_length - *k + 1, ' ') << *j;
  return s.str();
}


string
AlignmentBlock::maf_block_string() const {
  assert(src.size() == aln.size() && 
	 "src.size() != aln.size() in AlignmentBlock::tostring()");

  vector<size_t> info_lengths;
  size_t max_length = 0;
  for (vector<GenomicRegion>::const_iterator i = src.begin();
       i != src.end(); ++i) {
    ostringstream s;
    s << "s " 
      << i->get_name() << '.' << i->get_chrom() << " " 
      << i->get_start() << " " << i->get_end() - i->get_start() << " " 
      << i->get_strand() << " " << i->get_score();
    info_lengths.push_back(s.str().length());
    max_length = std::max(info_lengths.back(), max_length);
  }

  ostringstream s;
  vector<GenomicRegion>::const_iterator i = src.begin();
  vector<string>::const_iterator j = aln.begin();
  vector<size_t>::const_iterator k = info_lengths.begin();
  s << 'a' << endl;
  while (i != src.end()) {
    s << "s " 
      << i->get_name() << '.' << i->get_chrom() << " " 
      << i->get_start() << " " << i->get_end() - i->get_start() << " " 
      << i->get_strand() << " " << i->get_score()
      << string(max_length - *k + 1, ' ')
      << *j << endl;
    ++i; ++j; ++k;
  }
  return s.str();
}


ostream& 
operator<<(ostream& s, const AlignmentBlock& the_block) {
  return s << the_block.tostring();
}


bool
AlignmentBlock::contains(const GenomicRegion &query) const {
  return src.front().contains(query);
}


size_t
AlignmentBlock::count_add_gaps(const string& sequence, size_t offset) {
  size_t count = 0;
  const char *base = sequence.c_str();
  const char *lookup = base;
  while (count < offset)
    count += (*(lookup++) != gap_symbol);
  while (*lookup == gap_symbol) lookup++;
  return lookup - base;
}


size_t
AlignmentBlock::count_sub_gaps(const string& sequence, size_t offset) {
  size_t count = 0;
  const char *lookup = sequence.c_str();
  const char *limit = lookup + offset;
  while (lookup < limit)
    count += (*(lookup++) != gap_symbol);
  while (*lookup == gap_symbol) lookup++;
  return count;
}


void
AlignmentBlock::get_slice(const GenomicRegion &query,
			  vector<string> &alnseqs) const {
  assert(contains(query) && 
	 "Attempting to get slice that does not exist:");
  
  const size_t start =
    count_add_gaps(aln.front(), query.get_start() - src.front().get_start());
  const size_t width =
    count_add_gaps(aln.front(), query.get_end() - src.front().get_start()) - start;
  
  for (vector<string>::const_iterator i = aln.begin(); i != aln.end(); ++i)
    alnseqs.push_back(i->substr(start, width));
}


void
AlignmentBlock::get_slice(const GenomicRegion &query,
 			  vector<GenomicRegion> &species_info,
			  vector<string> &alnseqs) const {
  assert(contains(query) && "Attempting to get slice that does not exist");
  
  const size_t ref_start = 
    count_add_gaps(aln.front(), query.get_start() - src.front().get_start());
  
  const size_t ref_width = 
    count_add_gaps(aln.front(),
		   query.get_end() - src.front().get_start()) - ref_start;
  
  for (size_t i = 0; i < aln.size(); ++i) {
    alnseqs.push_back(aln[i].substr(ref_start, ref_width));
    species_info.push_back(src[i]);
    species_info.back().set_start(src[i].get_start() +
				  count_sub_gaps(aln[i], ref_start));
    species_info.back().set_end(src[i].get_start() +
				count_sub_gaps(aln[i], ref_start + ref_width));
  }
}

string
AlignmentBlock::get_sequence(const GenomicRegion &query) const {
  assert(contains(query) &&
	 "Attempting to get segment that does not exist");
  
  const size_t start =
    count_add_gaps(aln.front(), query.get_start() - src.front().get_start());
  
  const size_t width =
    count_add_gaps(aln.front(), 
		   query.get_end() - src.front().get_start()) - start;
  
  return aln.front().substr(start, width);
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////        GenomeAlignment class         //////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


GenomeAlignment::GenomeAlignment(const vector<vector<GenomicRegion> > &regions,
				 const vector<vector<string> > &seqs) {
  for (size_t i = 0; i < regions.size(); ++i)
    blocks.push_back(AlignmentBlock(regions[i], seqs[i]));

  transform(blocks.begin(), blocks.end(), back_inserter(block_index),
	    std::mem_fun_ref(&AlignmentBlock::get_region));
}

inline static void
parse_block_line(const char *s,
		 const size_t len,
		 string &spec, 
		 string &chrom, 
		 size_t &start, 
		 size_t &end,
		 char &strand, 
		 string &seq) {
  size_t i = 0;
  
  // the species
  while (isspace(s[i]) && i < len) ++i;
  size_t j = i;
  while (!isspace(s[i]) && i < len) ++i;
  if (j < i) spec = string(s + j, i - j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // the chromosome
  while (isspace(s[i]) && i < len) ++i;
  j = i;
  while (!isspace(s[i]) && i < len) ++i;
  if (j < i) chrom = string(s + j, i - j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // start of the region (a positive integer)
  while (isspace(s[i]) && i < len) ++i;
  j = i;
  while (!isspace(s[i]) && i < len) ++i;
  if (j < i) start = atoi(s + j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // end of the region (a positive integer)
  while (isspace(s[i]) && i < len) ++i;
  j = i;
  while (!isspace(s[i]) && i < len) ++i;
  if (j < i) end = atoi(s + j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // the strand, as single character ('+' or '-')
  while (isspace(s[i]) && i < len) ++i;
  j = i;
  while (!isspace(s[i]) && i < len) ++i;
  if (j < i) strand = s[j];
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // the actual alignment sequence
  while (isspace(s[i]) && i < len) ++i;
  j = i;
  while (!isspace(s[i]) && i < len) ++i;
  if (j < i)
    seq = string(s + j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
}

inline static bool
is_space_or_tab(const char &c) {
  return c == ' ' || c == '\t';
}

inline static void
parse_maf_block_line(const char *s,
		     const size_t len,
		     string &spec, 
		     string &chrom, 
		     size_t &start, 
		     size_t &end,
		     char &strand, 
		     string &seq) {
  size_t i = 2;
  
  // s name.chrom start length strand junk sequence
  
  // the species
  while (s[i] != '.' && i < len) ++i;
  if (s[i] != '.')
    throw GenomeAlignmentException("invalid line:\n" + string(s));
  spec = string(s + 2, i - 2);//s + i - i - j);
  ++i;
  // else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // the chromosome
  while (is_space_or_tab(s[i]) && i < len) ++i;
  size_t j = i;
  while (!is_space_or_tab(s[i]) && i < len) ++i;
  if (j < i) chrom = string(s + j, i - j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // start of the region (a positive integer)
  while (is_space_or_tab(s[i]) && i < len) ++i;
  j = i;
  while (!is_space_or_tab(s[i]) && i < len) ++i;
  if (j < i) start = atoi(s + j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // end of the region (a positive integer)
  while (is_space_or_tab(s[i]) && i < len) ++i;
  j = i;
  while (!is_space_or_tab(s[i]) && i < len) ++i;
  if (j < i) 
    end = start + atoi(s + j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // the strand, as single character ('+' or '-')
  while (is_space_or_tab(s[i]) && i < len) ++i;
  j = i;
  while (!is_space_or_tab(s[i]) && i < len) ++i;
  if (j < i) strand = s[j];
  else throw GenomeAlignmentException("invalid line:\n" + string(s));

  // the actual alignment sequence
  while (is_space_or_tab(s[i]) && i < len) ++i;
  j = i;
  while (!is_space_or_tab(s[i]) && i < len) ++i;
  if (j < i)
    string jnk = string(s + j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
  
  // the actual alignment sequence
  while (is_space_or_tab(s[i]) && i < len) ++i;
  j = i;
  while (!is_space_or_tab(s[i]) && i < len) ++i;
  if (j < i)
    seq = string(s + j);
  else throw GenomeAlignmentException("invalid line:\n" + string(s));
}


// GenomeAlignment::GenomeAlignment(const string filename) {
//   std::ifstream in(filename.c_str());
//   if (!in) 
//     throw GenomeAlignmentException("cannot open alignment file " + filename);
  
//   vector<GenomicRegion> regions;
//   vector<string> seqs;
//   while (!in.eof()) {
    
//     char buffer[input_buffer_size];
//     in.getline(buffer, input_buffer_size);
    
//     if (!(*buffer)) {
      
//       if (!block_index.empty() && regions.front() < block_index.back())
// 	throw GenomeAlignmentException("out of order:\n" +
// 				       regions.front().tostring() + "\n" +
// 				       block_index.back().tostring());
//       blocks.push_back(AlignmentBlock(regions, seqs));
//       block_index.push_back(regions.front());
      
//       regions.clear();
//       seqs.clear();
//     }
//     else {
//       string spec, chrom, seq;
//       size_t start = 0, end = 0;
//       char strand = '+'; 
//       assert(static_cast<size_t>(in.gcount()) == strlen(buffer) + 1);
//       parse_block_line(buffer, in.gcount() - 1,
// 		       spec, chrom, start, end, strand, seq);
//       regions.push_back(GenomicRegion(chrom, start, end, spec, 0, strand));
//       seqs.push_back(seq);
//     }
//   }
// }

GenomeAlignment::GenomeAlignment(const string filename) {
  // open the input file
  std::ifstream in(filename.c_str());
  if (!in) 
    throw GenomeAlignmentException("cannot open input file " + filename);
  
  vector<GenomicRegion> regions;
  vector<string> seqs;

  while (!in.eof()) {

    char buffer[input_buffer_size];
    in.getline(buffer, input_buffer_size);
    //     if (in.gcount() == input_buffer_size - 1)
    //       throw GenomeAlignmentException("Line in " + filename + 
    // 				     "\nexceeds max length: " +
    // 				     cread::toa(input_buffer_size - 1));
    
    const size_t buffer_size = in.gcount() - 1;
    // assert(strlen(buffer) == buffer_size);
    // correct for dos carriage returns before newlines
    if (buffer[buffer_size - 1] == '\r')
      buffer[buffer_size - 1] = '\0';
    // get all the block lines and split according to blocks
    
    if (buffer[0] == 'a') {
      if (!block_index.empty() && regions.front() < block_index.back())
	throw GenomeAlignmentException("out of order:\n" +
				       regions.front().tostring() + "\n" +
				       block_index.back().tostring());
      if (!regions.empty()) {
	blocks.push_back(AlignmentBlock(regions, seqs));
	block_index.push_back(regions.front());

	regions.clear();
	seqs.clear();
      }
    }
    else if (buffer[0] == 's') {
      string spec, chrom, seq;
      size_t start = 0, end = 0;
      char strand = '+'; 
      // assert(static_cast<size_t>(in.gcount()) == strlen(buffer) + 1);
      parse_maf_block_line(buffer, in.gcount() - 1,
			   spec, chrom, start, end, strand, seq);
      regions.push_back(GenomicRegion(chrom, start, end, spec, 0, strand));
      seqs.push_back(seq);
    }
    in.peek();
  }
  in.close();
  
  if (!regions.empty()) {
    blocks.push_back(AlignmentBlock(regions, seqs));
    block_index.push_back(regions.front());
  }
}


string 
GenomeAlignment::tostring() const {
  ostringstream ss;
  if (!blocks.empty()) {
    vector<AlignmentBlock>::const_iterator i;
    for (i = blocks.begin(); i != blocks.end() - 1; ++i)
      ss << *i << endl << endl;
    ss << blocks.back();
  }
  return ss.str();
}


size_t
GenomeAlignment::containing_block(const GenomicRegion &query) const {
  vector<GenomicRegion>::const_iterator i = find_closest(block_index, query);
  if (i != block_index.end() && i->contains(query))
    return (i - block_index.begin());
  else return not_found;
}


bool
GenomeAlignment::contains(const GenomicRegion &query) const {
  size_t curr_index = containing_block(GenomicRegion(query.get_chrom(),
						     query.get_start(),
						     (query.get_start() + 1)));
  const size_t end_index = containing_block(GenomicRegion(query.get_chrom(),
							  (query.get_end() - 1),
							  query.get_end()));

  if (curr_index == not_found || end_index == not_found)
    return false;
  
  while (curr_index != end_index && 
	 (block_index[curr_index].get_end() == 
	  block_index[curr_index + 1].get_start()))
    ++curr_index;
  return (curr_index == end_index);
}


void
GenomeAlignment::get_slice(const GenomicRegion& query,
			   vector<string> &segments) const {
  if (!contains(query)) return;
  
  size_t curr_index = containing_block(GenomicRegion(query.get_chrom(),
						     query.get_start(),
						     (query.get_start() + 1)));
  const size_t end_index = containing_block(GenomicRegion(query.get_chrom(),
							  (query.get_end() - 1),
							  query.get_end()));
  typedef std::map<string, size_t> OrderMap;
  
  OrderMap spec_order;
  size_t ref_len = 0;
  for (size_t i = curr_index; i <= end_index; ++i) {
    
    vector<string> chunks;
    vector<GenomicRegion> spec_info;
    blocks[i].get_slice(intersection(block_index[i], query), spec_info, chunks);
    
    for (size_t j = 0; j < chunks.size(); ++j) {
      const string spec_name(spec_info[j].get_name());
      OrderMap::iterator k = spec_order.find(spec_name);
      if (k == spec_order.end()) {
	spec_order[spec_name] = segments.size();
	segments.push_back(string(ref_len, gap_symbol));
	segments.back() += chunks[j];
      }
      else {
	const size_t spec_id = k->second;
	segments[spec_id].append(ref_len - segments[spec_id].length(), gap_symbol);
	segments[spec_id].append(chunks[j]);
      }
    }
    ref_len = segments.front().length();
  }
  for (OrderMap::iterator i = spec_order.begin(); i != spec_order.end(); ++i)
    segments[i->second].append(ref_len - segments[i->second].length(), gap_symbol);
}


void
GenomeAlignment::get_slice(const GenomicRegion& query,
			   vector<GenomicRegion> &regions,
			   vector<string> &segments) const {
  
  if (!contains(query)) return;

  size_t curr_index = containing_block(GenomicRegion(query.get_chrom(),
						     query.get_start(),
						     (query.get_start() + 1)));
  const size_t end_index = containing_block(GenomicRegion(query.get_chrom(),
							  (query.get_end() - 1),
							  query.get_end()));
  typedef std::map<string, size_t> OrderMap;

  OrderMap spec_order;
  size_t ref_len = 0;
  for (size_t i = curr_index; i <= end_index; ++i) {
    
    vector<string> chunks;
    vector<GenomicRegion> spec_info;
    blocks[i].get_slice(intersection(block_index[i], query), spec_info, chunks);
    
    for (size_t j = 0; j < chunks.size(); ++j) {
      const string spec_name(spec_info[j].get_name());
      OrderMap::iterator k = spec_order.find(spec_name);
      if (k == spec_order.end()) {
	spec_order[spec_name] = segments.size();
	segments.push_back(string(ref_len, gap_symbol));
	segments.back() += chunks[j];
	regions.push_back(spec_info[j]);
      }
      else {
	const size_t spec_id = k->second;
	segments[spec_id].append(ref_len - segments[spec_id].length(), gap_symbol);
	segments[spec_id].append(chunks[j]);
	if (regions[spec_id].get_chrom() == spec_info[j].get_chrom()) {
	  regions[spec_id].set_start(std::min(regions[spec_id].get_start(),
					      spec_info[j].get_start()));
	  regions[spec_id].set_end(std::max(regions[spec_id].get_end(),
					    spec_info[j].get_end()));
	}
      }
    }
    ref_len = segments.front().length();
  }
  for (OrderMap::iterator i = spec_order.begin(); i != spec_order.end(); ++i)
    segments[i->second].append(ref_len - segments[i->second].length(), gap_symbol);
}

string 
GenomeAlignment::get_sequence(const GenomicRegion& query) const {
  if (!contains(query)) return string();
  size_t curr_index = containing_block(GenomicRegion(query.get_chrom(),
						     query.get_start(),
						     (query.get_start() + 1)));
  const size_t end_index = containing_block(GenomicRegion(query.get_chrom(),
							  (query.get_end() - 1),
							  query.get_end()));
  string r;
  for (size_t i = curr_index; i <= end_index; ++i)
    r += blocks[i].get_sequence(intersection(block_index[i], query));
  return r;
}

ostream& 
operator<<(ostream& s, const GenomeAlignment& aln) {
  return s << aln.tostring();
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////////////////////  GenomeAlignmentOnDisk  /////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void
GenomeAlignmentOnDisk::ReadIndexFile(const string index_filename,
				     vector<size_t> &offsets) {
  static const size_t buffer_size = 1000; // Magic
  
  // open and check the file
  std::ifstream in(index_filename.c_str());
  if (!in)
    throw GenomeAlignmentOnDiskException("cannot open input file " + 
					 index_filename);
  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw GenomeAlignmentOnDiskException("Line too long in file: " + 
					   index_filename);
    const string line(buffer);
    
    vector<string> parts;
    cread::split_whitespace(line, parts);
    if (parts.size() < 5)
      throw GenomeAlignmentOnDiskException("ERROR: bad index file line: " + line);
    
    offsets.push_back(static_cast<size_t>(atol(parts[4].c_str())));
    in.peek();
  }
  in.close();
}


GenomeAlignmentOnDisk::GenomeAlignmentOnDisk(const string index_filename, 
					     const string aln_filename) :
  current_block_index(std::numeric_limits<size_t>::max()),
  alnfile(aln_filename.c_str(), std::ios_base::binary),
  alnfile_name(aln_filename) {
  ReadBEDFile(index_filename, block_index);
  ReadIndexFile(index_filename, offsets);
}
  
std::string
GenomeAlignmentOnDisk::tostring() const {
  ostringstream s;
  for (size_t i = 0; i < block_index.size(); ++i) {
    s << block_index[i].get_chrom() << "\t"
      << block_index[i].get_start() << "\t"
      << block_index[i].get_end() << "\t"
      << block_index[i].get_name() << "\t"
      << block_index[i].get_score() << "\t"
      << block_index[i].get_strand() << endl;
  }
  return s.str();
}


size_t
GenomeAlignmentOnDisk::containing_block(const GenomicRegion &query) const {
  vector<GenomicRegion>::const_iterator i = find_closest(block_index, query);
  if (i != block_index.end() && i->contains(query))
    return (i - block_index.begin());
  else return not_found;
}


bool
GenomeAlignmentOnDisk::contains(const GenomicRegion &query) const {
  size_t curr_index = containing_block(GenomicRegion(query.get_chrom(),
						     query.get_start(),
						     (query.get_start() + 1)));
  const size_t end_index = containing_block(GenomicRegion(query.get_chrom(),
							  (query.get_end() - 1),
							  query.get_end()));
  
  if (curr_index == not_found || end_index == not_found)
    return false;
  
  while (curr_index != end_index && 
	 (block_index[curr_index].get_end() == 
	  block_index[curr_index + 1].get_start()))
    ++curr_index;
  return (curr_index == end_index);
}

AlignmentBlock
GenomeAlignmentOnDisk::get_block(const size_t offset_in_file) const {
  alnfile.seekg(offset_in_file);
  
  vector<GenomicRegion> regions;
  vector<string> seqs;
  
  while (alnfile.peek() != '\n' && !alnfile.eof()) {
    char buffer[input_buffer_size];
    alnfile.getline(buffer, input_buffer_size);
    if (buffer[0] == 's') {
      if (alnfile.gcount() == input_buffer_size - 1)
	throw GenomeAlignmentOnDiskException("Line too long at offset: " +
					     cread::toa(offset_in_file) +
					     " in file: " +
					     cread::toa(alnfile_name));
      
      string spec, chrom, seq;
      size_t start = 0, end = 0;
      char strand = '+'; 
      // assert(static_cast<size_t>(alnfile.gcount()) == strlen(buffer) + 1);
      parse_maf_block_line(buffer, alnfile.gcount() - 1,
			   spec, chrom, start, end, strand, seq);
      regions.push_back(GenomicRegion(chrom, start, end, spec, 0, strand));
      seqs.push_back(seq);
    }
  }
  return AlignmentBlock(regions, seqs);
}

void
GenomeAlignmentOnDisk::get_slice(const GenomicRegion& query,
				 vector<string> &segments) const {
  if (!contains(query)) return;
  
  size_t curr_index = containing_block(GenomicRegion(query.get_chrom(),
						     query.get_start(),
						     (query.get_start() + 1)));
  const size_t end_index = containing_block(GenomicRegion(query.get_chrom(),
							  (query.get_end() - 1),
							  query.get_end()));
  typedef std::map<string, size_t> OrderMap;
  
  OrderMap spec_order;
  size_t ref_len = 0;
  for (size_t i = curr_index; i <= end_index; ++i) {
    
    vector<string> chunks;
    vector<GenomicRegion> spec_info;
    if (i != current_block_index) {
      const size_t offset_in_file = offsets[i];
      current_block = get_block(offset_in_file);
      current_block_index = i;
    }
    current_block.get_slice(intersection(block_index[i], query), spec_info, chunks);

    for (size_t j = 0; j < chunks.size(); ++j) {
      const string spec_name(spec_info[j].get_name());
      OrderMap::iterator k = spec_order.find(spec_name);
      if (k == spec_order.end()) {
	spec_order[spec_name] = segments.size();
	segments.push_back(string(ref_len, gap_symbol));
	segments.back() += chunks[j];
      }
      else {
	const size_t spec_id = k->second;
	segments[spec_id].append(ref_len - segments[spec_id].length(), gap_symbol);
	segments[spec_id].append(chunks[j]);
      }
    }
    ref_len = segments.front().length();
  }
  for (OrderMap::iterator i = spec_order.begin(); i != spec_order.end(); ++i)
    segments[i->second].append(ref_len - segments[i->second].length(), gap_symbol);
}


void
GenomeAlignmentOnDisk::get_slice(const GenomicRegion& query,
				 std::vector<GenomicRegion> &regions,
				 std::vector<std::string> &segments) const {
  
  if (!contains(query)) return;
  size_t curr_index = containing_block(GenomicRegion(query.get_chrom(),
						     query.get_start(),
						     (query.get_start() + 1)));
  const size_t end_index = containing_block(GenomicRegion(query.get_chrom(),
							  (query.get_end() - 1),
							  query.get_end()));
  typedef std::map<string, size_t> OrderMap;
  
  OrderMap spec_order;
  size_t ref_len = 0;
  for (size_t i = curr_index; i <= end_index; ++i) {
    
    vector<string> chunks;
    vector<GenomicRegion> spec_info;
    if (i != current_block_index) {
      const size_t offset_in_file = offsets[i];
      current_block = get_block(offset_in_file);
      current_block_index = i;
    }
    current_block.get_slice(intersection(block_index[i], query), 
			    spec_info, chunks);
    
    for (size_t j = 0; j < chunks.size(); ++j) {
      const string spec_name(spec_info[j].get_name());
      OrderMap::iterator k = spec_order.find(spec_name);
      if (k == spec_order.end()) {
	spec_order[spec_name] = segments.size();
	segments.push_back(string(ref_len, gap_symbol));
	segments.back() += chunks[j];
	regions.push_back(spec_info[j]);
      }
      else {
	const size_t spec_id = k->second;
	segments[spec_id].append(ref_len - segments[spec_id].length(), gap_symbol);
	segments[spec_id].append(chunks[j]);
	if (regions[spec_id].get_chrom() == spec_info[j].get_chrom()) {
	  regions[spec_id].set_start(std::min(regions[spec_id].get_start(),
					      spec_info[j].get_start()));
	  regions[spec_id].set_end(std::max(regions[spec_id].get_end(),
					    spec_info[j].get_end()));
	}
      }
    }
    ref_len = segments.front().length();
  }
  for (OrderMap::iterator i = spec_order.begin(); i != spec_order.end(); ++i)
    segments[i->second].append(ref_len - segments[i->second].length(), gap_symbol);
}


std::string
GenomeAlignmentOnDisk::get_sequence(const GenomicRegion& query) const {
  if (!contains(query)) return string();
  size_t curr_index = containing_block(GenomicRegion(query.get_chrom(),
						     query.get_start(),
						     (query.get_start() + 1)));
  const size_t end_index = containing_block(GenomicRegion(query.get_chrom(),
							  (query.get_end() - 1),
							  query.get_end()));
  string r;
  for (size_t i = curr_index; i <= end_index; ++i) {
    if (i != current_block_index) {
      const size_t offset_in_file = offsets[i];
      current_block = get_block(offset_in_file);
      current_block_index = i;
    }
    r += current_block.get_sequence(intersection(block_index[i], query));
  }
  return r;
}
