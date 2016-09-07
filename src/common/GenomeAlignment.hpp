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

#ifndef GENOME_ALIGNMENT_HPP
#define GENOME_ALIGNMENT_HPP

#include "cread.hpp"
#include "GenomicRegion.hpp"

/**
   An AlignmentBlock represents a traditional, order-preserving
   alignment. Alignments of multiple genomes are made of sets of these
   blocks.
*/
class AlignmentBlock {
  std::vector<GenomicRegion> src;
  std::vector<std::string> aln;
  
  static size_t count_add_gaps(const std::string& sequence, size_t offset);
  static size_t count_sub_gaps(const std::string& sequence, size_t offset);
  
public:
  explicit AlignmentBlock() {}
  explicit AlignmentBlock(const std::vector<GenomicRegion> &regions,
			  const std::vector<std::string> &seqs) :
    src(regions), aln(seqs) {}
  
  std::string tostring() const;
  std::string maf_block_string() const;
  
  size_t size() const {return src.size();}
  
  GenomicRegion get_region() const {return src.front();}
  std::vector<GenomicRegion> get_regions() const {return src;}
  std::vector<std::string> get_alignment() const {return aln;}
  
  bool contains(const GenomicRegion &query) const;
  void get_slice(const GenomicRegion &query,
		 std::vector<std::string> &segments) const;
  void get_slice(const GenomicRegion &query,
		 std::vector<GenomicRegion> &species_info,
		 std::vector<std::string> &alnseqs) const;
  std::string get_sequence(const GenomicRegion &query) const;

  static const char gap_symbol = '-';
};

std::ostream& operator<<(std::ostream& s, const AlignmentBlock& ab);

/**
   Each GenomeAlignment is essentially a set of AlignmentBlocks, with
   an index to help in searching the blocks.
*/
class GenomeAlignment {
  static const size_t input_buffer_size = 100000;
  
  std::vector<AlignmentBlock> blocks;
  std::vector<GenomicRegion> block_index;
  
  // size_t containing_block(const GenomicRegion &query) const;
  
public:
  explicit 
  GenomeAlignment(const std::vector<std::vector<GenomicRegion> > &regions,
		  const std::vector<std::vector<std::string> > &seqs);
  explicit 
  GenomeAlignment(const std::string filename);
  
  std::string tostring() const;
  
  bool contains(const GenomicRegion &query) const;
  void get_slice(const GenomicRegion& query,
		 std::vector<std::string> &segments) const;
  void get_slice(const GenomicRegion& query,
		 std::vector<GenomicRegion> &species_info,
		 std::vector<std::string> &segments) const;
  
  std::string get_sequence(const GenomicRegion& query) const;
  
  typedef std::vector<AlignmentBlock>::const_iterator const_iterator;
  const_iterator begin() const {return blocks.begin();}
  const_iterator end() const {return blocks.end();}

  size_t containing_block(const GenomicRegion &query) const;
  
  static const char gap_symbol = '-';
  static const size_t not_found = static_cast<size_t>(-1);
};

class GenomeAlignmentException : public CREADException {
public:
  GenomeAlignmentException(std::string m = "") : CREADException(m) {}
}; 

/**
   Each GenomeAlignment is essentially a set of AlignmentBlocks, with
   an index to help in searching the blocks.
*/
class GenomeAlignmentOnDisk {
  static const size_t input_buffer_size = 100000;

  mutable AlignmentBlock current_block;
  mutable size_t current_block_index;
  mutable std::ifstream alnfile;
  std::vector<GenomicRegion> block_index;
  std::vector<size_t> offsets;
  std::string alnfile_name;
  // size_t containing_block(const GenomicRegion &query) const;
  
public:
  explicit
  GenomeAlignmentOnDisk(const std::string index_filename, 
			const std::string aln_filename);
  
  std::string tostring() const;
  
  bool contains(const GenomicRegion &query) const;
  void get_slice(const GenomicRegion& query,
		 std::vector<std::string> &segments) const;
  void get_slice(const GenomicRegion& query,
		 std::vector<GenomicRegion> &species_info,
		 std::vector<std::string> &segments) const;
  
  std::string get_sequence(const GenomicRegion& query) const;
  
  //   typedef std::vector<AlignmentBlock>::const_iterator const_iterator;
  //   const_iterator begin() const {return blocks.begin();}
  //   const_iterator end() const {return blocks.end();}
  
  size_t containing_block(const GenomicRegion &query) const;
  AlignmentBlock get_block(const size_t offset_in_file) const;
  
  static const char gap_symbol = '-';
  static const size_t not_found = static_cast<size_t>(-1);

  static void ReadIndexFile(const std::string index_filename,
			    std::vector<size_t> &offsets);
};

void
MakeGenomeAlignmentIndex(const std::string alignment_file_name,
			 const std::string index_file_name);

std::ostream& operator<<(std::ostream& s, const GenomeAlignment& mafa);
std::ostream& operator<<(std::ostream& s, const GenomeAlignmentOnDisk& mafa);

class GenomeAlignmentOnDiskException : public CREADException {
public:
  GenomeAlignmentOnDiskException(std::string m = "") : CREADException(m) {}
}; 

#endif
