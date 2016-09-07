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

#include "SuffixTree.hpp"
#include <cassert>

using std::min;
using std::max;
using std::sort;
using std::fill;
using std::greater;
using std::pair;
using std::ptr_fun;
using std::priority_queue;
using std::vector;
using std::string;
using std::numeric_limits;
using std::copy;
using std::fill;
using std::ostringstream;
using std::ostream_iterator;
using std::endl;
using std::bind1st;
using std::make_pair;

using std::cerr;
using std::endl;

class SuffixNode {
public:

  SuffixNode(size_t s = 0, 
	     size_t e = 0, 
	     size_t pl = 0, 
	     size_t c = 0,
	     SuffixNode* p = 0, 
	     SuffixNode* l = 0, 
	     SuffixNode** ch = 0);
  ~SuffixNode();
  
  void BuildSuffixTree(std::string&, size_t);
  void CollectIndices();
  
  std::string tostring(const char*, size_t) const;
  void Insert(char*, size_t, size_t, size_t, size_t);
  
  static int b2i(char b) {
    static int btoi[20] = {
      //A, b, C, d, e, f, g, h, i, j, k, l, m, N, o, p, q, r, s, T
      0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1, 4,-1,-1,-1,-1,-1, 3 
    };
    return (b > 'T' || b < 'A') ? -1 : btoi[b - 'A'];
  }
  
  // Get the top score and index
  void top_score(const char*, SuffixNode **top_node, float**, size_t, 
		 size_t, float, float&) const;
  // Get the top score and index for specified sequence
  void top_score(const char*, SuffixNode **top_node, float**, size_t, 
		 size_t, float, float&, size_t, size_t) const;
  // Get top scores
  void top_scores(const char*, hit_queue&, size_t&, size_t, float**, 
		  size_t, size_t, float, float&) const;
  // Get top scores specifying sequence
  void top_scores(const char*, hit_queue&, size_t&, size_t, float**, 
		  size_t, size_t, float, float&, size_t, size_t) const;
  // Get scores greater than
  void scores_greater(const char*, std::vector<valnode>&, float, float**, size_t, 
		      size_t, float) const;
  // Get scores greater than, specifying sequence
  void scores_greater(const char*, std::vector<valnode>&, float, float**, size_t, 
		      size_t, float, size_t, size_t) const;
  void window_scores_greater(const char*, std::vector<valnode>&, float, float**,
			     size_t, size_t, float, size_t, size_t) const;

  
  /*************** GET REPEATS ****************/
  void get_kmer_counts(std::vector<std::pair<std::string, size_t> > &,
		       size_t, size_t, const char *, std::string&) const;

  size_t get_count(const char *w, const char *text) const;
  void get_locations(const char *w, const char *text, 
		     std::vector<size_t> &locations) const;

  friend class SuffixTree;

private:

  // instance variables
  size_t start;
  size_t end;
  size_t path_length;
  size_t count;
  SuffixNode* parent;
  std::vector<size_t> id;
  SuffixNode* link;
  SuffixNode** child;
  
  void allocate_children();
  
  size_t WalkDown(const char* text, size_t offset) const;
  SuffixNode* JumpUp(int &total_characters_jumped);
  SuffixNode* JumpDown(const char *text,
		       size_t position_in_text,
		       int &characters_above_split_point);

  SuffixNode* Split(SuffixNode*, size_t, const char*, 
		    size_t, size_t, size_t, size_t);
  SuffixNode* InsertUnder(const char* the_text,
			  size_t start_param, 
			  size_t end_param, 
			  size_t depth,
			  size_t index);

  size_t MatchBranch(const char *text, float **scoring_matrix, size_t width,
		     size_t depth, float score, float best_so_far,
		     float &temp_score) const;
  bool HasSequence(size_t min_seqid, size_t max_seqid) const;

  
  /*!
    \brief The number of bases represented by the edge above this node.
  */
  size_t edge_length() const {
    return end - start;
  }

  /*!
    \brief True iff this node has at least one child.
  */
  bool has_children() const {
    return (child != 0);
  }

  /*!
    \brief True iff this node has a child for a particular base.
  */
  bool has_child(const size_t child_index) const {
    return (child[child_index] != 0);
  }
  
  /*!
    \brief True iff the link out of this node has been set.
  */
  bool has_link() const {
    return (link != 0);
  }
};


inline size_t 
SuffixNode::MatchBranch(const char *text, float **scoring_matrix, size_t width,
			size_t depth, float score, float cutoff,
			float &temp_score) const {
  size_t j = 0;
  temp_score = score;
  const size_t limit = std::min(end - start, width - depth);
  const char *text_i = text + start;
  while (j < limit && temp_score > cutoff) {
    temp_score -= scoring_matrix[depth + j][b2i(*text_i)];
    text_i++;
    j++;
  }
  return j;
}

inline bool
SuffixNode::HasSequence(size_t min_seqid, size_t max_seqid) const {
  std::vector<size_t>::const_iterator k = 
    lower_bound(id.begin(), id.end(), min_seqid);
  return (k != id.end() && *k < max_seqid);
}


/****************************************************
 * ALLOCATE SPACE FOR THE CHILDREN OF A SUFFIX NODE
 ****************************************************/
void
SuffixNode::allocate_children() {
  child = new SuffixNode *[alphabet_size + 1];
  fill(child, child + alphabet_size + 1, static_cast<SuffixNode *>(0));
}

/***************
 * constructor
 ***************/
SuffixNode::SuffixNode(size_t s, size_t e, size_t pl, size_t c,
		       SuffixNode *p, SuffixNode *l, SuffixNode **ch) :
  start(s), end(e), path_length(pl), count(c), parent(p), link(l), child(ch) {}

/**************
 * destructor
 **************/
SuffixNode::~SuffixNode() {
  if (has_children()) {
    for (size_t i = 0; i < alphabet_size + 1; ++i) 
      if (has_child(i))
	delete child[i];
    delete[] child;
  }
}

string
SuffixNode::tostring(const char *text, size_t depth) const {
  ostringstream s;
  fill_n(ostream_iterator<char>(s), depth, ' ');
  s << "(" << start << "," << end << ") ";
  copy(text + start, text + end, ostream_iterator<char>(s));
  s << " [";
  copy(id.begin(), id.end(), ostream_iterator<size_t>(s, " "));
  s << "]";
  s << endl;
  for (size_t i = 0; has_children() && i < alphabet_size + 1; i++)
    if (has_child(i)) 
      s << child[i]->tostring(text, depth + 1);
  return s.str();
}

/**********************************************************************
 *  WalkDown returns the number of characters that are to be above the
 *  split.  Finds how far down the edge to split.  Compares characters
 *  as it goes.
 **********************************************************************/
size_t
SuffixNode::WalkDown(const char *text, const size_t offset) const {
  const size_t length = edge_length();
  size_t i = 0;
  while (i < length && (text + offset)[i] == (text + start)[i])
    i++;
  return i;
}

/*********************************************************************
 * JumpDown already knows how far down to travel before a split must
 * occur. It moves down from SuffNode to SuffNode, only comparing to
 * make sure the path is correct (one comparison per SuffNode).  This
 * implements the "skip-count trick" from Gusfield.
 *********************************************************************/
SuffixNode *
SuffixNode::JumpDown(const char *text, 
		     size_t position_in_text,
		     int &characters_above_split_point) {
  SuffixNode *current_node = this;
  size_t next_child_index = b2i(text[position_in_text]);
  while (current_node->has_children() && 
	 current_node->has_child(next_child_index) &&
	 (characters_above_split_point >
	  static_cast<int>(current_node->child[next_child_index]->edge_length()))) {
    current_node = current_node->child[next_child_index];
    position_in_text += current_node->edge_length();
    characters_above_split_point -= current_node->edge_length();
    next_child_index = b2i(text[position_in_text]);
  }
  return current_node;
}

/*************************************************************************
 * JumpUp finds and returns the lowest ancestor of the calling
 * SuffixNode that has a suffix link out from it.  The number of
 * characters on the path up is recorded so that they won't have to be
 * matched when the next insertion occurs that many characters below
 * the SuffixNode at the end of the link out of the SuffixNode where
 * the "jumping up" ends.
 *************************************************************************/ 
SuffixNode *
SuffixNode::JumpUp(int &total_characters_jumped) {
  SuffixNode *current_node = this;
  // while the current node has no long, jump to its parent
  while (!current_node->has_link()) {
    // record number of characters jumped
    total_characters_jumped += current_node->edge_length();
    current_node = current_node->parent;
  }
  return current_node;
}

/***********************
 * Split the edge above the SuffNode "old_child"
 * by adding a new SuffNode "new_child" above it
 * and a new leaf "new_leaf" below "new_child".
 ***************************/
SuffixNode *
SuffixNode::Split(SuffixNode *old_child, size_t depth, const char *text,
		  size_t s, size_t e, size_t parent_tpl, size_t index) {
  // Create the new nodes
  SuffixNode *new_child = new SuffixNode(old_child->start,
					 old_child->start + depth,
					 parent_tpl + depth, count, this);
  SuffixNode *new_leaf = new SuffixNode(s + depth, e, parent_tpl + e - s,
					1, new_child);
  new_leaf->id.push_back(index);
  // Modify the old child
  old_child->parent = new_child;
  old_child->start += depth;
  // allocate the child pointers for the new node and set them
  new_child->allocate_children();
  new_child->child[b2i(text[old_child->start])] = old_child;
  new_child->child[b2i(text[new_leaf->start])] = new_leaf;
  return new_child;
}

/*********************
 * All insertions are under a particular SuffNode,
 * even if it is just under the root
 * (which is given a string not in the alphabet).
 ******************************/
SuffixNode *
SuffixNode::InsertUnder(const char *text, size_t s, size_t e, 
			size_t depth, size_t index) {
  int i = b2i(text[s]);
  if (i == -1) return this;
  if (!has_children())
    allocate_children();
  if (!has_child(i)) {
    child[i] = new SuffixNode(s, e, path_length + e - s, 1, this);
    child[i]->id.push_back(index);
  }
  else {
    // !!! The following comment (referring to (depth == 0) condition)
    // might be wrong: 
    // --> This means that the previous "jumpUp" landed at the root <--
    if (depth == 0) 
      depth = child[i]->WalkDown(text, s); // So a "walkDown" is in order.
    if (depth < e - s) {
      if (depth < child[i]->edge_length()) {
	child[i] = Split(child[i], depth, text, s, e, path_length, index);
	/* This is commented out below so that all the ids are only
	   at the leaves */
	// child[i]->id.push_back(index);
	return child[i];
      }
      else {
	return child[i]->InsertUnder(text, s + child[i]->edge_length(), e,
				     depth - child[i]->edge_length(), index);
      }
    }
    // This next line is here for bounded depth trees, where if a
    // substring is inserted that is identical to the current
    // substring, then that substring will end at the exact same spot
    // as the previous, and require the 'index' to be added
    else if (depth == child[i]->edge_length()) 
      child[i]->id.push_back(index);
  }
  return this;
}

/*************************
 * This is the interface to building a suffix tree.
 ************************************/
void
SuffixNode::BuildSuffixTree(string &text, size_t max_depth) {
  
  // set the maximum (string) depth of the tree
  max_depth = min(max_depth, text.length());

  SuffixNode *current_node = this->link = this;
  SuffixNode *previous_node = 0;
  
  // loop to insert each suffix
  for (size_t i = 0; i < text.length(); i++) {
    // jump up until a node with a link is found (k keeps track of how
    // many nodes were jumped).
    int distance_above_insertion = 0; // Keeps track of nodes jumped up
    current_node = current_node->JumpUp(distance_above_insertion);
    // if the current_node is the root... 
    if (current_node == current_node->link) {
      assert(this == current_node && 
	     "current_node has self link, but is not root!");
      // this clause should only be satisfied 3 times (alphabet_size
      // - 1) because it corresponds to insertion of a suffix that
      // starts with a completely new symbol.
      if (distance_above_insertion > 0)
	--distance_above_insertion;
    }
    // Jump down the same number of nodes we jumped up.
    current_node = 
      current_node->link->JumpDown(text.c_str(), 
				   i + current_node->link->path_length, 
				   distance_above_insertion);
    // We have found the parent of the next insertion point, so
    // insert under it.
    const size_t start = i + current_node->path_length;
    const size_t end = min(text.length(), i + max_depth);
    current_node = current_node->InsertUnder(text.c_str(), start, end,
					     distance_above_insertion, i);
    // If a link should be made, make it.
    if (previous_node &&
	previous_node->path_length == current_node->path_length + 1)
      previous_node->link = current_node;
    // make the previous node point to the current_node node in
    // preparation for the next iteration.
    previous_node = current_node;
  }
}

void
SuffixNode::CollectIndices() {
  if (child) {
    // Need to iterate over _ALL_ characters (even the 'N') because
    // otherwise some indexes won't get counted. This would happen for
    // every suffix that completely matches some prefix of an earlier
    // suffix, and therefore has a split where one of the children is
    // just an 'N'.
    for (size_t i = 0; i <= alphabet_size; ++i)
      if (child[i]) {
	child[i]->CollectIndices();
	const size_t previous_end = id.size();
	id.insert(id.end(), child[i]->id.begin(), child[i]->id.end());
	std::inplace_merge(id.begin(), id.begin() + previous_end, id.end());
      }
  }
  /* The commented-out code below should not be needed, because all
     the ids should be at the leaves, and no id should be in more than
     one leaf. */
  // id.erase(std::unique(id.begin(), id.end()), id.end());
}

/******************* SUBSTRING NODE QUERY METHODS ***********/

void
SuffixNode::top_score(const char *text, SuffixNode **top_node,
		      float **scoremat, size_t width, size_t depth,
		      float score, float &cutoff) const {
  if (child){
    for (size_t i = 0; i < alphabet_size; i++){
      if (child[i]) {
	float temp_score = 0;
	size_t j = child[i]->MatchBranch(text, scoremat, width, depth, score, 
					 cutoff, temp_score);
	if (temp_score > cutoff){
	  if (depth + j < width){
	    child[i]->top_score(text, top_node, scoremat, width, 
				depth + j, temp_score, cutoff);
          } 
          
	  else if (temp_score > cutoff) {
	    *top_node = child[i];
	    cutoff = temp_score;
          }
        }
      }   
    }
  }	  
}
void
SuffixNode::top_score(const char *text, SuffixNode **top_node, 
		      float **scoremat, size_t width, size_t depth, 
		      float score, float &cutoff,
		      size_t min_seqid, size_t max_seqid) const {
  if (child){
    for (size_t i = 0; i < alphabet_size; i++){
      if (child[i]) {
	float temp_score = 0;
	size_t j = child[i]->MatchBranch(text, scoremat, width, depth, score, 
					 cutoff, temp_score);
	if (temp_score > cutoff && child[i]->HasSequence(min_seqid, max_seqid)){
	  if (depth + j < width){
	    child[i]->top_score(text, top_node, scoremat, width, depth + j,
				temp_score, cutoff, min_seqid, max_seqid);
          }
	  else if (temp_score > cutoff) {
	    *top_node = child[i];
	    cutoff = temp_score;
	  }
        }
      }
    }
  }
}

// GET SCORE GREATER THAN A THRESHOLD
void
SuffixNode::scores_greater(const char *text, vector<valnode> &V, float cutoff, 
			   float **scoremat, size_t width, size_t depth, 
			   float score) const {
  if (child){
    for (size_t i = 0; i < alphabet_size; i++) {
      if (child[i]) {
	float temp_score = 0;
	size_t j = child[i]->MatchBranch(text, scoremat, width, depth, score, 
					 cutoff, temp_score);
	if (temp_score >= cutoff){
	  if (depth + j < width)
	    child[i]->scores_greater(text, V, cutoff, scoremat,
				     width, depth + j, temp_score);
	  else V.push_back(valnode(temp_score, child[i]));
        }
      }
    }
  }
}

// Get scores greater than, specifying sequence
void
SuffixNode::scores_greater(const char *text, vector<valnode> &V, float cutoff, 
			   float **scoremat, size_t width, size_t depth, 
			   float score, size_t min_seqid, size_t max_seqid) const {
  if (child){
    for (size_t i = 0; i < alphabet_size; i++){
      if (child[i]) {
	float temp_score = 0;
	size_t j = child[i]->MatchBranch(text, scoremat, width, depth, score, 
					 cutoff, temp_score);
	if (temp_score >= cutoff && child[i]->HasSequence(min_seqid, max_seqid)){
	  if (depth + j < width) 
	    child[i]->scores_greater(text, V, cutoff, scoremat, width, 
				     depth + j, temp_score, min_seqid, max_seqid);
	  else V.push_back(valnode(temp_score, child[i]));
        }
      }
    }
  }
}

// TODO: record the size of a hit_queue in terms of the size of the nodes
// GET TOP SCORES
void 
SuffixNode::top_scores(const char *text, hit_queue &PQ, size_t &pqsize, size_t n_top, 
		       float **scoremat, size_t width, size_t depth, 
		       float score, float &cutoff) const {
  if (child){
    for (size_t i = 0; i < alphabet_size; i++){
      if (child[i]) {
	float temp_score = 0;
	size_t j = child[i]->MatchBranch(text, scoremat, width, depth, score, 
					 cutoff, temp_score);
	if (temp_score >= cutoff){
	  if (depth + j < width)
	    child[i]->top_scores(text, PQ, pqsize, n_top, scoremat, width,
				 depth + j, temp_score, cutoff);
	  else {
	    PQ.push(valnode(temp_score, child[i]));
	    pqsize += child[i]->id.size();
	    while (pqsize - PQ.top().second->id.size() >= n_top) {
	      pqsize -= PQ.top().second->id.size();
	      PQ.pop();
	    }
	    if (PQ.size() >= n_top) cutoff = PQ.top().first;
	  }
        }
      }
    }
  }
}

void 
SuffixNode::top_scores(const char *text, hit_queue &PQ, size_t &pqsize, size_t n_top, 
		       float **scoremat, size_t width, size_t depth, float score,
		       float &cutoff, size_t min_seqid, size_t max_seqid) const {
  if (child){
    for (size_t i = 0; i < alphabet_size; i++){
      if (child[i]) {
	float temp_score = 0;
	size_t j = child[i]->MatchBranch(text, scoremat, width, depth, score, 
					 cutoff, temp_score);
	if (temp_score >= cutoff && child[i]->HasSequence(min_seqid, max_seqid)){
	  if (depth + j < width)
	    child[i]->top_scores(text, PQ, pqsize, n_top, scoremat, width, depth + j, 
				 temp_score, cutoff, min_seqid, max_seqid);
	  else {
	    PQ.push(valnode(temp_score, child[i]));
	    pqsize += child[i]->id.size();
	    while (pqsize - PQ.top().second->id.size() >= n_top) {
	      pqsize -= PQ.top().second->id.size();
	      PQ.pop();
	    }
	    if (PQ.size() >= n_top) cutoff = PQ.top().first;
	  }
        }
      }
    }
  }
}


void
SuffixNode::get_kmer_counts(vector<pair<string, size_t> > &kmers, size_t kmer,
			    size_t depth, const char *text, 
			    string& current_prefix) const {
  for (size_t i = 0; i < alphabet_size; i++) {
    if (child[i]) {
      const size_t local_start = child[i]->start;
      const size_t limit = min(kmer - depth, child[i]->end - local_start);
      for (size_t j = 0; j < limit; ++j)
	current_prefix[depth + j] = text[local_start + j];
      if (depth + child[i]->edge_length() >= kmer)
	kmers.push_back(make_pair(current_prefix, child[i]->id.size()));
      else if (child[i]->has_children()) {
	child[i]->get_kmer_counts(kmers, kmer, 
				  depth + child[i]->end - local_start,
				  text, current_prefix);
      }
    }
  }
}

size_t
SuffixNode::get_count(const char *word, const char *text) const {
  if (*word == '\0') return id.size();
  SuffixNode *temp;
  if (child && (temp = child[b2i(*word)])) {
    size_t j = 0;
    const size_t lim = temp->edge_length();
    while (j < lim && word[j] == text[temp->start + j]) j++;
    return temp->get_count(word + j, text);
  }
  else return 0;
}

void
SuffixNode::get_locations(const char *word, const char *text,
			  vector<size_t> &locations) const {
  if (*word == '\0')
    copy(id.begin(), id.end(), back_inserter(locations));
  else {
    SuffixNode *temp;
    if (child && (temp = child[b2i(*word)])) {
      size_t j = 0;
      const size_t lim = temp->edge_length();
      while (j < lim && word[j] == text[temp->start + j]) ++j;
      temp->get_locations(word + j, text, locations);
    }
  }
}

/**********************************************************************
 **********************************************************************
 *************************** SUFFIX TREE ******************************
 **********************************************************************
 **********************************************************************/

pair<size_t, size_t>
SuffixTree::index2seq_offset(size_t index) const {
  // TODO: this should be a binary search -- these offsets are in
  // non-decreasing order.
  size_t j = 0;
  while (j < offset.size() && offset[j] <= index) j++;
  return make_pair(j - 1, index - offset[j - 1]);
}

SuffixTree::SuffixTree(const string& s, int d) : sequence(s) {
  sequence.append("N");
  offset.push_back(0);
  offset.push_back(sequence.length() + 1);
  if (d == -1) depth = sequence.length();
  else depth = d;
  root = new SuffixNode();
  root->BuildSuffixTree(sequence, depth);
  root->CollectIndices();
}

SuffixTree::SuffixTree(const vector<string>& s, int d) {
  size_t total = 0;
  vector<string>::const_iterator i(s.begin());
  for (; i != s.end(); ++i) {
    offset.push_back(total);
    sequence.append(*i);
    sequence.append("N");
    total += i->length() + 1;
  }
  offset.push_back(total);
  if (d == -1)
    depth = sequence.length();
  else depth = d;
  root = new SuffixNode();
  root->BuildSuffixTree(sequence, depth);
  root->CollectIndices();
}

/******* FUNCTIONS DEALING WITH QUERIES *******/

float**
SuffixTree::PrepareScoringMatrix(const ScoringMatrix& matrix, 
				 float &max_score) {
  float **scoremat = new float *[matrix.get_width()];
  max_score = 0.0;
  for (size_t i = 0; i < matrix.get_width(); i++) {
    scoremat[i] = new float[alphabet_size + 1];
    float max_column_score = 0.0;
    for (size_t j = 0; j < alphabet_size; j++) {
      scoremat[i][j] = matrix[i][j];
      max_column_score = max(scoremat[i][j], max_column_score);
    }
    for (size_t j = 0; j < alphabet_size; j++)
      scoremat[i][j] = max_column_score - scoremat[i][j];
    scoremat[i][alphabet_size] = numeric_limits<float>::max();
    max_score += max_column_score;
  }
  return scoremat;
}

float
SuffixTree::top_score(const ScoringMatrix& matrix, int seqid) const {
  float max_score = 0;
  float **scoremat = PrepareScoringMatrix(matrix, max_score);
  float r = -numeric_limits<float>::max();
  SuffixNode *top_node = 0;
  if (seqid == -1) 
    root->top_score(sequence.c_str(), &top_node, scoremat, 
		    min(matrix.get_width(), depth), 0, max_score, r);
  else root->top_score(sequence.c_str(), &top_node, scoremat,
		       min(matrix.get_width(), depth), 0, max_score, r,
		       offset[seqid], offset[seqid + 1]);
  for (size_t i = 0; i < matrix.get_width(); i++) 
    delete[] scoremat[i];
  delete[] scoremat;
  return r;
}

float
SuffixTree::top_score(const ScoringMatrix& matrix, 
		      size_t start, size_t end, int seqid) const {
  float max_score = 0;
  float **scoremat = PrepareScoringMatrix(matrix, max_score);
  float r = -numeric_limits<float>::max();
  SuffixNode *top_node = 0;
  if (seqid == -1) 
    root->top_score(sequence.c_str(), &top_node, scoremat,
		    min(matrix.get_width(), depth), 0, max_score, r,
		    start, end);
  else root->top_score(sequence.c_str(), &top_node, scoremat,
		       min(matrix.get_width(), depth), 0, max_score, r,
		       offset[seqid] + start, 
		       offset[seqid] + end);
  for (size_t i = 0; i < matrix.get_width(); i++)
    delete[] scoremat[i];
  delete[] scoremat;
  return r;
}

void 
SuffixTree::top_scores(vector<float> &A, size_t n_top, 
		       const ScoringMatrix& matrix, int seqid) const {

  // TODO: define an exception for this
  if (n_top == 0) return;

  float max_score = 0;
  float **scoremat = PrepareScoringMatrix(matrix, max_score);
  float cutoff = -numeric_limits<float>::max();
  hit_queue PQ;
  size_t pqsize = 0;
  if (seqid == -1) 
    root->top_scores(sequence.c_str(), PQ, pqsize, n_top, scoremat, 
		     min(matrix.get_width(), depth), 0, max_score, cutoff);
  else 
    root->top_scores(sequence.c_str(), PQ, pqsize, n_top, scoremat,
		     min(matrix.get_width(), depth), 0, max_score, cutoff,
		     offset[seqid], offset[seqid + 1]);
  for (size_t i = 0; i < matrix.get_width(); i++) 
    delete[] scoremat[i];
  delete[] scoremat;
  // TODO: make sure this is sorted correctly before popping off
  while (PQ.size() > 0) {
    for (size_t i = 0; i < PQ.top().second->id.size(); ++i) {
      pair<size_t, size_t> ind_off(index2seq_offset(PQ.top().second->id[i]));
      if (seqid == -1 || static_cast<int>(ind_off.first) == seqid)
	A.push_back(PQ.top().first);
    }
    PQ.pop();
  }
  std::reverse(A.begin(), A.end());
  A.erase(A.begin() + n_top, A.end());
}

void
SuffixTree::scores_greater(vector<float> &V, float threshold, 
			   const ScoringMatrix& matrix, int seqid) const {
  float max_score = 0;
  float **scoremat = PrepareScoringMatrix(matrix, max_score);
  vector<valnode> B;
  if (seqid == -1) 
    root->scores_greater(sequence.c_str(), B, threshold, scoremat, 
			 min(matrix.get_width(), depth), 0, max_score);
  else root->scores_greater(sequence.c_str(), B, threshold, scoremat, 
			    min(matrix.get_width(), depth), 0, max_score, 
			    offset[seqid], offset[seqid + 1]);
  for (size_t i = 0; i < matrix.get_width(); i++) 
    delete[] scoremat[i];
  delete[] scoremat;
  for (vector<valnode>::iterator i = B.begin(); i != B.end(); ++i)
    for (size_t j = 0; j < i->second->id.size(); ++j) {
      pair<size_t, size_t> ind_off(index2seq_offset(i->second->id[j]));
      if (seqid == -1 || static_cast<int>(ind_off.first) == seqid)
	V.push_back(i->first);
    }
}

seq_pos_score
SuffixTree::top_score_index(const ScoringMatrix& matrix, int seqid) const {
  float max_score = 0;
  float **scoremat = PrepareScoringMatrix(matrix, max_score);
  float cutoff = -numeric_limits<float>::max();
  SuffixNode *top_node = 0;
  if (seqid == -1) 
    root->top_score(sequence.c_str(), &top_node, scoremat, 
		    matrix.get_width(), 0, max_score, cutoff);
  else 
    root->top_score(sequence.c_str(), &top_node, scoremat, 
		    matrix.get_width(), 0, max_score, cutoff,
		    offset[seqid], offset[seqid + 1]);
  for (size_t i = 0; i < matrix.get_width(); i++) 
    delete[] scoremat[i];
  delete[] scoremat;
  return seq_pos_score(index2seq_offset(top_node->id.front()), cutoff);
}

seq_pos_score
SuffixTree::top_score_index(const ScoringMatrix& matrix, 
			    size_t lower, size_t upper, int seqid) const {
  float max_score = 0;
  float **scoremat = PrepareScoringMatrix(matrix, max_score);
  float cutoff = -numeric_limits<float>::max();
  SuffixNode *top_node = 0;
  if (seqid == -1) 
    root->top_score(sequence.c_str(), &top_node, scoremat, 
		    matrix.get_width(), 0, max_score, cutoff,
		    lower, upper);
  else 
    root->top_score(sequence.c_str(), &top_node, scoremat, 
		    matrix.get_width(), 0, max_score, cutoff,
		    offset[seqid] + lower, offset[seqid] + upper);
  for (size_t i = 0; i < matrix.get_width(); i++)
    delete[] scoremat[i];
  delete[] scoremat;
  return seq_pos_score(index2seq_offset(top_node->id.front()), cutoff);
}

void 
SuffixTree::top_scores_indices(vector<seq_pos_score> &A, size_t n_top, 
			       const ScoringMatrix& matrix, 
			       int seqid) const {
  
  // TODO: define an exception for this
  if (n_top == 0) return;
  
  float max_score = 0;
  float **scoremat = PrepareScoringMatrix(matrix, max_score);
  float cutoff = -numeric_limits<float>::max();
  hit_queue PQ;
  size_t pqsize = 0;
  if (seqid == -1) 
    root->top_scores(sequence.c_str(), PQ, pqsize, n_top, scoremat,
		     matrix.get_width(), 0, max_score, cutoff);
  else
    root->top_scores(sequence.c_str(), PQ, pqsize, n_top, scoremat,
		     matrix.get_width(), 0, max_score, cutoff,
		     offset[seqid], offset[seqid + 1]);
  for (size_t i = 0; i < matrix.get_width(); i++)
    delete[] scoremat[i];
  delete[] scoremat;
  // TODO: make sure this is sorted correctly before popping off
  while (PQ.size() > 0) {
    for (size_t i = 0; i < PQ.top().second->id.size(); ++i) {
      pair<size_t, size_t> ind_off(index2seq_offset(PQ.top().second->id[i]));
      if (seqid == -1 || static_cast<int>(ind_off.first) == seqid)
	A.push_back(seq_pos_score(ind_off, PQ.top().first));
    }
    PQ.pop();
  }
  std::reverse(A.begin(), A.end());
  A.erase(A.begin() + min(n_top, A.size()), A.end());
}

void
SuffixTree::scores_greater_indices(float threshold,
				   const ScoringMatrix& matrix, 
				   vector<seq_pos_score> &A,
				   int seqid) const {
  float max_score = 0;
  float **scoremat = PrepareScoringMatrix(matrix, max_score);
  vector<valnode> nodes_above_threshold;
  if (seqid == -1) 
    root->scores_greater(sequence.c_str(), nodes_above_threshold, 
			 threshold, scoremat, matrix.get_width(), 0, max_score);
  else 
    root->scores_greater(sequence.c_str(), nodes_above_threshold,
			 threshold, scoremat, matrix.get_width(), 0, max_score,
			 offset[seqid], offset[seqid + 1]);
  for (size_t i = 0; i < matrix.get_width(); i++) 
    delete[] scoremat[i];
  delete[] scoremat;
  // TODO: make sure nodes_above_threshold is sorted properly
  sort(nodes_above_threshold.begin(), nodes_above_threshold.end());
  for (vector<valnode>::iterator i = nodes_above_threshold.begin(); 
       i != nodes_above_threshold.end(); ++i)
    for (size_t j = 0; j < i->second->id.size(); ++j) {
      pair<size_t, size_t> ind_off(index2seq_offset(i->second->id[j]));
      if (seqid == -1 || static_cast<int>(ind_off.first) == seqid)
	A.push_back(seq_pos_score(ind_off, i->first));
    }
}

void
SuffixTree::window_scores_greater_indices(std::vector<seq_pos_score> &A, float threshold,
                                          const ScoringMatrix& matrix, 
					  size_t start, size_t stop, size_t seqid) const{
  float max_score = 0;
  float **scoremat = PrepareScoringMatrix(matrix, max_score);
  vector<valnode> B;
  root->scores_greater(sequence.c_str(), B, threshold, scoremat,
		       matrix.get_width(), 0, max_score,
		       start, stop);
  for (size_t i = 0; i < matrix.get_width(); ++i)
    delete[] scoremat[i];
  delete[] scoremat;
  
  sort(B.begin(), B.end());
  for (vector<valnode>::iterator i = B.begin(); i != B.end(); ++i){
    for (size_t j = 0; j < i->second->id.size(); ++j){
      // check to make sure in window
      if (i->second->id[j] > start && i->second->id[j] < stop)
	A.push_back(seq_pos_score(pair<size_t, size_t>(seqid, i->second->id[j]), i->first));
    }
  }
}

bool
SuffixTree::no_overlaps(vector<pair<size_t, size_t> >& sites,
			size_t p, size_t w) {
  for (size_t c = 0; c < sites.size(); ++c)
    if ((p >= sites[c].first && p <= sites[c].second) || 
	(p + w >=sites[c].first && p + w <= sites[c].second) ||
	(p <= sites[c].first && p+w >= sites[c].second))
      return false;
  return true;
}

void
SuffixTree::window_scores_greater_indices(std::vector<seq_pos_score> &A, 
					  float threshold,
                                          const ScoringMatrix& matrix,
					  size_t start, 
					  size_t stop, 
					  size_t seqid,
					  vector<pair<size_t, size_t> > &recorded_sites) const{
  float max_score = 0;
  float **scoremat = PrepareScoringMatrix(matrix, max_score);
  vector<valnode> B;
  root->scores_greater(sequence.c_str(), B, threshold, scoremat,
		       matrix.get_width(), 0, max_score, start, stop);
  for (size_t i = 0; i < matrix.get_width(); ++i)
    delete[] scoremat[i];
  delete[] scoremat;
  
  sort(B.begin(), B.end());
  for (vector<valnode>::iterator i = B.begin(); i != B.end(); ++i)
    for (size_t j = 0; j < i->second->id.size(); ++j)
      // check to make sure doesn't overlap sites in recorded_sites
      if (no_overlaps(recorded_sites, i->second->id[j], matrix.get_width()))
	// check to make sure in window
	if (i->second->id[j] > start && i->second->id[j] < stop)
	  A.push_back(seq_pos_score(make_pair(seqid, i->second->id[j]), i->first));
}


string
SuffixTree::tostring() const {
  return sequence + string("\n") + root->tostring(sequence.c_str(), 0);
}

void
SuffixTree::get_kmer_counts(vector<pair<string, size_t> > &kmers,
			    size_t kmer) const {
  string s(kmer, 'N');
  root->get_kmer_counts(kmers, kmer, 0, sequence.c_str(), s);
}

void
SuffixTree::get_locations(string w, vector<pair<size_t, size_t> > &L) const {
  vector<size_t> locations;
  root->get_locations(w.c_str(), sequence.c_str(), locations);
  transform(locations.begin(), locations.end(), back_inserter(L),
	    bind1st(mem_fun(&SuffixTree::index2seq_offset), this));
}
// TO DO: CHECK FOR MATRIX WIDTH

size_t
SuffixTree::get_count(string w) const {
  return root->get_count(w.c_str(), sequence.c_str());
}

SuffixTree::~SuffixTree() {
  delete root;
}
