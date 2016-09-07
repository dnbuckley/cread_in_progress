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

#include "IndexedSuffixTree.hpp"
#include <cassert>

using std::vector;
using std::string;

static const size_t alpha_plus_two = alphabet_size + 2;


// class IndexedSuffixNode;

typedef std::pair<float, IndexedSuffixNode*> valnode;
typedef std::priority_queue<valnode, std::vector<valnode>, 
			    std::greater<valnode> > hit_queue;


class IndexedSuffixNode {

  size_t start;
  size_t end;
  size_t path_length;
  size_t count;
  IndexedSuffixNode* parent;
  IndexedSuffixNode* link;
  IndexedSuffixNode** child;
  
  void allocate_children();
  
  size_t WalkDown(const char* text, size_t offset) const;
  IndexedSuffixNode* JumpUp(int &total_characters_jumped);
  IndexedSuffixNode* JumpDown(const char *text,
			      size_t position_in_text,
			      int &characters_above_split_point);
  
  IndexedSuffixNode* Split(IndexedSuffixNode*, size_t, const char*, 
			   size_t, size_t, size_t, size_t);
  IndexedSuffixNode* InsertUnder(const char* the_text,
				 size_t start_param, 
				 size_t end_param, 
				 size_t depth,
				 size_t index);

  void thread_leaves(IndexedSuffixNode** prev_leaf);
  void adjust_links();
  void collect_leaf_indices(std::vector<size_t>& indices) const;
  
  bool HasSequence(size_t min_seqid, size_t max_seqid) const;
  
  size_t edge_length() const {
    return end - start;
  }
  bool has_children() const {
    return (child != 0);
  }
  bool is_leaf() const {
    return (child == 0);
  }
  bool has_child(const size_t child_index) const {
    return (child[child_index] != 0);
  }
  bool has_link() const {
    return (link != 0);
  }

  // Get scores greater than
  void scores_greater(const char*, std::vector<valnode>&, float, float**, 
		      size_t, size_t, float) const;
  void top_scores(const char*, hit_queue&, size_t&, size_t, float**,
		  size_t, size_t, float, float&) const;
  
  
  size_t MatchBranch(const char *text, float **scoring_matrix, size_t width,
		     size_t depth, float score, float best_so_far,
		     float &temp_score) const;
  
public:

  IndexedSuffixNode(size_t = 0, size_t = 0, size_t = 0, size_t = 0,
		    IndexedSuffixNode* = 0, IndexedSuffixNode* = 0, IndexedSuffixNode** = 0);
  ~IndexedSuffixNode();
  
  void BuildSuffixTree(std::string&);

  float top_score(const ScoringMatrix&) const;
  
  std::string tostring(const char*, size_t) const;
  void Insert(char*, size_t, size_t, size_t, size_t);
  
  static char Int2Base(size_t i) {
    static char itob[4] = { 'A', 'C', 'G', 'T' };
    return (i > 3) ? 'N' : itob[i];
    //original line = "return (i < 0 || i > 3) ? 'N' : itob[i];" 
  }
  static int b2i(char b) {
    static int btoi[20] = {
      //A, b, C, d, e, f, g, h, i, j, k, l, m, N, o, p, q, r, s, T
      0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1, 4,-1,-1,-1,-1,-1, 3 
    };
    if (b == '$') return alphabet_size + 1;
    return (b > 'T' || b < 'A') ? -1 : btoi[b - 'A'];
  }
  static bool valid(char b) {return b2i(b) != -1;}

  friend class IndexedSuffixTree;
};



inline size_t 
IndexedSuffixNode::MatchBranch(const char *text, float **scoring_matrix, size_t width,
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


/****************************************************
 * ALLOCATE SPACE FOR THE CHILDREN OF A SUFFIX NODE
 ****************************************************/
void
IndexedSuffixNode::allocate_children() {
  child = new IndexedSuffixNode *[alpha_plus_two];
  std::fill(child, child + alpha_plus_two, static_cast<IndexedSuffixNode *>(0));
}


/***************
 * constructor
 ***************/
IndexedSuffixNode::IndexedSuffixNode(size_t s, size_t e, size_t pl, size_t c,
				     IndexedSuffixNode *p, IndexedSuffixNode *l, IndexedSuffixNode **ch) :
  start(s), end(e), path_length(pl), count(c), parent(p), link(l), child(ch) {
}


/**************
 * destructor
 **************/
IndexedSuffixNode::~IndexedSuffixNode() {
  if (has_children()) {
    for (size_t i = 0; i < alpha_plus_two; ++i) 
      if (has_child(i))
	delete child[i];
    delete[] child;
  }
}


string
IndexedSuffixNode::tostring(const char *text, size_t depth) const {
  std::ostringstream s;
  fill_n(std::ostream_iterator<char>(s), depth, ' ');
  s << "<" << this << "> ";
  s << "<" << parent << "> ";
  s << "<" << link << "> ";
  s << "(" << start << "," << end << ") __";
  copy(text + start, text + end, std::ostream_iterator<char>(s));
  s << "__ [" << path_length << "]";
  if (has_children()) {
    vector<size_t> indices;
    collect_leaf_indices(indices);
    copy(indices.begin(), indices.end(),
	 std::ostream_iterator<size_t>(s, ","));
  }
  else s << end;
  s << "\t" << "|" << count << "|" << std::endl;
  if (has_children()) {
    for (size_t i = 0; i < alpha_plus_two; i++)
      if (has_child(i)) 
	s << child[i]->tostring(text, depth + 1);
  }
  return s.str();
}


/**********************************************************************
 *  WalkDown returns the number of characters that are to be above the
 *  split.  Finds how far down the edge to split.  Compares characters
 *  as it goes.
 **********************************************************************/
size_t
IndexedSuffixNode::WalkDown(const char *text, const size_t offset) const {
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
IndexedSuffixNode *
IndexedSuffixNode::JumpDown(const char *text, 
			    size_t position_in_text,
			    int &characters_above_split_point) {
  IndexedSuffixNode *curr_node = this;
  size_t next_child_index = b2i(text[position_in_text]);
  while (curr_node->has_children() && 
	 curr_node->has_child(next_child_index) &&
	 (characters_above_split_point >
	  static_cast<int>(curr_node->child[next_child_index]->edge_length()))) {
    curr_node = curr_node->child[next_child_index];
    position_in_text += curr_node->edge_length();
    characters_above_split_point -= curr_node->edge_length();
    next_child_index = b2i(text[position_in_text]);
  }
  return curr_node;
}


/*************************************************************************
 * JumpUp finds and returns the lowest ancestor of the calling
 * IndexedSuffixNode that has a suffix link out from it.  The number
 * of characters on the path up is recorded so that they won't have to
 * be matched when the next insertion occurs that many characters
 * below the IndexedSuffixNode at the end of the link out of the
 * IndexedSuffixNode where the "jumping up" ends.
 *************************************************************************/ 
IndexedSuffixNode *
IndexedSuffixNode::JumpUp(int &total_characters_jumped) {
  IndexedSuffixNode *current_node = this;
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
IndexedSuffixNode *
IndexedSuffixNode::Split(IndexedSuffixNode *old_child, size_t depth, const char *text,
			 size_t s, size_t e, size_t parent_tpl, size_t index) {
  // Create the new nodes
  IndexedSuffixNode *new_child = new IndexedSuffixNode(old_child->start,
						       old_child->start + depth,
						       parent_tpl + depth, count, this);
  IndexedSuffixNode *new_leaf = new IndexedSuffixNode(s + depth, e, parent_tpl + e - s,
						      1, new_child);
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
IndexedSuffixNode *
IndexedSuffixNode::InsertUnder(const char *text, size_t s, size_t e,
			       size_t depth, size_t index) {
  int i = b2i(text[s]);
  if (i == -1) return this;
  if (!has_children())
    allocate_children();
  if (!has_child(i)) {
    child[i] = new IndexedSuffixNode(s, e, path_length + e - s, 1, this);
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
	return child[i];
      }
      else {
	return child[i]->InsertUnder(text, s + child[i]->edge_length(), e,
				     depth - child[i]->edge_length(), index);
      }
    }
  }
  return this;
}


/*************************
 * This is the interface to building a suffix tree.
 ************************************/
void
IndexedSuffixNode::BuildSuffixTree(string &text) {
  
  IndexedSuffixNode *current_node = this->link = this;
  IndexedSuffixNode *previous_node = 0;
  
  // loop to insert each suffix
  for (size_t i = 0; i < text.length(); i++) {
    if (text[i] == 'N' || text[i] == '$') {
      previous_node = 0;
      current_node = this;
    }
    else {
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
      const size_t end = text.length();
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
}


void
IndexedSuffixNode::adjust_links() {
  if (!has_children()) return;
  size_t leftmost = 0, rightmost = 0;
  bool found_leftmost = false;
  for (size_t i = 0; i < alpha_plus_two; ++i)
    if (has_child(i)) {
      child[i]->adjust_links();
      count += child[i]->count;
      if (!found_leftmost) {
	leftmost = i;
	found_leftmost = true;
      }
      rightmost = i;
    }
  parent = (child[leftmost]->has_children()) ?
    child[leftmost]->parent : child[leftmost];
  link = (child[rightmost]->has_children()) ?
    child[rightmost]->link : child[rightmost];
}


void
IndexedSuffixNode::thread_leaves(IndexedSuffixNode** prev_leaf) {
  if (!has_children()) {
    if (*prev_leaf)
      (*prev_leaf)->link = this;
    *prev_leaf = this;
    count = 1;
  }
  else for (size_t i = 0; i < alpha_plus_two; ++i)
    if (child[i])
      child[i]->thread_leaves(prev_leaf);
}


void
IndexedSuffixNode::collect_leaf_indices(vector<size_t>& indices) const {
  for (const IndexedSuffixNode *curr = parent; curr != link; curr = curr->link)
    indices.push_back(curr->end - curr->path_length);
  indices.push_back(link->end - link->path_length);
}


void
IndexedSuffixNode::scores_greater(const char *text, vector<valnode> &matches, 
				  float cutoff, float **scoremat, size_t width, 
				  size_t depth, float score) const {
  if (child)
    for (size_t i = 0; i < alphabet_size; i++) {
      if (child[i]) {
	float temp_score = 0;
	size_t j = child[i]->MatchBranch(text, scoremat, width, depth, score,
					 cutoff, temp_score);
    if (temp_score >= cutoff){ //ab dnb
	  if (depth + j < width)
	    child[i]->scores_greater(text, matches, cutoff, scoremat,
				     width, depth + j, temp_score);
      else {
          matches.push_back(valnode(temp_score, child[i]));
            }
        } //ab dnb
      }
    }
}


void 
IndexedSuffixNode::top_scores(const char *text, hit_queue &top_matches, size_t &pqsize,
			      size_t n_top, float **scoremat, size_t width, 
			      size_t depth, float score, float &cutoff) const {
  if (child)
    for (size_t i = 0; i < alphabet_size; i++)
      if (child[i]) {
	float temp_score = 0;
	size_t j = child[i]->MatchBranch(text, scoremat, width, depth, score, 
					 cutoff, temp_score);
    if (temp_score >= cutoff){ //ab dnb
	  if (depth + j < width)
	    child[i]->top_scores(text, top_matches, pqsize, n_top, scoremat, 
				 width, depth + j, temp_score, cutoff);
	  else {
	    top_matches.push(valnode(temp_score, child[i]));
        pqsize += child[i]->count;
	    while (pqsize - top_matches.top().second->count >= n_top) {
	      pqsize -= top_matches.top().second->count;
	      top_matches.pop();
	    }
	    if (top_matches.size() >= n_top)
	      cutoff = top_matches.top().first;
      }
	  } //ab dnb
      }
}


IndexedSuffixTree::~IndexedSuffixTree() {
  delete root;
}

IndexedSuffixTree::IndexedSuffixTree(const string& s) : sequence(s) {
  sequence.append("$");
  root = new IndexedSuffixNode();
  root->BuildSuffixTree(sequence);
  IndexedSuffixNode **prev_leaf = new IndexedSuffixNode *;
  *prev_leaf = static_cast<IndexedSuffixNode*>(0);
  root->thread_leaves(prev_leaf);
  root->adjust_links();
}


string
IndexedSuffixTree::tostring() const {
  return sequence + string("\n") + root->tostring(sequence.c_str(), 0);
}


float**
IndexedSuffixTree::PrepareScoringMatrix(const ScoringMatrix& matrix, 
					float &max_score) {
  float **scoremat = new float *[matrix.get_width()];
  max_score = 0.0;
  for (size_t i = 0; i < matrix.get_width(); i++) {
    scoremat[i] = new float[alpha_plus_two];
    float max_column_score = 0.0;
    for (size_t j = 0; j < alphabet_size; j++) {
      scoremat[i][j] = matrix[i][j];
      max_column_score = std::max(scoremat[i][j], max_column_score);
    }
    for (size_t j = 0; j < alphabet_size; j++)
      scoremat[i][j] = max_column_score - scoremat[i][j];
    scoremat[i][alphabet_size] = std::numeric_limits<float>::max();
    scoremat[i][alphabet_size + 1] = std::numeric_limits<float>::max();
    max_score += max_column_score;
  }
  return scoremat;
}


void
IndexedSuffixTree::destroy_st_scoring_matrix(float **sm, const size_t w) {
  for (size_t i = 0; i < w; i++)
    delete[] sm[i];
  delete[] sm;
}


void
IndexedSuffixTree::scores_greater_indices(const ScoringMatrix& sm,
					  const float cutoff,
					  vector<pos_score> &final_matches) const {
  float max_score = 0;
  float **scoremat = PrepareScoringMatrix(sm, max_score);
  vector<valnode> matching_nodes;
  root->scores_greater(sequence.c_str(), matching_nodes, 
		       cutoff, scoremat, sm.get_width(), 0, max_score);
  destroy_st_scoring_matrix(scoremat, sm.get_width());
  // TODO: make sure matching_nodes is sorted properly
  vector<valnode>::iterator i;
  for (i = matching_nodes.begin(); i != matching_nodes.end(); ++i) {
    vector<size_t> indices;
    if (i->second->is_leaf())
      final_matches.push_back(pos_score(i->second->end - 
					i->second->path_length, i->first));
    else {
      i->second->collect_leaf_indices(indices);
      for (size_t j = 0; j < indices.size(); ++j)
	final_matches.push_back(pos_score(indices[j], i->first));
    }
  }
}


void 
IndexedSuffixTree::top_scores_indices(const ScoringMatrix& matrix,
				      const size_t n_top,
				      vector<pos_score> &final_matches) const {
  
  // TODO: define an exception for this
  if (n_top == 0) return;
  
  float max_score = 0;
  float **scoremat = PrepareScoringMatrix(matrix, max_score);
  float cutoff = -std::numeric_limits<float>::max();
  hit_queue top_matching_nodes;
  size_t pqsize = 0;
  root->top_scores(sequence.c_str(), top_matching_nodes, pqsize, n_top, 
		   scoremat, matrix.get_width(), 0, max_score, cutoff);
  destroy_st_scoring_matrix(scoremat, matrix.get_width());
  while (top_matching_nodes.size() > 0) {
    vector<size_t> indices;
    valnode vn(top_matching_nodes.top());
    if (vn.second->is_leaf())
      final_matches.push_back(pos_score(vn.second->end - 
					vn.second->path_length, vn.first));
    else {
      vn.second->collect_leaf_indices(indices);
      for (size_t j = 0; j < indices.size(); ++j)
	final_matches.push_back(pos_score(indices[j], vn.first));
    }
    top_matching_nodes.pop();
  }
  std::reverse(final_matches.begin(), final_matches.end());
  final_matches.erase(final_matches.begin() +
		      std::min(n_top, final_matches.size()),
		      final_matches.end());
}
