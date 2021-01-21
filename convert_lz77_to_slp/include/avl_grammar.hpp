#ifndef __AVL_GRAMMAR_HPP_INCLUDED
#define __AVL_GRAMMAR_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>


//=============================================================================
// A class used to represent the nonterminal in the AVL grammar.
//=============================================================================
template<typename char_type>
struct avl_grammar_node {
  typedef avl_grammar_node<char_type> node_type;

  // Class members.
  char_type m_char;
  std::uint8_t m_height;
  std::uint64_t m_exp_len;
  const node_type * const m_left;
  const node_type * const m_right;

  // Default constructor.
  avl_grammar_node() :
      m_left(NULL),
      m_right(NULL) {
    m_char = (char_type)0;
    m_height = 0;
    m_exp_len = 1;
  }

  // Constructor for a node expanding to a single symbol.
  avl_grammar_node(const char_type c) :
      m_left(NULL),
      m_right(NULL) {
    m_char = c;
    m_height = 0;
    m_exp_len = 1;
  }

  // Constructor for a standard nonterminal
  // (expanding to two other nonterminals).
  avl_grammar_node(
      const node_type * const left,
      const node_type * const right) :
        m_left(left),
        m_right(right) {
    m_char = (char_type)0;
    m_height = std::max(left->m_height, right->m_height) + 1;
    m_exp_len = left->m_exp_len + right->m_exp_len;
  }

  // Print the string encoded by the grammar.
  void print_expansion() const {
    if (m_height == 0)
      fprintf(stderr, "%c", (char)m_char);
    else {
      m_left->print_expansion();
      m_right->print_expansion();
    }
  }
};

//=============================================================================
// Given the AVL grammars expanding to strings X and Y, compute the
// AVL grammar expanding to XY.
//=============================================================================
template<typename char_type>
avl_grammar_node<char_type> *merge_avl_grammars(
    std::vector<avl_grammar_node<char_type> *> &nonterminals,
    const avl_grammar_node<char_type> * const left,
    const avl_grammar_node<char_type> * const right) {
  
  // Declare type.
  typedef avl_grammar_node<char_type> node_type;

  // Consider two cases, depending on whether
  // left of right nonterminal is taller.
  if (left->m_height >= right->m_height) {
    if (left->m_height - right->m_height <= 1) {

      // Height are close. Just merge and return.
      node_type *newroot = new node_type(left, right);
      nonterminals.push_back(newroot);
      return newroot;
    } else {
      node_type *newright =
        merge_avl_grammars<char_type>(nonterminals, left->m_right, right);
      if (newright->m_height > left->m_left->m_height &&
          newright->m_height - left->m_left->m_height > 1) {
        
        // Rebalancing needed.
        if (newright->m_left->m_height > newright->m_right->m_height) {

          // Double (right-left) rotation.
          node_type *X =
            new node_type(left->m_left, newright->m_left->m_left);
          node_type *Z =
            new node_type(newright->m_left->m_right, newright->m_right);
          node_type *Y = new node_type(X, Z);
          nonterminals.push_back(X);
          nonterminals.push_back(Y);
          nonterminals.push_back(Z);
          return Y;
        } else {

          // Single (left) rotation.
          node_type *X = new node_type(left->m_left, newright->m_left);
          node_type *Y = new node_type(X, newright->m_right);
          nonterminals.push_back(X);
          nonterminals.push_back(Y);
          return Y;
        }
      } else {

        // No need to rebalance.
        node_type *newroot = new node_type(left->m_left, newright);
        nonterminals.push_back(newroot);
        return newroot;
      }
    }
  } else {
    if (right->m_height - left->m_height <= 1) {

      // Heights are close. Just merge and return.
      node_type *newroot = new node_type(left, right);
      nonterminals.push_back(newroot);
      return newroot;
    } else {
      node_type *newleft =
        merge_avl_grammars<char_type>(nonterminals, left, right->m_left);
      if (newleft->m_height > right->m_right->m_height &&
          newleft->m_height - right->m_right->m_height > 1) {

        // Rebalancing needed.
        if (newleft->m_right->m_height > newleft->m_left->m_height) {

          // Double (left-right) rotation.
          node_type *X =
            new node_type(newleft->m_left, newleft->m_right->m_left);
          node_type *Z =
            new node_type(newleft->m_right->m_right, right->m_right);
          node_type *Y = new node_type(X, Z);
          nonterminals.push_back(X);
          nonterminals.push_back(Y);
          nonterminals.push_back(Z);
          return Y;
        } else {

          // Single (right) rotation.
          node_type *Y = new node_type(newleft->m_right, right->m_right);
          node_type *X = new node_type(newleft->m_left, Y);
          nonterminals.push_back(X);
          nonterminals.push_back(Y);
          return X;
        }
      } else {

        // No need to rebalance.
        node_type *newroot = new node_type(newleft, right->m_right);
        nonterminals.push_back(newroot);
        return newroot;
      }
    }
  }
}

//=============================================================================
// Given the AVL grammar expanding to string T, compute the AVL
// grammar expanding to T[begin..end).
//=============================================================================
template<typename char_type>
avl_grammar_node<char_type> *create_substring_avl_grammar(
    std::vector<avl_grammar_node<char_type> *> &nonterminals,
    const avl_grammar_node<char_type> * const root,
    const std::uint64_t begin,
    const std::uint64_t end) {

  // Check input correctness.
  if (begin > end ||
      end > root->m_exp_len) {
    fprintf(stderr, "\nError: extract: end > root->m_exp_len!\n");
    std::exit(EXIT_FAILURE);
  }
  
  // Declare types.
  typedef avl_grammar_node<char_type> node_type;

  // Handle boundary case.
  if (begin == end)
    return (node_type *)NULL;

  // Find the deepest node in the parse tree containing the range
  // [begin..end).
  node_type *x = root;
  std::uint64_t cur_range_beg = 0;
  std::uint64_t cur_range_end = x->m_exp_len;
  while (cur_range_end - cur_range_beg > end - begin &&
      (end <= cur_range_beg + x->m_left->m_exp_len ||
       begin >= cur_range_beg + x->m_left->m_exp_len)) {
    if (end <= cur_range_beg + x->m_left->m_exp_len) {
      cur_range_end = cur_range_beg + x->m_left->m_exp_len;
      x = x->m_left;
    } else {
      cur_range_beg += x->m_left->m_exp_len;
      x = x->m_right;
    }
  }

  // Check if the range of x is exactly [begin..end).
  if (cur_range_beg == begin && cur_range_end == end) {

    // If yes, return x as the answer.
    return x;
  } else {

    // Otherwise, we perform two traversals in the tree. If by G we
    // denote the AVL grammar x->m_left, we first create the AVL
    // grammar expanding to suffix of length (cur_range_beg +
    // x->m_left->m_exp_len - begin) of the string exp(G). The
    // root of the resulting grammar is stored in left_grammar.
    node_type *left_grammar = NULL;
    {
      std::uint64_t left_range_beg = cur_range_beg;
      std::uint64_t left_range_end = cur_range_beg + x->m_left->m_exp_len;
      std::uint64_t suffix_length = left_range_end - begin;
      node_type *y = x->m_left;

      // The tricky thing here is that we cannot merge the nodes right
      // away. For the merge to be O(log n) time, we need to merge
      // grammar short-to-tall and hence we first collect their roots
      // in a vector.
      std::vector<node_type*> grammars_to_merge;
      while (suffix_length > 0) {
        if (y->m_exp_len == suffix_length) {
          grammars_to_merge.push_back(y);
          suffix_length -= y->m_exp_len;
        } else if (suffix_length > y->m_right->m_exp_len) {
          grammars_to_merge.push_back(y->m_right);
          suffix_length -= y->m_right->m_exp_len;
          y = y->m_left;
        } else y = y->m_right;
      }

      // Merge AVL grammars in grammars_to_merge into a single AVL grammar.
      while (grammars_to_merge.size() > 1) {
        node_type *left = grammars_to_merge.back();
        grammars_to_merge.pop_back();
        node_type *right = grammars_to_merge.back();
        grammars_to_merge.pop_back();
        grammars_to_merge.push_back(
            merge_avl_grammars<char_type>(nonterminals, left, right));
      }
      left_grammar = grammars_to_merge.back();
    }

    // Perform the analogous operation for the right side.
    node_type *right_grammar = NULL;
    {
      std::uint64_t right_range_beg = cur_range_beg + x->m_left->m_exp_len;
      std::uint64_t right_range_end = cur_range_end;
      std::uint64_t prefix_length = end - right_range_beg;
      node_type *y = x->m_right;
      
      // Collect the roots of grammars to merge.
      std::vector<node_type*> grammars_to_merge;
      while (prefix_length > 0) {
        if (y->m_exp_len == prefix_length) {
          grammars_to_merge.push_back(y);
          prefix_length -= y->m_exp_len;
        } else if (prefix_length > y->m_left->m_exp_len) {
          grammars_to_merge.push_back(y->m_left);
          prefix_length -= y->m_left->m_exp_len;
          y = y->m_right;
        } else y = y->m_left;
      }

      // Merge AVL grammars in grammars_to_merge into a single AVL grammar.
      while (grammars_to_merge.size() > 1) {
        node_type *right = grammars_to_merge.back();
        grammars_to_merge.pop_back();
        node_type *left = grammars_to_merge.back();
        grammars_to_merge.pop_back();
        grammars_to_merge.push_back(
            merge_avl_grammars<char_type>(nonterminals, left, right));
      }
      right_grammar = grammars_to_merge.back();
    }

    // Merge left_grammar with right_grammar and return.
    // Both are guaranteed to be non-NULL.
    node_type *final_grammar =
      merge_avl_grammars<char_type>(nonterminals,
          left_grammar, right_grammar);

    // Return the result.
    return final_grammar;
  }
}

#endif  // __AVL_GRAMMAR_HPP_INCLUDED
