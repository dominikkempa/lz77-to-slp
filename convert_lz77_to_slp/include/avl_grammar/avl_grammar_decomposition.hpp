#ifndef __AVL_GRAMMAR_DECOMPOSITION_HPP_INCLUDED
#define __AVL_GRAMMAR_DECOMPOSITION_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "avl_grammar_node.hpp"


//=============================================================================
// Given the nonterminal of the AVL grammar expanding to string S,
// return the sequence of nonterminals expanding to S[begin..end).
//=============================================================================
template<typename char_type>
std::vector<const avl_grammar_node<char_type>*> decomposition(
    const avl_grammar_node<char_type> * const root,
    const std::uint64_t begin,
    const std::uint64_t end) {

  // Check input correctness.
  if (begin > end ||
      end > root->m_exp_len) {
    fprintf(stderr, "\nError: decomposition: invalid range!\n");
    std::exit(EXIT_FAILURE);
  }
  
  // Declare types.
  typedef avl_grammar_node<char_type> node_type;

  // Declare the vector storing the result.
  std::vector<const node_type*> decomposition;

  // Handle boundary case.
  if (begin == end)
    return decomposition;

  // Find the deepest node in the parse tree containing the range
  // [begin..end).
  const node_type *x = root;
  std::uint64_t cur_range_beg = 0;
  std::uint64_t cur_range_end = x->m_exp_len;
  while (x->m_height > 0 &&
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
    decomposition.push_back(x);
  } else {

    // Otherwise, we perform two traversals in the tree. If by G we
    // denote the AVL grammar x->m_left, we first create the AVL
    // grammar expanding to suffix of length (cur_range_beg +
    // x->m_left->m_exp_len - begin) of the string exp(G). The
    // root of the resulting grammar is stored in left_grammar.
    {
      const std::uint64_t left_range_end =
        cur_range_beg + x->m_left->m_exp_len;
      std::uint64_t suffix_length = left_range_end - begin;
      const node_type *y = x->m_left;

      // The tricky thing here is that we cannot merge the nodes right
      // away. For the merge to be O(log n) time, we need to merge
      // grammar short-to-tall and hence we first collect their roots
      // in a vector.
      while (suffix_length > 0) {
        if (y->m_exp_len == suffix_length) {
          decomposition.push_back(y);
          suffix_length -= y->m_exp_len;
        } else if (suffix_length > y->m_right->m_exp_len) {
          decomposition.push_back(y->m_right);
          suffix_length -= y->m_right->m_exp_len;
          y = y->m_left;
        } else y = y->m_right;
      }
    }

    // Reverse the first sequence of nonterminals
    // collected during the left downward traversal.
    std::reverse(decomposition.begin(), decomposition.end());

    // Perform the analogous operation for the right side.
    {
      const std::uint64_t right_range_beg =
        cur_range_beg + x->m_left->m_exp_len;
      std::uint64_t prefix_length = end - right_range_beg;
      const node_type *y = x->m_right;
      
      // Collect the roots of grammars to merge.
      while (prefix_length > 0) {
        if (y->m_exp_len == prefix_length) {
          decomposition.push_back(y);
          prefix_length -= y->m_exp_len;
        } else if (prefix_length > y->m_left->m_exp_len) {
          decomposition.push_back(y->m_left);
          prefix_length -= y->m_left->m_exp_len;
          y = y->m_right;
        } else y = y->m_left;
      }
    }
  }

  // Return the result.
  return decomposition;
}

#endif  // __AVL_GRAMMAR_DECOMPOSITION_HPP_INCLUDED
