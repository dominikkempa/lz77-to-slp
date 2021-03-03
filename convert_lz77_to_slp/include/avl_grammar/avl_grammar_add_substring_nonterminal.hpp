#ifndef __AVL_GRAMMAR_ADD_SUBSTRING_NONTERMINAL_HPP_INCLUDED
#define __AVL_GRAMMAR_ADD_SUBSTRING_NONTERMINAL_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "avl_grammar_node.hpp"
#include "avl_grammar_add_concat_nonterminal.hpp"
#include "../utils/hash_table.hpp"


//=============================================================================
// Given the nonterminal of the AVL grammar expanding to string S,
// create and return the nonterminal expanding to S[begin..end).
// The root of the input grammar remains unchanged.
//=============================================================================
template<typename char_type>
const avl_grammar_node<char_type> *add_substring_nonterminal(
    hash_table<std::uint64_t, const avl_grammar_node<char_type>*> &hashes,
    std::vector<const avl_grammar_node<char_type> *> &nonterminals,
    const avl_grammar_node<char_type> * const root,
    const std::uint64_t begin,
    const std::uint64_t end,
    const std::uint64_t hash_variable,
    const std::uint64_t mersenne_prime_exponent) {

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
    return x;
  } else {

    // Otherwise, we perform two traversals in the tree. If by G we
    // denote the AVL grammar x->m_left, we first create the AVL
    // grammar expanding to suffix of length (cur_range_beg +
    // x->m_left->m_exp_len - begin) of the string exp(G). The
    // root of the resulting grammar is stored in left_grammar.
    const node_type * left_grammar = NULL;
    {
      const std::uint64_t left_range_end =
        cur_range_beg + x->m_left->m_exp_len;
      std::uint64_t suffix_length = left_range_end - begin;
      const node_type *y = x->m_left;

      // The tricky thing here is that we cannot merge the nodes right
      // away. For the merge to be O(log n) time, we need to merge
      // grammar short-to-tall and hence we first collect their roots
      // in a vector.
      std::vector<const node_type*> grammars_to_merge;
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
        const node_type * const left = grammars_to_merge.back();
        grammars_to_merge.pop_back();
        const node_type * const right = grammars_to_merge.back();
        grammars_to_merge.pop_back();
        grammars_to_merge.push_back(
            add_concat_nonterminal<char_type>(
              hashes, nonterminals, left, right,
              hash_variable, mersenne_prime_exponent));
      }
      left_grammar = grammars_to_merge.back();
    }

    // Perform the analogous operation for the right side.
    const node_type *right_grammar = NULL;
    {
      const std::uint64_t right_range_beg =
        cur_range_beg + x->m_left->m_exp_len;
      std::uint64_t prefix_length = end - right_range_beg;
      const node_type *y = x->m_right;
      
      // Collect the roots of grammars to merge.
      std::vector<const node_type*> grammars_to_merge;
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
        const node_type * const right = grammars_to_merge.back();
        grammars_to_merge.pop_back();
        const node_type * const left = grammars_to_merge.back();
        grammars_to_merge.pop_back();
        grammars_to_merge.push_back(
            add_concat_nonterminal<char_type>(
              hashes, nonterminals, left, right,
              hash_variable, mersenne_prime_exponent));
      }
      right_grammar = grammars_to_merge.back();
    }

    // Merge left_grammar with right_grammar and return.
    // Both are guaranteed to be non-NULL.
    const node_type * const final_grammar =
      add_concat_nonterminal<char_type>(
          hashes, nonterminals, left_grammar, right_grammar,
          hash_variable, mersenne_prime_exponent);

    // Return the result.
    return final_grammar;
  }
}

#endif  // __AVL_GRAMMAR_ADD_SUBSTRING_NONTERMINAL_HPP_INCLUDED
