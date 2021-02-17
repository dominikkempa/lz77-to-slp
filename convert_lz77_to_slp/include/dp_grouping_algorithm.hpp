#ifndef __DP_GROUPING_ALGORITHM_HPP_INCLUDED
#define __DP_GROUPING_ALGORITHM_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <vector>

#include "avl_grammar_node.hpp"
#include "hash_table.hpp"
#include "karp_rabin_hashing.hpp"


template<typename char_type>
std::vector<const avl_grammar_node<char_type>*> dp_grouping_algorithm(
    const hash_table<std::uint64_t, const avl_grammar_node<char_type>*> &hashes,
    const std::vector<const avl_grammar_node<char_type>*> &seq,
    const std::uint64_t hash_variable,
    const std::uint64_t mersenne_prime_exponent) {

#if 1
  // Create the vector to hold the solution.
  std::vector<const avl_grammar_node<char_type>*> ret;

  // Handle special case.
  if (seq.empty())
    return ret;

  // Declate typedefs
  typedef const avl_grammar_node<char_type>* ptr_type;

  // Allocate the DP array.
  std::uint64_t length = seq.size();
  std::uint64_t **dp = new std::uint64_t*[length];
  std::uint64_t **dp_sol = new std::uint64_t*[length];
  ptr_type **dp_nonterm = new ptr_type*[length];
  for (std::uint64_t i = 0; i < length; ++i) {
    dp[i] = new std::uint64_t[length];
    dp_sol[i] = new std::uint64_t[length];
    dp_nonterm[i] = new ptr_type[length];
  }

  // Fill in the DP array.
  // Start with len = 1.
  for (std::uint64_t i = 0; i < length; ++i) {
    dp[i][i] = 1;
    dp_sol[i][i] = 1;
    dp_nonterm[i][i] = seq[i];
  }

  // Solve for subarray of length > 1.
  for (std::uint64_t len = 2; len <= length; ++len) {
    for (std::uint64_t beg = 0; beg <= length - len; ++beg) {
      const std::uint64_t end = beg + len - 1;

      // Initialize to solution to the initial choice.
      dp[beg][end] = 1 + dp[beg + 1][end];
      dp_sol[beg][end] = 1;
      dp_nonterm[beg][end] = seq[beg];

      // Try all other possible choices.
      std::uint64_t h = seq[beg]->m_kr_hash;
      for (std::uint64_t leftlen = 2; leftlen <= len; ++leftlen) {
        const std::uint64_t last = beg + leftlen - 1;
        
        // Update rolling hash.
        h = mod_mersenne(
            mul_mod_meresenne(
              h,
              pow_mod_mersenne(
                hash_variable,
                seq[last]->m_exp_len,
                mersenne_prime_exponent),
              mersenne_prime_exponent) + seq[last]->m_kr_hash,
            mersenne_prime_exponent);

        // If a nonterminal expanding to the substring with the
        // same hash already exists in the grammar, check if
        // creating the block results in a better solution than
        // the currently optimal solution.
        typedef const avl_grammar_node<char_type>* value_type;
        const value_type *nonterm = hashes.find(h);
        if (nonterm != NULL) {
          std::uint64_t sol_cost = 1;
          if (leftlen < len)
            sol_cost += dp[last + 1][end];
          if (sol_cost < dp[beg][end]) {
            dp[beg][end] = sol_cost;
            dp_sol[beg][end] = leftlen;
            dp_nonterm[beg][end] = *nonterm;
          }
        }
      }
    }
  }

  // Restore the optimal solution.
  std::uint64_t prefix_length = 0;
  while (prefix_length < length) {
    ret.push_back(dp_nonterm[prefix_length][length - 1]);
    prefix_length += dp_sol[prefix_length][length - 1];
  }

  // Clean up.
  for (std::uint64_t i = 0; i < length; ++i) {
    delete[] dp[i];
    delete[] dp_sol[i];
    delete[] dp_nonterm[i];
  }
  delete[] dp;
  delete[] dp_sol;
  delete[] dp_nonterm;

  // Return the result.
  return ret;
#else
  return seq;
#endif
}

#endif  // __DP_GROUPING_ALGORITHM_HPP_INCLUDED

