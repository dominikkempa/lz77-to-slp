#ifndef __CONVERT_LZ77_TO_AVL_GRAMMAR_HPP_INCLUDED
#define __CONVERT_LZ77_TO_AVL_GRAMMAR_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "hash_table.hpp"
#include "karp_rabin_hashing.hpp"
#include "avl_grammar_node.hpp"
#include "avl_grammar_add_concat_nonterminal.hpp"
#include "avl_grammar_add_substring_nonterminal.hpp"
#include "avl_grammar.hpp"


//=============================================================================
// Given the LZ77 parsing (consisting of z phrases) of text T, compute
// the AVL grammar of size O(z log n) expanding to T.  The parsing is
// given as a sequence of pairs (pos, len), where either len > 0 and
// pos encodes the position of the previous occurrence in the string,
// or len = 0 and then pos contain the text symbol.
// TODO: the grammar at this point is guaranteed to be of size O(z log
//       n), but there might be some unused nonterminals. They should
//       be removed. Anyway, at this point, I just want to test the
//       correctness of the conversion.
//=============================================================================
template<
  typename char_type,
  typename text_offset_type>
avl_grammar<char_type> *convert_lz77_to_avl_grammar(
    const std::vector<
      std::pair<text_offset_type, text_offset_type> > &parsing) {

  // Declare types.
  typedef avl_grammar_node<char_type> node_type;
  typedef avl_grammar<char_type> grammar_type;

  // Set hashing variables.
  const std::uint64_t mersenne_prime_exponent = 61;
  const std::uint64_t hash_variable =
    rand_mod_mersenne(mersenne_prime_exponent);

  // Compute the AVL grammar expanding to T.
  grammar_type *grammar = new grammar_type(
      hash_variable, mersenne_prime_exponent);
  std::uint64_t prefix_length = 0;
  for (std::uint64_t phrase_id = 0; phrase_id < parsing.size();
      ++phrase_id) {
    std::pair<text_offset_type, text_offset_type> p = parsing[phrase_id];
    std::uint64_t pos = p.first;
    std::uint64_t len = p.second;
    
    // Compute the AVL grammar expanding to phrase p.
    const node_type *phrase_root = NULL;
    if (len == 0) {

      // If this is a literal phrase, create a trivial grammar.
      phrase_root = new node_type((char_type)pos, mersenne_prime_exponent);
      grammar->m_nonterminals.push_back(phrase_root);
    } else {

      // We proceed differently, depending on whether
      // the phrase is self-overlapping. This part is
      // unadressed in the original Rytter's paper. The
      // solution is described in the proof of Theorem 6.1
      // in https://arxiv.org/abs/1910.10631v3.
      if (pos + len > prefix_length) {

        // If the phase is self-overlapping, we create the
        // nonterminal expanding to text[pos..prefix_length).
        const node_type * const suffix_nonterminal =
          add_substring_nonterminal<char_type>(
              grammar->m_nonterminals, grammar->m_root, pos, prefix_length,
              hash_variable, mersenne_prime_exponent);

        // Square the above nonterminal until
        // it reaches length >= len.
        const node_type *suffix_pow_nonterminal = suffix_nonterminal;
        std::uint64_t curlen = prefix_length - pos;
        while (curlen < len) {
          const node_type * const square =
            new node_type(suffix_pow_nonterminal, suffix_pow_nonterminal,
                hash_variable, mersenne_prime_exponent);
          grammar->m_nonterminals.push_back(square);
          curlen <<= 1;
          suffix_pow_nonterminal = square;
        }

        // Create a nonterminal expanding to the prefix
        // of exp(suffix_pow_nonterminal) of length len.
        phrase_root = add_substring_nonterminal<char_type>(
            grammar->m_nonterminals, suffix_pow_nonterminal, 0, len,
            hash_variable, mersenne_prime_exponent);
      } else {

        // Add the nonterminal expanding to phrase p.
        std::uint64_t begin = pos;
        std::uint64_t end = begin + len;
        phrase_root = add_substring_nonterminal<char_type>(
            grammar->m_nonterminals, grammar->m_root, begin, end,
            hash_variable, mersenne_prime_exponent);
      }
    }

    // Update prefix_grammar to encode the longer prefix.
    if (grammar->m_root == NULL)
      grammar->m_root = phrase_root;
    else
      grammar->m_root =
        add_concat_nonterminal<char_type>(
            grammar->m_nonterminals, grammar->m_root, phrase_root,
            hash_variable, mersenne_prime_exponent);

    // Update prefix_length.
    prefix_length += std::max(len, (std::uint64_t)1);
  }

  // Return the result.
  return grammar;
}

#endif  // __CONVERT_LZ77_TO_AVL_GRAMMAR_HPP_INCLUDED
