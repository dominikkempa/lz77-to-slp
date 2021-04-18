#ifndef __CONVERT_LZ77_TO_AVL_GRAMMAR_HPP_INCLUDED
#define __CONVERT_LZ77_TO_AVL_GRAMMAR_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "../utils/karp_rabin_hashing.hpp"
#include "avl_grammar_node.hpp"
#include "avl_grammar_add_concat_nonterminal.hpp"
#include "avl_grammar.hpp"


//=============================================================================
// Given the LZ77-like parsing of size z of text T, compute the AVL
// grammar of size O(z log n) expanding to T. The parsing is given
// as a sequence of pairs (pos, len), where either len > 0 and
// pos encodes the position of the previous occurrence in the string,
// or len = 0 and then pos contain the text symbol.
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

  // Init Karp-Rabin hashing.
  karp_rabin_hashing::init();

  // Compute the AVL grammar expanding to T.
  grammar_type *grammar = new grammar_type();
  std::uint64_t prefix_length = 0;
  for (std::uint64_t phrase_id = 0;
      phrase_id < parsing.size(); ++phrase_id) {

    // Get the next phrase len and pos values.
    std::pair<text_offset_type, text_offset_type> p =
      parsing[phrase_id];
    std::uint64_t pos = p.first;
    std::uint64_t len = p.second;
    
    // Compute the AVL grammar expanding to phrase p.
    const node_type *phrase_root = NULL;
    if (len == 0) {

      // If this is a literal phrase, create a trivial grammar.
      phrase_root = new node_type((char_type)pos);
      grammar->add_nonterminal(phrase_root);
    } else {

      // Self-overlapping phrases are unadressed in Rytter's paper.
      // We follow Theorem 6.1 in https://arxiv.org/abs/1910.10631v3.
      if (pos + len > prefix_length) {

        // If the phase is self-overlapping, we create the
        // nonterminal expanding to text[pos..prefix_length).
        const node_type * const suffix_nonterm =
          grammar->add_substring_nonterminal(pos, prefix_length);

        // Square the nonterminal until it reaches length >= len.
        const node_type *suffix_pow_nonterm = suffix_nonterm;
        std::uint64_t curlen = prefix_length - pos;
        while (curlen < len) {
          const node_type * const square = new node_type(
              suffix_pow_nonterm, suffix_pow_nonterm);
          grammar->add_nonterminal(square);
          curlen <<= 1;
          suffix_pow_nonterminal = square;
        }

        // Create a nonterminal expanding to the prefix length len.
        phrase_root = grammar->add_substring_nonterminal(
            suffix_pow_nonterm, 0, len);
      } else {

        // Add the nonterminal expanding to phrase p.
        std::uint64_t begin = pos;
        std::uint64_t end = begin + len;
        phrase_root = grammar->add_substring_nonterminal(begin, end);
      }
    }

    // Update prefix length and add new root to the grammar.
    prefix_length += std::max(len, (std::uint64_t)1);
    if (grammar->m_root == NULL) grammar->m_root = phrase_root;
    else grammar->m_root =
      add_concat_nonterminal<char_type>(
          grammar->m_hashes, grammar->m_nonterminals,
          grammar->m_root, phrase_root);
  }

  // Return the result.
  return grammar;
}

#endif  // __CONVERT_LZ77_TO_AVL_GRAMMAR_HPP_INCLUDED
