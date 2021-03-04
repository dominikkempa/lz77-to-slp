#ifndef __CONVERT_LZ77_TO_AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
#define __CONVERT_LZ77_TO_AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "../utils/karp_rabin_hashing.hpp"
#include "avl_grammar_node.hpp"
#include "avl_grammar_multiroot.hpp"


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
avl_grammar_multiroot<char_type> *convert_lz77_to_avl_grammar_multiroot(
    const std::vector<
      std::pair<text_offset_type, text_offset_type> > &parsing) {

  // Declare types.
  typedef avl_grammar_node<char_type> node_type;
  typedef avl_grammar_multiroot<char_type> grammar_type;

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
    std::vector<const node_type*> phrase_roots;
    if (len == 0) {

      // If this is a literal phrase, create a trivial grammar.
      phrase_root = new node_type((char_type)pos);
      grammar->m_nonterminals.push_back(phrase_root);
      phrase_roots.push_back(phrase_root);
    } else {

      // We proceed differently, depending on whether
      // the phrase is self-overlapping. This part is
      // unadressed in the original Rytter's paper. The
      // solution is described in the proof of Theorem 6.1
      // in https://arxiv.org/abs/1910.10631v3.
      if (pos + len > prefix_length) {
        fprintf(stderr, "\nSelf-overlapping phrase detected!\n");
        std::exit(EXIT_FAILURE);

        // TODO: implement handling of
        // self-overlapping phrases.
      } else {

        // Add the nonterminal expanding to phrase p.
        std::uint64_t begin = pos;
        std::uint64_t end = begin + len;
        phrase_roots =
          grammar->add_substring_nonterminal(begin, end);
      }
    }

    // Update prefix length and add new root to the grammar.
    //prefix_length += phrase_len;
    //grammar->m_roots[prefix_length] = phrase_root;
    for (std::uint64_t t = 0; t < phrase_roots.size(); ++t) {
      prefix_length += phrase_roots[t]->m_exp_len;
      grammar->m_roots[prefix_length] = phrase_roots[t];
    }
  }

  // Return the result.
  return grammar;
}

#endif  // __CONVERT_LZ77_TO_AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
