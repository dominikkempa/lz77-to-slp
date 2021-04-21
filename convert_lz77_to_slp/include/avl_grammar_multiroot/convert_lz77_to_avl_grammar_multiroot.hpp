#ifndef __CONVERT_LZ77_TO_AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
#define __CONVERT_LZ77_TO_AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "../utils/karp_rabin_hashing.hpp"
#include "avl_grammar_multiroot.hpp"


//=============================================================================
// Given the LZ77-like parsing of size z of text T, compute the AVL
// multi-root grammar of size O(z log n) expanding to T. The parsing
// is given as a sequence of pairs (pos, len), where either len > 0 and
// pos encodes the position of the previous occurrence in the string,
// or len = 0 and then pos contain the text symbol.
//=============================================================================
template<
  typename char_type,
  typename text_offset_type>
avl_grammar_multiroot<char_type, text_offset_type>*
convert_lz77_to_avl_grammar_multiroot(
    const std::vector<
      std::pair<text_offset_type, text_offset_type> > &parsing) {

  // Declare types.
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  typedef avl_grammar_multiroot<char_type, text_offset_type> grammar_type;

  // Init Karp-Rabin hashing.
  karp_rabin_hashing::init();

  // Compute the AVL grammar expanding to T.
  grammar_type *grammar = new grammar_type();
  std::uint64_t prefix_length = 0;
  for (std::uint64_t phrase_id = 0;
      phrase_id < parsing.size(); ++phrase_id) {

    // Check if we need to run garbage collector.
    grammar->check_gargage_collector();

    // Get the next phrase len and pos values.
    std::pair<text_offset_type, text_offset_type> p =
      parsing[phrase_id];
    std::uint64_t pos = p.first;
    std::uint64_t len = p.second;
    
    // Compute the AVL grammar expanding to phrase p.
    std::vector<text_offset_type> phrase_roots;
    if (len == 0) {

      // If this is a literal phrase, create a trivial grammar.
      const nonterminal_type root((char_type)pos);
      const std::uint64_t root_id = grammar->add_nonterminal(root);
      phrase_roots.push_back((text_offset_type)root_id);
    } else {

      // Self-overlapping phrases are unadressed in Rytter's paper.
      // We follow Theorem 6.1 in https://arxiv.org/abs/1910.10631v3.
      if (pos + len > prefix_length) {
        std::uint64_t left = len;
        std::uint64_t exist = prefix_length - pos;
        while (left > 0) {

          std::uint64_t next = std::min(left, exist);
          grammar->merge_enclosed_roots(pos, pos + next);
          std::vector<text_offset_type> v =
            grammar->decomposition(pos, pos + next);
          v = grammar->find_equivalent_seq(v);

          for (std::uint64_t t = 0; t < v.size(); ++t) {
            const std::uint64_t id = v[t];
            const std::uint64_t exp_len = grammar->get_exp_len(id);
            prefix_length += exp_len;
            grammar->push_root(prefix_length, v[t]);
          }
          exist += next;
          left -= next;
        }
        continue;
      } else {

        // Add the nonterminal expanding to phrase p.
        std::uint64_t begin = pos;
        std::uint64_t end = begin + len;
        grammar->merge_enclosed_roots(begin, end);
        phrase_roots = grammar->decomposition(begin, end);
        phrase_roots = grammar->find_equivalent_seq(phrase_roots);
      }
    }

    // Update prefix length and add new roots to the grammar.
    for (std::uint64_t t = 0; t < phrase_roots.size(); ++t) {
      const std::uint64_t id = phrase_roots[t];
      const std::uint64_t exp_len = grammar->get_exp_len(id);
      prefix_length += exp_len;
      grammar->push_root(prefix_length, phrase_roots[t]);
    }
  }

  // Return the result.
  return grammar;
}

#endif  // __CONVERT_LZ77_TO_AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
