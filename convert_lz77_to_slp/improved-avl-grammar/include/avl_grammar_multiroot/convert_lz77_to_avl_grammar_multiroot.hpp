#ifndef __CONVERT_LZ77_TO_AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
#define __CONVERT_LZ77_TO_AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>

#include "../utils/karp_rabin_hashing.hpp"
#include "../utils/packed_pair.hpp"
#include "../utils/space_efficient_vector.hpp"
#include "../io/async_stream_reader.hpp"
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
    const std::string parsing_filename,
    bool use_kr_hashing,
    long double kr_hashing_prob) {

  // Start the timer.
  long double start = utils::wclock();

  // Declare types.
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  typedef text_offset_type ptr_type;
  typedef avl_grammar_multiroot<char_type, text_offset_type> grammar_type;
  typedef async_stream_reader<text_offset_type> reader_type;
  typedef packed_pair<ptr_type, text_offset_type> pair_type;

  // Initialize the parsing reader.
  const std::uint64_t bufsize = (1 << 19);
  const std::uint64_t n_buffers = 4;
  reader_type *parsing_reader =
    new reader_type(parsing_filename, bufsize, n_buffers);

  // Compute parsing size.
  const std::uint64_t parsing_size =
    utils::file_size(parsing_filename) / (2 * sizeof(text_offset_type));

  // Init Karp-Rabin hashing.
  karp_rabin_hashing::init();

  // Compute the AVL grammar expanding to T.
  grammar_type *grammar = new grammar_type(use_kr_hashing, kr_hashing_prob);
  std::uint64_t prefix_length = 0;
  for (std::uint64_t phrase_id = 0; phrase_id < parsing_size; ++phrase_id) {

    if (((phrase_id + 1) & ((1 << 16) - 1)) == 0) {
      long double elapsed = utils::wclock() - start;
      fprintf(stderr, "\rInfo: elapsed = %.2Lfs, "
          "progress = %lu (%.2Lf%%) phrases (%.2LfMiB prefix), "
          "peak RAM = %.2LfMiB",
          elapsed, phrase_id + 1,
          (100.L * (phrase_id + 1)) / parsing_size,
          (1.L * prefix_length) / (1 << 20),
          (1.L * utils::get_peak_ram_allocation()) / (1UL << 20));
    }

    // Check if we need to run garbage collector.
    grammar->check_gargage_collector();

    // Get the next phrase len and pos values.
    std::uint64_t pos = parsing_reader->read();
    std::uint64_t len = parsing_reader->read();
    
    // Compute the AVL grammar expanding to phrase p.
    space_efficient_vector<pair_type> phrase_roots;
    if (len == 0) {

      // If this is a literal phrase, create a trivial grammar.
      const nonterminal_type root((char_type)pos);
      const std::uint64_t root_id = grammar->add_nonterminal(root);
      phrase_roots.push_back(pair_type(
            (text_offset_type)root_id,
            (text_offset_type)1));
    } else {

      // Self-overlapping phrases are unadressed in Rytter's paper.
      // We follow Theorem 6.1 in https://arxiv.org/abs/1910.10631v3.
      if (pos + len > prefix_length) {
        std::uint64_t left = len;
        std::uint64_t exist = prefix_length - pos;
        while (left > 0) {

          std::uint64_t next = std::min(left, exist);
          space_efficient_vector<pair_type> v;
          grammar->merge_enclosed_roots(pos, pos + next);
          grammar->decomposition(pos, pos + next, v);
          grammar->find_equivalent_seq(v);

          for (std::uint64_t t = 0; t < v.size(); ++t) {
            const std::uint64_t id = v[t].first;
            const std::uint64_t exp_len = v[t].second;
            prefix_length += exp_len;
            grammar->push_root(prefix_length, id);
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
        grammar->decomposition(begin, end, phrase_roots);
        grammar->find_equivalent_seq(phrase_roots);
      }
    }

    // Update prefix length and add new roots to the grammar.
    for (std::uint64_t t = 0; t < phrase_roots.size(); ++t) {
      const std::uint64_t id = phrase_roots[t].first;
      const std::uint64_t exp_len = phrase_roots[t].second;
      prefix_length += exp_len;
      grammar->push_root(prefix_length, id);
    }
  }

  long double total_time = utils::wclock() - start;
  fprintf(stderr, "\rInfo: elapsed = %.2Lfs, "
      "progress = %lu (100.00%%) phrases (%.2LfMiB prefix), "
      "peak RAM = %.2LfMiB",
      total_time, parsing_size, (1.L * prefix_length) / (1 << 20),
      (1.L * utils::get_peak_ram_allocation()) / (1UL << 20));

  // Clean up.
  parsing_reader->stop_reading();
  delete parsing_reader;

  // Return the result.
  return grammar;
}

#endif  // __CONVERT_LZ77_TO_AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
