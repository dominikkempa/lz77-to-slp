#ifndef __CONVERT_LZ77_TO_AVL_GRAMMAR_HPP_INCLUDED
#define __CONVERT_LZ77_TO_AVL_GRAMMAR_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>

#include "../utils/karp_rabin_hashing.hpp"
#include "../io/async_stream_reader.hpp"
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
avl_grammar<char_type, text_offset_type>*
convert_lz77_to_avl_grammar(const std::string parsing_filename) {

  // Start the timer.
  long double start = utils::wclock();

  // Declare types.
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  typedef text_offset_type ptr_type;
  typedef avl_grammar<char_type, text_offset_type> grammar_type;
  typedef async_stream_reader<text_offset_type> reader_type;

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
  grammar_type *grammar = new grammar_type();
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

    // Get the next phrase len and pos values.
    std::uint64_t pos = parsing_reader->read();
    std::uint64_t len = parsing_reader->read();
    
    // Compute the AVL grammar expanding to phrase p.
    ptr_type phrase_root = std::numeric_limits<text_offset_type>::max();
    if (len == 0) {

      // If this is a literal phrase, create a trivial grammar.
      const nonterminal_type root((char_type)pos);
      const std::uint64_t root_id = grammar->add_nonterminal(root);
      phrase_root = root_id;
    } else {

      // Self-overlapping phrases are unadressed in Rytter's paper.
      // We follow Theorem 6.1 in https://arxiv.org/abs/1910.10631v3.
      if (pos + len > prefix_length) {

        // If the phase is self-overlapping, we create the
        // nonterminal expanding to text[pos..prefix_length).
        ptr_type suffix_nonterm =
          grammar->add_substring_nonterminal(pos, prefix_length);

        // Square the nonterminal until it reaches length >= len.
        ptr_type suffix_pow_nonterm = suffix_nonterm;
        std::uint64_t curlen = prefix_length - pos;
        while (curlen < len) {
          const ptr_type square =
              grammar->add_nonterminal(
                  suffix_pow_nonterm,
                  suffix_pow_nonterm);
          curlen <<= 1;
          suffix_pow_nonterm = square;
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
    const std::uint64_t exp_len = grammar->get_exp_len(phrase_root);
    prefix_length += exp_len;
    if (phrase_id == 0) {
      grammar->set_root(phrase_root);
    } else {
      const ptr_type old_root = grammar->get_root();
      const ptr_type new_root =
        grammar->add_concat_nonterminal(old_root, phrase_root);
      grammar->set_root(new_root);
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

#endif  // __CONVERT_LZ77_TO_AVL_GRAMMAR_HPP_INCLUDED
