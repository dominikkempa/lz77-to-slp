#ifndef __NAIVE_LZ77_ENCODE_HPP_INCLUDED
#define __NAIVE_LZ77_ENCODE_HPP_INCLUDED

#include <cstdint>
#include <vector>
#include <algorithm>


//=============================================================================
// Compute the LZ77 parsing naively in O(n^2) time.
//=============================================================================
template<
  typename char_type,
  typename text_offset_type>
void naive_lz77_encode(
    const char_type * const text,
    const std::uint64_t text_length,
    std::vector<std::pair<text_offset_type,
      text_offset_type> > &parsing) {

  // Compute the parsing.
  std::uint64_t cur_pos = 0;
  while (cur_pos < text_length) {
    std::uint64_t max_lcp = 0;
    std::uint64_t prev_pos;

    // Try all positions for previous factor.
    for (std::uint64_t pos = 0; pos < cur_pos; ++pos) {
      std::uint64_t lcp = 0;
      while (cur_pos + lcp < text_length &&
          text[pos + lcp] == text[cur_pos + lcp])
        ++lcp;

      // Update max previous factor.
      if (lcp > max_lcp) {
        max_lcp = lcp;
        prev_pos = pos;
      }
    }

    // Store the phrase.
    if (max_lcp == 0) {
      parsing.push_back(
          std::make_pair(
            (text_offset_type)text[cur_pos],
            (text_offset_type)0));
      ++cur_pos;
    } else {
      parsing.push_back(
          std::make_pair(
            (text_offset_type)prev_pos,
            (text_offset_type)max_lcp));
      cur_pos += max_lcp;
    }
  }
}

#endif  // __NAIVE_LZ77_ENCODE_HPP_INCLUDED
