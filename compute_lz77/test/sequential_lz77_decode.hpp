#ifndef __SEQUENTIAL_LZ77_DECODE_HPP_INCLUDED
#define __SEQUENTIAL_LZ77_DECODE_HPP_INCLUDED

#include <vector>
#include <algorithm>
#include <cstdint>


//=============================================================================
// Sequential LZ77 decoding.
//=============================================================================
template<
  typename char_type,
  typename text_offset_type>
char_type *sequential_lz77_decode(
    const std::uint64_t text_length,
    std::vector<std::pair<text_offset_type,
      text_offset_type> > &parsing) {

  // Get parsing size.
  const std::uint64_t n_phrases = parsing.size();

  // Allocate the text.
  char_type * const text = new char_type[text_length];

  // Decode the text.
  std::uint64_t j = 0;
  for (std::uint64_t i = 0; i < n_phrases; ++i) {
    const std::uint64_t pos = parsing[i].first;
    const std::uint64_t len = parsing[i].second;
    if (len == 0) text[j++] = (char_type)pos;
    else {
      for (std::uint64_t t = 0; t < len; ++t)
        text[j + t] = text[pos + t];
      j += len;
    }
  }

  // Return the result.
  return text;
}

#endif  // __SEQUENTIAL_LZ77_DECODE_HPP_INCLUDED
