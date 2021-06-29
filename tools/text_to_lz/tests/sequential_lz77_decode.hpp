/**
 * @file    sequential_lz77_decode.hpp
 * @section LICENCE
 *
 * This file is part of Lazy-AVLG v0.1.0
 * See: https://github.com/dominikkempa/lz77-to-slp
 *
 * Copyright (C) 2021
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __SEQUENTIAL_LZ77_DECODE_HPP_INCLUDED
#define __SEQUENTIAL_LZ77_DECODE_HPP_INCLUDED

#include <cstdint>
#include <vector>
#include <algorithm>


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
