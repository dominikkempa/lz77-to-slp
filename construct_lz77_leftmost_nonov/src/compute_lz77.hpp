/**
 * @file    compute_lz77.hpp
 * @section LICENCE
 *
 * Copyright (C) 2017
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

#ifndef __COMPUTE_LZ77_HPP_INCLUDED
#define __COMPUTE_LZ77_HPP_INCLUDED

#include <cstdint>
#include <limits>
#include <algorithm>
#include <vector>

#include "compute_sa.hpp"
#include "byte_rank.hpp"
#include "rmq_tree.hpp"


//=============================================================================
// Excluding the output parsing, but including
// the input text, the function uses 3.128n +
// sizeof(text_offset_type) bytes, e.g., assuming
// 32-bit integers, it uses 7.128n bytes.
//=============================================================================
template<typename text_offset_type>
void compute_lz77(
    const std::uint8_t * const text,
    const std::uint64_t text_length,
    std::vector<std::pair<text_offset_type, text_offset_type> > &parsing) {

  // Handle special case.
  if (text_length == 0)
    return;

  // Check if text is non NULL.
  if (text == NULL) {
    fprintf(stderr, "\nError: text == NULL\n");
    std::exit(EXIT_FAILURE);
  }

  // Check if text_offset_type is sufficiently large.
  const std::uint64_t max_text_offset_type = 
    (std::uint64_t)std::numeric_limits<text_offset_type>::max();
  if (max_text_offset_type < text_length - 1) {
    fprintf(stderr, "\nError: text_offset_type is too small\n");
    std::exit(EXIT_FAILURE);
  }

  // Allocate SA.
  text_offset_type * const sa =
    new text_offset_type[text_length];

  // Allocate BWT.
  std::uint8_t * const bwt =
    new std::uint8_t[text_length];

  // Temporarily copy text to bwt array.
  std::uint8_t * const rev_text = bwt;
  std::copy(text, text + text_length, rev_text);
  std::reverse(rev_text, rev_text + text_length);
  
  // Compute SA of reversed string.
  compute_sa(rev_text, text_length, sa);
  for (std::uint64_t i = 0; i < text_length; ++i)
    sa[i] = text_length - 1 - (std::uint64_t)sa[i];

  // Compute BWT and isa0.
  bwt[0] = text[0];
  std::uint64_t isa0 = 0;
  {
    std::uint64_t ptr = 1;
    for (std::uint64_t j = 0; j < text_length; ++j)
      if ((std::uint64_t)sa[j] + 1 < text_length)
        bwt[ptr++] = text[(std::uint64_t)sa[j] + 1];
      else isa0 = j;
  }

  // Create RMQ over SA.
  typedef rmq_tree<text_offset_type> rmq_type;
  rmq_type * const rmq = new rmq_type(sa, text_length);

  // Compute the count array over BWT.
  static const std::uint64_t k_sigma = 256;
  std::vector<std::uint64_t> count(k_sigma, (std::uint64_t)0);
  for (std::uint64_t j = 0; j < text_length; ++j)
    ++count[(std::uint64_t)bwt[j]];

  // Do the exclusive prefix sum over count array.
  {
    std::uint64_t sum = 0;
    for (std::uint64_t j = 0; j < k_sigma; ++j) {
      const std::uint64_t temp = count[j];
      count[j] = sum;
      sum += temp;
    }
  }

  // Create rank over BWT.
  typedef byte_rank<> rank_type;
  rank_type * const rank =
    new rank_type(bwt, text_length);

  // Parse.
  std::uint64_t parsed_length = 0;
  while (parsed_length < text_length) {
    std::uint64_t next_phrase_length = 0;
    std::uint64_t beg = 0;
    std::uint64_t end = text_length;

    // Refine the range as long as it is nonempty and
    // it contains are lease one non-overlapping occurrence.
    while (parsed_length + next_phrase_length < text_length) {

      // Perform a single step of backwards search.
      // Invariant: current range in SA is [beg..end).
      const std::uint8_t c = text[parsed_length + next_phrase_length];
      const std::uint64_t newbeg = count[c] +
        rank->query(beg + (next_phrase_length > 0) - (beg > isa0), c);
      const std::uint64_t newend = count[c] +
        rank->query(end + 1 - (end > isa0), c);

      // Exit is the new range does not
      // fulfill the above conditions.
      if (newbeg >= newend ||
          !rmq->less(newbeg, newend, parsed_length))
        break;
      else {
        beg = newbeg;
        end = newend;
      }

      // Increase length of current phrase.
      ++next_phrase_length;
    }

    // Add the new phrase to the parsing.
    if (next_phrase_length == 0) {
      const text_offset_type pos = text[parsed_length];
      const text_offset_type len = (std::uint64_t)0;
      parsing.push_back(std::make_pair(pos, len));
      parsed_length += 1;
    } else {
      const std::uint64_t posidx = rmq->rmq(beg, end);
      const text_offset_type pos =
        (std::uint64_t)sa[posidx] + 1 - next_phrase_length;
      const text_offset_type len = next_phrase_length;
      parsing.push_back(std::make_pair(pos, len));
      parsed_length += next_phrase_length;
    }
  }

  // Clean up.
  delete rank;
  delete rmq;
  delete[] bwt;
  delete[] sa;
}

#endif  // __COMPUTE_LZ77_HPP_INCLUDED
