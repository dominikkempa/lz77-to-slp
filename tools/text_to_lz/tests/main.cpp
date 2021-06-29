/**
 * @file    main.cpp
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

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <string>
#include <ctime>
#include <unistd.h>

#include "../include/uint40.hpp"
#include "../include/uint48.hpp"
#include "../include/utils.hpp"
#include "../include/compute_sa.hpp"
#include "../include/compute_lz77.hpp"
#include "sequential_lz77_decode.hpp"
#include "naive_lz77_encode.hpp"


//=============================================================================
// Run tests for text up to a given length on a given number of cases.
//=============================================================================
template<typename char_type, typename text_offset_type>
void test(
    const std::uint64_t max_text_length,
    const std::uint64_t testcases) {

  // Print initial message.
  fprintf(stderr, "TEST, max_length = %lu, testcases = %lu\n",
      max_text_length, testcases);

  // Declate types.
  typedef std::pair<text_offset_type, text_offset_type> pair_type;

  // Allocate text and SA.
  char_type * const text = new char_type[max_text_length];
  text_offset_type * const sa = new text_offset_type[max_text_length];

  // Run tests.
  for (std::uint64_t testid = 0; testid < testcases; ++testid) {

    // Print progress message.
    if (testid % 10 == 0)
      fprintf(stderr, "%.2Lf%%\r", (100.L * testid) / testcases);

    // Generate the text.
    const std::uint64_t text_length =
      utils::random_int<std::uint64_t>(
          (std::uint64_t)1,
          (std::uint64_t)max_text_length);
    for (std::uint64_t i = 0; i < text_length; ++i)
      text[i] = 'a' + utils::random_int<std::uint64_t>(0UL, 3UL);

    // Compute SA of text.
    compute_sa(text, text_length, sa);

    // Compute the parsin using the tested algorithm.
    std::vector<pair_type> parsing;
    compute_lz77::kkp2n(text, text_length, sa, parsing);

    // Compute the parsing using naive algorithm.
    std::vector<pair_type> correct_parsing;
    naive_lz77_encode(text, text_length, correct_parsing);

    // Compare parsing sizes.
    bool ok = true;
    if (parsing.size() != correct_parsing.size()) ok = false;
    else {
      for (std::uint64_t i = 0; i < parsing.size(); ++i)
        if (parsing[i].second != correct_parsing[i].second)
          ok = false;
    }

    // Decode computed parsing.
    char_type *decoded_text = sequential_lz77_decode
      <char_type, text_offset_type>(text_length, parsing);

    // Compare decoded text to the original.
    if (!ok || !std::equal(text, text + text_length, decoded_text)) {
      fprintf(stderr, "\nError:\n");
      fprintf(stderr, "  text_length = %lu\n", text_length);
      fprintf(stderr, "  text = ");
      for (std::uint64_t i = 0; i < text_length; ++i)
        fprintf(stderr, "%c", text[i]);
      fprintf(stderr, "\n");
      fprintf(stderr, "  correct parsing size = %lu\n",
          correct_parsing.size());
      fprintf(stderr, "  compute parsing size = %lu\n",
          parsing.size());
      std::exit(EXIT_FAILURE);
    }

    // Clean up.
    delete[] decoded_text;
  }

  // Clean up.
  delete[] text;
  delete[] sa;
}

int main() {

  // Init random number generator.
  srand(time(0) + getpid());

#ifdef NDEBUG
  static const std::uint64_t text_length_limit = (1 << 10);
  static const std::uint64_t n_tests = 10000;
#else
  static const std::uint64_t text_length_limit = (1 << 9);
  static const std::uint64_t n_tests = 1000;
#endif

  typedef std::uint8_t char_type;
  typedef uint48 text_offset_type;

  // Run tests.
  for (std::uint64_t max_text_length = 1;
      max_text_length <= text_length_limit; max_text_length *= 2)
    test<char_type, text_offset_type>(max_text_length, n_tests);

  // Print summary.
  fprintf(stderr, "All tests passed.\n");
}

