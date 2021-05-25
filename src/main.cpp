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
#include <getopt.h>

// #define TEST_CORRECTNESS
// #define ADDITIONAL_TESTS

#include "../include/types/uint40.hpp"
#include "../include/utils/utils.hpp"
#include "../include/lazy_avl_grammar/lazy_avl_grammar.hpp"
#include "../include/lazy_avl_grammar/lz_to_lazy_avl_grammar.hpp"


//=============================================================================
// Decode text from streamed LZ77 parsing, and
// then compare to the output of grammar.
//=============================================================================
template<
  typename char_type,
  typename text_offset_type>
void check_correctness(
    const std::string parsing_filename,
    lazy_avl_grammar<char_type, text_offset_type> * const grammar) {

  // Print initial message.
  fprintf(stderr, "\nRun correctness tests:\n");

  // Declare typedefs.
  typedef async_stream_reader<text_offset_type> reader_type;

  // Compute text length (needed to allocate text).
  std::uint64_t text_length = 0;
  std::uint64_t n_phrases = 0;
  {

    // Start the timer.
    fprintf(stderr, "  Compute text length... ");
    long double start = utils::wclock();

    // Initialize the parsing reader.
    const std::uint64_t bufsize = (1 << 19);
    const std::uint64_t n_buffers = 4;
    reader_type *parsing_reader = new reader_type(
        parsing_filename, bufsize, n_buffers);

    // Stream the parsing.
    while (!parsing_reader->empty()) {
      (void) parsing_reader->read();
      const std::uint64_t len = parsing_reader->read();
      std::uint64_t phrase_len = std::max((std::uint64_t)1, len);
      text_length += phrase_len;
      ++n_phrases;
    }

    // Clean up.
    parsing_reader->stop_reading();
    delete parsing_reader;

    // Print summary.
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
    fprintf(stderr, "  Text length = %lu (%.2LfMiB)\n",
        text_length, (1.L * text_length * sizeof(char_type) / (1 << 20)));
    fprintf(stderr, "  Average phrase length = %.2Lf\n",
        (1.L * text_length) / n_phrases);
  }

  // Allocate text.
  char_type * const text = utils::allocate_array<char_type>(text_length);

  // Decode text from streamed LZ77 parsing.
  {

    // Start the timer.
    fprintf(stderr, "  Decode text from LZ77... ");
    long double start = utils::wclock();

    // Initialize the parsing reader.
    const std::uint64_t bufsize = (1 << 19);
    const std::uint64_t n_buffers = 4;
    reader_type *parsing_reader = new reader_type(
        parsing_filename, bufsize, n_buffers);

    // Stream the parsing.
    std::uint64_t decoded_prefix_length = 0;
    while (!parsing_reader->empty()) {
      const std::uint64_t pos = parsing_reader->read();
      const std::uint64_t len = parsing_reader->read();
      if (len == 0)
        text[decoded_prefix_length++] = (char_type)pos;
      else {
        for (std::uint64_t i = 0; i < len; ++i)
          text[decoded_prefix_length++] = text[pos + i];
      }
    }

    // Clean up.
    parsing_reader->stop_reading();
    delete parsing_reader;

    // Print summary.
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Compare grammar output to text.
  {

    // Start the timer.
    fprintf(stderr, "  Compare grammar output to text... ");
    long double start = utils::wclock();

    // Run the comparison.
    bool eq = grammar->compare_expansion_to_text(text, text_length);

    // Print summary.
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lf ", elapsed);
    fprintf(stderr, "(%s)\n", eq ? "OK" : "FAILED");
  }

  // Test AVL property of the grammar.
  {

    // Start the timer.
    fprintf(stderr, "  Test AVL property... ");
    long double start = utils::wclock();

    // Run the test.
    bool result = grammar->test_avl_property();

    // Print summary.
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs ", elapsed);
    fprintf(stderr, "%s\n", result ? "(OK)" : "(FAILED)");
  }

  // Clean up.
  utils::deallocate(text);
}

template<
  typename char_type,
  typename text_offset_type>
void additional_tests(
    lazy_avl_grammar<char_type, text_offset_type> * const grammar) {

  // Print initial message.
  fprintf(stderr, "\nRun additional tests:\n");

  // Check the number of different Mersenne KR hashes.
  // These hashes are much better, so there should no collisions.
  {
    fprintf(stderr, "  Compute different Karp-Rabin hashes (method #1)... ");
    std::vector<std::uint64_t> mersenne_hashes;
    long double start = utils::wclock();
    grammar->collect_mersenne_karp_rabin_hashes(mersenne_hashes);
    std::sort(mersenne_hashes.begin(), mersenne_hashes.end());
    mersenne_hashes.erase(std::unique(mersenne_hashes.begin(),
          mersenne_hashes.end()), mersenne_hashes.end());
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs ", elapsed);
    fprintf(stderr, "(%lu)\n", mersenne_hashes.size());
  }

  // Check the number of different Mersenne hashes (method #2).
  {
    fprintf(stderr, "  Compute different Karp-Rabin hashes (method #2)... ");
    long double start = utils::wclock();
    std::vector<text_offset_type> pointers;
    std::vector<std::uint64_t> hashes;
    grammar->collect_nonterminal_pointers(pointers);
    for (std::uint64_t i = 0; i < pointers.size(); ++i)
      hashes.push_back(grammar->get_kr_hash(pointers[i]));
    std::sort(hashes.begin(), hashes.end());
    hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs ", elapsed);
    fprintf(stderr, "(%lu)\n", hashes.size());
  }

  // Compute the number of nodes in the pruned grammar.
  {
    fprintf(stderr, "  Count nodes in the pruned grammar... ");
    long double start = utils::wclock();
    hash_table<text_offset_type, std::uint64_t> hashes;
    hash_table<std::uint64_t, bool> seen_hashes;
    std::uint64_t count = 0;
    grammar->collect_mersenne_karp_rabin_hashes_2(hashes);
    grammar->count_nonterminals_in_pruned_grammar(hashes, seen_hashes, count);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs ", elapsed);
    fprintf(stderr, "(%lu)\n", count);
  }

  // Check the number of different reachable nonterminals.
  {
    fprintf(stderr, "  Collect reachable nonterminals... ");
    long double start = utils::wclock();
    std::vector<text_offset_type> pointers;
    grammar->collect_nonterminal_pointers(pointers);
    std::sort(pointers.begin(), pointers.end());
    pointers.erase(std::unique(pointers.begin(),
          pointers.end()), pointers.end());
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs ", elapsed);
    fprintf(stderr, "(%lu)\n", pointers.size());
  }
}

template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
void test_conversion(
    std::string parsing_filename,
    std::string output_filename,
    bool use_kr_hashing,
    long double kr_hashing_prob) {

  // Declare types.
  typedef std::pair<text_offset_type, text_offset_type> phrase_type;
  typedef lazy_avl_grammar<char_type, text_offset_type> grammar_type;

  // Turn paths absolute.
  parsing_filename = utils::absolute_path(parsing_filename);
  output_filename = utils::absolute_path(output_filename);

  // Obtain some basic statistics about input.
  std::uint64_t parsing_size =
      utils::file_size(parsing_filename) / sizeof(phrase_type);

  // Print parameters.
  fprintf(stderr, "Convert LZ77 to SLP\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "Parsing filename = %s\n", parsing_filename.c_str());
  fprintf(stderr, "Output filename = %s\n", output_filename.c_str());
  fprintf(stderr, "Number of LZ77 phrases = %lu (%.2LfMiB)\n",
      parsing_size, (2.L * parsing_size * sizeof(text_offset_type)) / (1 << 20));
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n", sizeof(text_offset_type));
  fprintf(stderr, "Use KR hashing = %s\n", use_kr_hashing ? "TRUE" : "FALSE");
  if (use_kr_hashing)
    fprintf(stderr, "KR hashing probability = %.2Lf%%\n",
        100.L * kr_hashing_prob);
  fprintf(stderr, "\n\n");

  // Convert LZ77 to AVL grammar.
  grammar_type *grammar = NULL;
  std::uint64_t text_length = 0;
  std::uint64_t n_phrases = 0;
  long double conversion_time = 0.0;
  {
    fprintf(stderr, "Convert LZ77 to SLP...\n");
    long double start = utils::wclock();
    grammar =
      lz_to_lazy_avl_grammar<char_type, text_offset_type>(
          parsing_filename, use_kr_hashing, kr_hashing_prob,
          n_phrases, text_length);

    // Print summary.
    conversion_time = utils::wclock() - start;
    fprintf(stderr, "\n\n");
  }

  // Obtain statistics.
  const std::uint64_t n_nonterminals = grammar->size();
  const std::uint64_t n_roots = grammar->number_of_roots();
  const std::uint64_t total_rhs_length = grammar->total_rhs_length();
  const std::uint64_t grammar_size = total_rhs_length + n_roots;
  const std::uint64_t ram_use = grammar->ram_use();
  long double avoided_merges = 100.L * grammar->get_avoided_merges();

  // Print info. Note that the grammar may
  // still contain unused nonterminals.
  fprintf(stderr, "Statistics:\n");
  fprintf(stderr, "  Text length = %lu\n", text_length);
  fprintf(stderr, "  Number of nonterminals = %lu\n", n_nonterminals);
  fprintf(stderr, "  Number of roots = %lu\n", n_roots);
  fprintf(stderr, "  Grammar size = %lu (%.2Lfelems/phrase)\n",
      grammar_size, (1.L * grammar_size) / n_phrases);
  fprintf(stderr, "  Conversion time = %.2Lfs (%.2Lfns/char)\n",
      conversion_time, (1000000000.L * conversion_time) / text_length);
  fprintf(stderr, "  Grammar RAM use = %.2LfMiB (%.2Lfbytes/phrase)\n",
      (1.L * ram_use) / (1UL << 20), (1.L * ram_use) / n_phrases);
  fprintf(stderr, "  Avoided merges: %.2Lf%%\n", avoided_merges);
  fprintf(stderr, "\n");

  // Store stats useful for experiments.
#if 0
  {
    std::string stats_filename = parsing_filename + ".stats";
    std::FILE *f = std::fopen(stats_filename.c_str(), "a+");
    fprintf(f, "GramSiz\t%.2Lf\t%.2Lf\n",
        kr_hashing_prob, (1.L * grammar_size) / n_phrases);
    fprintf(f, "Runtime\t%.2Lf\t%.2Lf\n",
        kr_hashing_prob, (1000000000.L * conversion_time) / text_length);
    fprintf(f, "RamUse\t%.2Lf\t%.2Lf\n",
        kr_hashing_prob, (1.L * ram_use) / n_phrases);
    fprintf(f, "AvMerge\t%.2Lf\t%.2Lf\n",
        kr_hashing_prob, avoided_merges);
    std::fclose(f);
  }
#endif

  // Print RAM use.
  grammar->print_stats();

#ifdef TEST_CORRECTNESS
  check_correctness<char_type, text_offset_type>(parsing_filename, grammar);
#endif
#ifdef ADDITIONAL_TESTS
  additional_tests<char_type, text_offset_type>(grammar);
#endif

  // Write grammar to file.
  {
    fprintf(stderr, "\nWrite grammar to file... ");
    long double start = utils::wclock();
    grammar->write_to_file(output_filename);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "DONE (%.2Lfs)\n", elapsed);
  }

  // Clean up.
  delete grammar;

  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  RAM allocation: cur = %lu bytes, peak = %.2LfMiB\n",
      utils::get_current_ram_allocation(),
      (1.L * utils::get_peak_ram_allocation()) / (1UL << 20));
}

int main(int argc, char **argv) {
  if (argc != 2 && argc != 3)
    std::exit(EXIT_FAILURE);

  // Initialize runtime statistics.
  utils::initialize_stats();

  // Declare types.
  typedef std::uint8_t char_type;
  typedef uint40 text_offset_type;

  // Obtain filenames.
  std::string parsing_filename = argv[1];
  std::string output_filename = parsing_filename + ".slg";
  bool use_kr_hashing = true;
  long double kr_hashing_prob = 0.125;

  if (argc == 3) {
    std::string p = std::string(argv[2]);
    if (p == std::string("0")) {
      use_kr_hashing = false;
      kr_hashing_prob = 0.0;
    } else
      kr_hashing_prob = atof(p.c_str());
  }

  // Run the algorithm.
  test_conversion<char_type, text_offset_type>(
      parsing_filename, output_filename,
      use_kr_hashing, kr_hashing_prob);
}

