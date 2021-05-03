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
#include "../include/avl_grammar/avl_grammar.hpp"
#include "../include/avl_grammar/convert_lz77_to_avl_grammar.hpp"


//=============================================================================
// Decode text from streamed LZ77 parsing, and
// then compare to the output of grammar.
//=============================================================================
template<
  typename char_type,
  typename text_offset_type>
void check_correctness(
    const std::string parsing_filename,
    avl_grammar<char_type, text_offset_type> * const grammar) {

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
    avl_grammar<char_type, text_offset_type> * const grammar) {

  // Print initial message.
  fprintf(stderr, "\nRun additional tests:\n");

  // Check the number of different Mersenne KR hashes.
  // These hashes are much better, so there should no collisions.
  {
    fprintf(stderr, "  Compute different Karp-Rabin hashes... ");
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
    std::uint64_t total_rhs_length = 0;
    for (std::uint64_t i = 0; i < pointers.size(); ++i) {
      if (grammar->is_literal(pointers[i]))
        total_rhs_length += 1;
      else total_rhs_length += 2;
    }
    fprintf(stderr, "(#nonterminals = %lu, total rhs len = %lu), ",
        pointers.size(), total_rhs_length);
  }
}

template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
void test_conversion(
    std::string parsing_filename) {

  // Declare types.
  typedef std::pair<text_offset_type, text_offset_type> phrase_type;
  typedef avl_grammar<char_type, text_offset_type> grammar_type;

  // Turn paths absolute.
  parsing_filename = utils::absolute_path(parsing_filename);

  // Obtain some basic statistics about input.
  std::uint64_t parsing_size =
      utils::file_size(parsing_filename) / sizeof(phrase_type);

  // Print parameters.
  fprintf(stderr, "Convert LZ77 to SLP\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "Parsing filename = %s\n", parsing_filename.c_str());
  fprintf(stderr, "Number of LZ77 phrases = %lu (%.2LfMiB)\n",
      parsing_size, (2.L * parsing_size * sizeof(text_offset_type)) / (1 << 20));
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n", sizeof(text_offset_type));
  fprintf(stderr, "\n\n");

  // Convert LZ77 to AVL grammar.
  grammar_type *grammar = NULL;
  {
    fprintf(stderr, "Convert LZ77 to SLP...\n");
    long double start = utils::wclock();
    grammar =
      convert_lz77_to_avl_grammar<char_type, text_offset_type>(parsing_filename);

    // Print summary.
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "\nConversion time: %.2Lfs\n", elapsed);
  }

  // Obtain statistics.
  const std::uint64_t n_nonterminals = grammar->size();
  const std::uint64_t total_rhs_length = grammar->total_rhs_length();
  const std::uint64_t grammar_size = total_rhs_length;

  // Print info. Note that the grammar may
  // still contain unused nonterminals.
  fprintf(stderr, "Number of nonterminals = %lu\n", n_nonterminals);
  fprintf(stderr, "Grammar size = %lu\n", grammar_size);
  fprintf(stderr, "\n");

  // Print RAM use.
  // grammar->print_stats();

#ifdef TEST_CORRECTNESS
  check_correctness<char_type, text_offset_type>(parsing_filename, grammar);
#endif
#ifdef ADDITIONAL_TESTS
  additional_tests<char_type, text_offset_type>(grammar);
#endif

  // Clean up.
  delete grammar;

  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  RAM allocation: cur = %lu bytes, peak = %.2LfMiB\n",
      utils::get_current_ram_allocation(),
      (1.L * utils::get_peak_ram_allocation()) / (1UL << 20));
}

int main(int argc, char **argv) {
  if (argc != 2)
    std::exit(EXIT_FAILURE);

  // Initialize runtime statistics.
  utils::initialize_stats();

  // Declare types.
  typedef std::uint8_t char_type;
  typedef uint40 text_offset_type;

  // Obtain filenames.
  std::string parsing_filename = argv[1];

  // Run the algorithm.
  test_conversion<char_type, text_offset_type>(parsing_filename);
}

