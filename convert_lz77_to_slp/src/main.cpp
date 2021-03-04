#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <string>
#include <ctime>
#include <unistd.h>
#include <getopt.h>

#define MULTIROOT

#include "../include/types/uint40.hpp"
#include "../include/utils/utils.hpp"
#include "../include/utils/hash_table.hpp"
#include "../include/avl_grammar/avl_grammar_node.hpp"
#ifdef MULTIROOT
#include "../include/avl_grammar/avl_grammar_multiroot.hpp"
#include "../include/avl_grammar/convert_lz77_to_avl_grammar_multiroot.hpp"
#else
#include "../include/avl_grammar/avl_grammar.hpp"
#include "../include/avl_grammar/convert_lz77_to_avl_grammar.hpp"
#endif


template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
void test_conversion(
    std::string text_filename,
    std::string parsing_filename) {

  // Declare types.
  typedef std::pair<text_offset_type, text_offset_type> phrase_type;
#ifdef MULTIROOT
  typedef avl_grammar_multiroot<char_type> grammar_type;
#else
  typedef avl_grammar<char_type> grammar_type;
#endif

  // Turn paths absolute.
  text_filename = utils::absolute_path(text_filename);
  parsing_filename = utils::absolute_path(parsing_filename);

  // Obtain some basic statistics about input.
  std::uint64_t text_length =
    utils::file_size(text_filename) / sizeof(char_type);
  std::uint64_t parsing_size =
      utils::file_size(parsing_filename) / sizeof(phrase_type);

  // Print parameters.
  fprintf(stderr, "Convert LZ77 to SLP\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "Text filename = %s\n", text_filename.c_str());
  fprintf(stderr, "Parsing filename = %s\n", parsing_filename.c_str());
  fprintf(stderr, "Text length = %lu (%.2LfMiB)\n",
      text_length, (1.L * text_length * sizeof(char_type)) / (1 << 20));
  fprintf(stderr, "Number of LZ77 phrases = %lu\n", parsing_size);
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n",
      sizeof(text_offset_type));
  fprintf(stderr, "\n\n");

  // Read parsing.
  std::vector<phrase_type> parsing;
  {
    fprintf(stderr, "Read parsing... ");
    long double start = utils::wclock();
    parsing.resize(parsing_size);
    utils::read_from_file(parsing.data(), parsing_size, parsing_filename);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Convert LZ77 to AVL grammar.
  grammar_type *grammar = NULL;
  {
    fprintf(stderr, "Convert LZ77 to SLP... ");
    long double start = utils::wclock();
#ifdef MULTIROOT
    grammar =
      convert_lz77_to_avl_grammar_multiroot<char_type, text_offset_type>(parsing);
#else
    grammar =
      convert_lz77_to_avl_grammar<char_type, text_offset_type>(parsing);
#endif
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Print info. Note that the grammar may
  // still contain unused nonterminals.
  fprintf(stderr, "Grammar size = %lu\n", grammar->size());
  fprintf(stderr, "\n");


  // Run tests of correctness.
  fprintf(stderr, "Tests of correctness:\n");

  // Check if the resulting grammar indeed expands to the text.
  // For this, we first read the original text. Note: the text
  // could be streamed, or I could simple decoded it from LZ77.
  char_type *text = NULL;
  {
    fprintf(stderr, "  Read text... ");
    long double start = utils::wclock();
    text = new char_type[text_length];
    utils::read_from_file(text, text_length, text_filename);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Obtain the string to which the grammar is expanding.
  std::uint64_t decoded_text_length = 0;
  char_type *decoded_text = NULL;
  {
    fprintf(stderr, "  Decode text from grammar... ");
    long double start = utils::wclock();
    grammar->decode(decoded_text, decoded_text_length);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Compare text and decoded text.
  {
    fprintf(stderr, "  Compare the texts... ");
    long double start = utils::wclock();
    bool eq = true;
    if (text_length != decoded_text_length) eq = false;
    else {
      if (!std::equal(text, text + text_length, decoded_text))
        eq = false;
    }
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs ", elapsed);
    fprintf(stderr, "%s\n", eq ? "(OK)" : "(FAILED)");
  }

  // Test AVL property.
  {
    fprintf(stderr, "  Test AVL property... ");
    long double start = utils::wclock();
    bool result = grammar->test_avl_property();
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs ", elapsed);
    fprintf(stderr, "%s\n", result ? "(OK)" : "(FAILED)");
  }

  // Collect various statistic about grammar.
  fprintf(stderr, "\n");
  fprintf(stderr, "Additional statistics:\n");

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
    typedef avl_grammar_node<char_type> node_type;
    long double start = utils::wclock();
    std::vector<const node_type*> pointers;
    std::vector<std::uint64_t> hashes;
    grammar->collect_nonterminal_pointers(pointers);
    for (std::uint64_t i = 0; i < pointers.size(); ++i)
      hashes.push_back(pointers[i]->m_kr_hash);
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
    typedef avl_grammar_node<char_type> node_type;
    hash_table<const node_type*, std::uint64_t> hashes;
    hash_table<std::uint64_t, bool> seen_hashes;
    std::uint64_t count = 0;
    grammar->collect_mersenne_karp_rabin_hashes_2(hashes);
    grammar->count_nodes_in_pruned_grammar(hashes, seen_hashes, count);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs ", elapsed);
    fprintf(stderr, "(%lu)\n", count);
  }

  // Check the number of different reachable nonterminals.
  {
    fprintf(stderr, "  Collect reachable nonterminals... ");
    typedef avl_grammar_node<char_type> node_type;
    long double start = utils::wclock();
    std::vector<const node_type*> pointers;
    grammar->collect_nonterminal_pointers(pointers);
    std::sort(pointers.begin(), pointers.end());
    pointers.erase(std::unique(pointers.begin(), pointers.end()), pointers.end());
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs ", elapsed);
    fprintf(stderr, "(%lu)\n", pointers.size());
  }

#ifdef MULTIROOT
  {
    fprintf(stderr, "  Number of roots = %lu\n",
        grammar->m_roots.size());
  }
#endif

  // Clean up.
  delete[] text;
  delete[] decoded_text;
  for (std::uint64_t i = 0; i < grammar->m_nonterminals.size(); ++i)
    delete grammar->m_nonterminals[i];
}

int main(int argc, char **argv) {
  if (argc != 3)
    std::exit(EXIT_FAILURE);

  // Declare types.
  typedef std::uint8_t char_type;
  typedef uint40 text_offset_type;

  // Obtain filenames.
  std::string text_filename = argv[1];
  std::string parsing_filename = argv[2];

  // Run the algorithm.
  test_conversion<char_type, text_offset_type>(
      text_filename, parsing_filename);
}

