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

#include "../include/uint40.hpp"
#include "../include/utils.hpp"
#include "../include/hash_table.hpp"
#include "../include/avl_grammar_node.hpp"
#ifdef MULTIROOT
#include "../include/avl_grammar_multiroot.hpp"
#include "../include/convert_lz77_to_avl_grammar_multiroot.hpp"
#else
#include "../include/avl_grammar.hpp"
#include "../include/convert_lz77_to_avl_grammar.hpp"
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

  // Print parameters.
  fprintf(stderr, "Convert LZ77 to SLP\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "Text filename = %s\n", text_filename.c_str());
  fprintf(stderr, "Parsing filename = %s\n", parsing_filename.c_str());
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n",
      sizeof(text_offset_type));
  fprintf(stderr, "\n\n");

  // Read parsing.
  std::vector<phrase_type> parsing;
  {
    fprintf(stderr, "Read parsing... ");
    long double start = utils::wclock();
    std::uint64_t parsing_size =
      utils::file_size(parsing_filename) / sizeof(phrase_type);
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
  fprintf(stderr, "Number of phrases = %lu\n", parsing.size());
  fprintf(stderr, "Grammar size = %lu\n", grammar->size());

  // Check if the resulting grammar indeed expands to the text.
  // For this, we first read the original text. Note: the text
  // could be streamed, or I could simple decoded it from LZ77.
  std::uint64_t text_length = 0;
  char_type *text = NULL;
  {
    fprintf(stderr, "Read text... ");
    long double start = utils::wclock();
    text_length = utils::file_size(text_filename) / sizeof(char_type);
    text = new char_type[text_length];
    utils::read_from_file(text, text_length, text_filename);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Obtain the string to which the grammar is expanding.
  std::uint64_t decoded_text_length = 0;
  char_type *decoded_text = NULL;
  {
    fprintf(stderr, "Decode text from grammar... ");
    long double start = utils::wclock();
    grammar->decode(decoded_text, decoded_text_length);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lf\n", elapsed);
  }

  // Compare text and decoded text.
  {
    fprintf(stderr, "Compare texts... ");
    long double start = utils::wclock();
    bool eq = true;
    if (text_length != decoded_text_length) eq = false;
    else {
      if (!std::equal(text, text + text_length, decoded_text))
        eq = false;
    }
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
    fprintf(stderr, "Result: %s\n", eq ? "OK" : "ERROR");
  }

  // Test AVL property.
  {
    fprintf(stderr, "Test AVL property... ");
    long double start = utils::wclock();
    bool result = grammar->test_avl_property();
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
    fprintf(stderr, "AVL property = %s\n", result ? "TRUE" : "FALSE");
  }

#if 0
  // Check the number of different KR hashes.
  {
    fprintf(stderr, "Collect Karp-Rabin hashes... ");

    // To make sure there is no colission, we use three
    // sets of different Karp-Rabin hashes.
    const std::uint64_t p1 = 1000000007;
    const std::uint64_t p2 = 1000000009;
    const std::uint64_t p3 = 1000000013;
    const std::uint64_t a1 = 435725636;
    const std::uint64_t a2 = 947758375;
    const std::uint64_t a3 = 279584594;
    std::vector<std::uint64_t> hashes1;
    std::vector<std::uint64_t> hashes2;
    std::vector<std::uint64_t> hashes3;
    long double start = utils::wclock();
    grammar->collect_karp_rabin_hashes(hashes1, a1, p1);
    grammar->collect_karp_rabin_hashes(hashes2, a2, p2);
    grammar->collect_karp_rabin_hashes(hashes3, a3, p3);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
    std::vector<std::pair<std::pair<std::uint64_t, std::uint64_t>, std::uint64_t> > good_hashes;
    for (std::uint64_t i = 0; i < hashes1.size(); ++i)
      good_hashes.push_back(std::make_pair(std::make_pair(hashes1[i], hashes2[i]), hashes3[i]));
    std::sort(good_hashes.begin(), good_hashes.end());
    good_hashes.erase(std::unique(good_hashes.begin(), good_hashes.end()), good_hashes.end());
    fprintf(stderr, "Number of unique hashes = %lu\n", good_hashes.size());
  }
#endif

  // Check the number of different Mersenne KR hashes.
  // These hashes are much better, so there should no collisions.
  {
    fprintf(stderr, "Collect Mersenne Karp-Rabin hashes... ");
    std::vector<std::uint64_t> mersenne_hashes;
    long double start = utils::wclock();
    grammar->collect_mersenne_karp_rabin_hashes(mersenne_hashes);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
    std::sort(mersenne_hashes.begin(), mersenne_hashes.end());
    mersenne_hashes.erase(std::unique(mersenne_hashes.begin(),
          mersenne_hashes.end()), mersenne_hashes.end());
    fprintf(stderr, "Number of unique hashes = %lu\n",
        mersenne_hashes.size());
  }

  // Compute the number of nodes in the pruned grammar.
  {
    fprintf(stderr, "Count nodes in the pruned grammar... ");
    long double start = utils::wclock();
    typedef avl_grammar_node<char_type> node_type;
    hash_table<const node_type*, std::uint64_t> hashes;
    hash_table<std::uint64_t, bool> seen_hashes;
    std::uint64_t count = 0;
    grammar->collect_mersenne_karp_rabin_hashes_2(hashes);
    grammar->count_nodes_in_pruned_grammar(hashes, seen_hashes, count);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
    fprintf(stderr, "Number of nodes in the pruned grammar = %lu\n", count);
  }

  // Check the number of different reachable nonterminals.
  {
    fprintf(stderr, "Collect reachable nonterminals... ");
    typedef avl_grammar_node<char_type> node_type;
    long double start = utils::wclock();
    std::vector<const node_type*> pointers;
    grammar->collect_nonterminal_pointers(pointers);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
    std::sort(pointers.begin(), pointers.end());
    pointers.erase(std::unique(pointers.begin(), pointers.end()), pointers.end());
    fprintf(stderr, "Number of unique nonterminals = %lu\n", pointers.size());
  }

  // Check the number of different Mersenne hashes (method #2).
  {
    fprintf(stderr, "Collect Mersenne hashes (method #2)... ");
    typedef avl_grammar_node<char_type> node_type;
    long double start = utils::wclock();
    std::vector<const node_type*> pointers;
    std::vector<std::uint64_t> hashes;
    grammar->collect_nonterminal_pointers(pointers);
    for (std::uint64_t i = 0; i < pointers.size(); ++i)
      hashes.push_back(pointers[i]->m_kr_hash);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
    std::sort(hashes.begin(), hashes.end());
    hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
    fprintf(stderr, "Number of unique Mersenne hashes (method #2) = %lu\n", hashes.size());
  }

#ifdef MULTIROOT
  {
    fprintf(stderr, "Number of roots = %lu\n",
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

