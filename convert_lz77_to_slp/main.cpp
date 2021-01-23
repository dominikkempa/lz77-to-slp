#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <string>
#include <ctime>
#include <unistd.h>
#include <getopt.h>

#include "include/uint40.hpp"
#include "include/utils.hpp"
#include "include/avl_grammar.hpp"


void quick_test_merge() {
  typedef std::uint8_t char_type;
  typedef avl_grammar_node<char_type> node_type;

  // Manually create an AVL grammar encoding the string
  // Fib_7 = abaababaabaab (the 7-th fibonacci string).
  // Nonterminals are names as in the example in Rytter's paper.
  std::vector<const node_type*> nonterminals;
  node_type *X1 = new node_type((char_type)'b');
  nonterminals.push_back(X1);
  node_type *X2 = new node_type((char_type)'a');
  nonterminals.push_back(X2);
  node_type *X3 = new node_type(X2, X1);
  nonterminals.push_back(X3);
  node_type *X4 = new node_type(X3, X2);
  nonterminals.push_back(X4);
  node_type *X5 = new node_type(X4, X3);
  nonterminals.push_back(X5);
  node_type *X6 = new node_type(X5, X4);
  nonterminals.push_back(X6);
  node_type *X7 = new node_type(X6, X5);
  nonterminals.push_back(X7);

  // Print expansions.
  for (std::uint64_t i = 0; i < nonterminals.size(); ++i) {
    fprintf(stderr, "exp(X%lu) = ", i + 1);
    nonterminals[i]->print_expansion();
    fprintf(stderr, ", height = %lu, explen = %lu\n",
        (std::uint64_t)nonterminals[i]->m_height,
        nonterminals[i]->m_exp_len);
  }

  // Test merging.
  {
    const node_type *X_7_6 =
      add_concat_nonterminal<char_type>(nonterminals, X7, X6);
    X_7_6->print_expansion();
    fprintf(stderr, "\n");
  }

  // Clean up.
  for (std::uint64_t i = 0; i < nonterminals.size(); ++i)
    delete nonterminals[i];
}

template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
void test_conversion(
    std::string text_filename,
    std::string parsing_filename) {

  // Declare types.
  typedef std::pair<text_offset_type, text_offset_type> phrase_type;
  typedef avl_grammar<char_type> grammar_type;

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

  // debug //
  // fprintf(stderr, "Parsing:\n");
  // for (std::uint64_t i = 0; i < parsing.size(); ++i)
  //   fprintf(stderr, "\t(%lu, %lu)\n",
  //       (std::uint64_t)parsing[i].first,
  //       (std::uint64_t)parsing[i].second);
  // fprintf(stderr, "\n");
  ///////////

  // Convert LZ77 to AVL grammar.
  grammar_type *grammar = NULL;
  {
    fprintf(stderr, "Convert LZ77 to SLP... ");
    long double start = utils::wclock();
    grammar =
      convert_lz77_to_avl_grammar<char_type, text_offset_type>(parsing);
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

