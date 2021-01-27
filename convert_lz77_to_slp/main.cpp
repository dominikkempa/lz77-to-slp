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

  // Test AVL property.
  {
    fprintf(stderr, "Test AVL property... ");
    long double start = utils::wclock();
    bool result = grammar->test_avl_property();
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
    fprintf(stderr, "AVL property = %s\n", result ? "TRUE" : "FALSE");
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

