#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>

#include "../include/utils.hpp"
#include "../include/compute_leftmost_lz77.hpp"


int main(int argc, char **argv) {

  // Exit if no filename is given.
  if (argc != 2)
    std::exit(EXIT_FAILURE);

  // Parse text filename and get text size.
  std::string text_filename = argv[1];
  std::uint64_t text_length = utils::file_size(text_filename);

  // Print info.
  fprintf(stderr, "Text filename = %s\n", text_filename.c_str());
  fprintf(stderr, "Text length = %lu\n", text_length);

  // Declare types.
  typedef std::uint8_t char_type;
  typedef std::uint32_t text_offset_type;

  // Read text.
  char_type *text = new char_type[text_length];
  utils::read_from_file(text, text_length, text_filename);

  // Parse.
  std::vector<std::pair<text_offset_type, text_offset_type> > parsing;
  {

    // Print initial message.
    fprintf(stderr, "Parsing: ");
    long double parsing_start = utils::wclock();

    // Compute parsing.
    compute_leftmost_lz77<text_offset_type>(text, text_length, parsing);

    // Print summary.
    long double parsing_time = utils::wclock() - parsing_start;
    fprintf(stderr, "%.2Lfs\n", parsing_time);
  }

  // Print summary.
  fprintf(stderr, "Number of phrases = %lu\n",
      (std::uint64_t)parsing.size());
}
