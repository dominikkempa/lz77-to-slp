#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <string>
#include <ctime>
#include <unistd.h>

#include "include/utils.hpp"
#include "include/compute_sa.hpp"
#include "include/compute_lz77.hpp"


//=============================================================================
// Computes the LZ77 parsing of the file given
// as an argument and write to output file.
//=============================================================================
template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
void compute_lz77_and_write_to_file(
    std::string text_filename,
    std::string output_filename) {

  // Declare typedefs.
  typedef std::pair<text_offset_type, text_offset_type> pair_type;

  // Turn paths absolute.                                                       
  text_filename = utils::absolute_path(text_filename);                          
  output_filename = utils::absolute_path(output_filename);

  // Get filesize.
  std::uint64_t text_length = utils::file_size(text_filename);

  // Print parameters.
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "Text filename = %s\n", text_filename.c_str());
  fprintf(stderr, "Output filename = %s\n", output_filename.c_str());
  fprintf(stderr, "Text length = %lu\n", text_length);
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n",
      sizeof(text_offset_type));
  fprintf(stderr, "\n\n");

  // Allocate text.
  char_type *text = new char_type[text_length];

  // Read text.
  {
    fprintf(stderr, "Read text... ");
    long double start = utils::wclock();
    utils::read_from_file(text, text_length, text_filename);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Allocate SA.
  text_offset_type *sa = new text_offset_type[text_length];

  // Compute SA.
  {
    fprintf(stderr, "Compute SA... ");
    long double start = utils::wclock();
    compute_sa(text, text_length, sa);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Compute parsing.
  std::vector<pair_type> parsing;
  {
    fprintf(stderr, "Compute LZ77... ");
    long double start = utils::wclock();
    compute_lz77::kkp2n(text, text_length, sa, parsing);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Print parsing size.
  fprintf(stderr, "Parsing size = %lu\n", parsing.size());

  // Write parsing to file.
  {
    fprintf(stderr, "Write parsing to file... ");
    long double start = utils::wclock();
    const pair_type * const parsing_data = parsing.data();
    utils::write_to_file<pair_type>(parsing_data,
        parsing.size(), output_filename);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Clean up.
  delete[] sa;
  delete[] text;

  // Print final message.
  fprintf(stderr, "Computation finished\n");
}
    

int main(int argc, char **argv) {

  // Check input parameters.
  if (argc != 3)
    std::exit(EXIT_FAILURE);

  // Set types.
  typedef std::uint8_t char_type;
  typedef uint40 text_offset_type;

  // Get filename and filesize.
  std::string text_filename = argv[1];
  std::string output_filename = argv[2];

  // Run the parsing algorithm.
  compute_lz77_and_write_to_file<char_type, text_offset_type>(
      text_filename, output_filename);
}

