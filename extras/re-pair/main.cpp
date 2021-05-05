#include <cstdio>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <map>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <sys/time.h>
#include <unistd.h>

#include "aux/utils.hpp"
#include "repair.hpp"
#include "uint40.hpp"


template<
  typename char_type,
  typename text_offset_type>
void test_repair_on_file(std::string text_filename) {

  // Init basic stats.
  utils::initialize_stats();
  srand(time(0) + getpid());
  std::uint64_t text_length = utils::file_size(text_filename);

  // Turn paths absolute.
  text_filename = utils::absolute_path(text_filename);

  // Print basic info.
  fprintf(stderr, "Running Re-Pair algorithm\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "Text filename = %s\n", text_filename.c_str());
  fprintf(stderr, "Text length = %lu (%.2LfMiB)\n", text_length,
      (1.L * text_length * sizeof(char_type)) / (1UL << 20));
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n",
      sizeof(text_offset_type));

  // Allocate text.
  char_type *text = utils::allocate_array<char_type>(text_length);

  // Read text.
  {
    fprintf(stderr, "Reading text... ");
    utils::read_from_file<char_type>(text, text_length, text_filename);
    fprintf(stderr, "DONE\n");
  }

  // Run the algorithm.
  std::uint64_t grammar_size = 0;
  {
    fprintf(stderr, "Running algorithm...\n");
    long double start = utils::wclock();
    grammar_size = repair<char_type, text_offset_type>(text, text_length);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "Elapsed: %.2Lfs\n", elapsed);
  }

  // Clean up.
  utils::deallocate(text);

  // Print summary.
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "Grammar size = %lu\n", grammar_size);
  fprintf(stderr, "  RAM allocation: cur = %lu bytes, peak = %.2LfMiB\n",
      utils::get_current_ram_allocation(),
      (1.L * utils::get_peak_ram_allocation()) / (1UL << 20));
}

int main(int argc, char **argv) {
  if (argc < 2)
    std::exit(EXIT_FAILURE);
  
  // Declare types.
  typedef std::uint8_t char_type;
  typedef uint40 text_offset_type;

  // Parse text filename.
  std::string text_filename = std::string(argv[1]);

  // Run the algorithm.
  test_repair_on_file<char_type, text_offset_type>(text_filename);
}

