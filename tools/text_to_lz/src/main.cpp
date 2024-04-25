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

#include "../include/types/uint40.hpp"
#include "../include/utils/utils.hpp"
#include "../include/utils/compute_sa.hpp"
#include "../include/compute_lz77.hpp"


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
  const std::uint64_t text_length = utils::file_size(text_filename);

  // Print parameters.
  fprintf(stderr, "Construct LZ77 parsing\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "Text filename = %s\n", text_filename.c_str());
  fprintf(stderr, "Output filename = %s\n", output_filename.c_str());
  fprintf(stderr, "Text length = %lu\n", text_length);
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n",
      sizeof(text_offset_type));
  fprintf(stderr, "\n\n");

  // Allocate text.
  char_type * const text = new char_type[text_length];

  // Read text.
  {
    fprintf(stderr, "Read text... ");
    long double start = utils::wclock();
    utils::read_from_file(text, text_length, text_filename);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Allocate SA.
  text_offset_type * const sa = new text_offset_type[text_length];

  // Compute SA.
  {
    fprintf(stderr, "Compute SA... ");
    const long double start = utils::wclock();
    compute_sa(text, text_length, sa);
    const long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Compute parsing.
  std::vector<pair_type> parsing;
  {
    fprintf(stderr, "Compute LZ77... ");
    const long double start = utils::wclock();
    compute_lz77::kkp2n(text, text_length, sa, parsing);
    const long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Print parsing size.
  fprintf(stderr, "Parsing size = %lu\n", parsing.size());

  // Write parsing to file.
  {
    fprintf(stderr, "Write parsing to file... ");
    const long double start = utils::wclock();
    const pair_type * const parsing_data = parsing.data();
    utils::write_to_file<pair_type>(parsing_data,
        parsing.size(), output_filename);
    const long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Clean up.
  delete[] sa;
  delete[] text;

  // Print final message.
  fprintf(stderr, "Computation finished\n");
}

//=============================================================================
// Print usage instructions and exit.
//=============================================================================
void usage(
    const char * const program_name,
    const int status) {
  printf(

"Usage: %s [OPTION]... FILE\n"
"Construct the LZ77 parsing of text stored in FILE.\n"
"\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -h, --help              display this help and exit\n"
"  -o, --output=OUTFILE    specify output filename. Default: FILE.lz77\n",

    program_name);

  std::exit(status);
}

int main(int argc, char **argv) {

  // Initial setup.
  srand(time(0) + getpid());
  const char * const program_name = argv[0];

  // Declare flags.
  static struct option long_options[] = {
    {"help",     no_argument,       NULL, 'h'},
    {"output",   required_argument, NULL, 'o'},
    {NULL,       0,                 NULL, 0}
  };

  // Initialize output filename.
  std::string output_filename("");

  // Parse command-line options.
  int c;
  while ((c = getopt_long(argc, argv, "ho:",
          long_options, NULL)) != -1) {
    switch(c) {
      case 'h':
        usage(program_name, EXIT_FAILURE);
        break;
      case 'o':
        output_filename = std::string(optarg);
        break;
      default:
        usage(program_name, EXIT_FAILURE);
        break;
    }
  }

  // Print error if there is not file.
  if (optind >= argc) {
    fprintf(stderr, "Error: FILE not provided\n\n");
    usage(program_name, EXIT_FAILURE);
  }

  // Parse the text filename.
  const std::string text_filename = std::string(argv[optind++]);
  if (optind < argc) {
    fprintf(stderr, "Warning: multiple input files provided. "
    "Only the first will be processed.\n");
  }

  // Set default output filename (if not provided).
  if (output_filename.empty())
    output_filename = text_filename + ".lz77";

  // Check for the existence of text.
  if (!utils::file_exists(text_filename)) {
    fprintf(stderr, "Error: input file (%s) does not exist\n\n",
        text_filename.c_str());
    usage(program_name, EXIT_FAILURE);
  }

  // Check if output file exists.
  if (utils::file_exists(output_filename)) {

    // Output file exists, should we proceed?
    char *line = new char[80];
    std::uint64_t len = 0;

    // Obtain the answer.
    do {
      printf("Output file (%s) exists. Overwrite? [y/n]: ",
          output_filename.c_str());
      char ch = 0;
      for (len = 0L; ((ch = std::getc(stdin)) != EOF) && (ch != '\n'); ++len)
        line[len] = ch;
    } while (len != 1 || (line[0] != 'y' && line[0] != 'n'));

    // If not, then exit.
    if (line[0] == 'n') {
      delete[] line;
      std::exit(EXIT_FAILURE);
    }

    // Otherwise, we proceed.
    delete[] line;
  }

  // Set types.
  typedef std::uint8_t char_type;
  typedef uint40 text_offset_type;

  // Run the parsing algorithm.
  compute_lz77_and_write_to_file<char_type, text_offset_type>(
      text_filename, output_filename);
}

