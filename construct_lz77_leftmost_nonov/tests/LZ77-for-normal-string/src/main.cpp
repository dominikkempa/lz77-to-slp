#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <ctime>
#include <unistd.h>

#include "../include/compute_leftmost_lz77.hpp"


//=============================================================================
// Compute leftmost non-verlapping LZ77 (pair version) naively.
//=============================================================================
template<typename text_offset_type>
void naive_leftmost_lz77_nonov(
    const std::uint8_t * const text,
    const std::uint64_t text_length,
    std::vector<std::pair<text_offset_type,
        text_offset_type> > &parsing) {

  // Parse.
  std::uint64_t parsed_length = 0;
  while (parsed_length < text_length) {
    std::uint64_t max_lcp = 0;
    std::uint64_t max_lcp_pos = 0;
    for (std::uint64_t j = 0; j < parsed_length; ++j) {
      std::uint64_t lcp = 0;
      while (j + lcp < parsed_length &&
          parsed_length + lcp < text_length &&
          text[j + lcp] == text[parsed_length + lcp])
        ++lcp;
      if (lcp > max_lcp) {
        max_lcp = lcp;
        max_lcp_pos = j;
      }
    }

    if (max_lcp == 0) {
      parsing.push_back(std::make_pair(
            text[parsed_length],
            (text_offset_type)((std::uint64_t)0)));
      ++parsed_length;
    } else {
      parsing.push_back(std::make_pair(
            (text_offset_type)max_lcp_pos,
            (text_offset_type)max_lcp));
      parsed_length += max_lcp;
    }
  }
}

//=============================================================================
// Compute leftmost LZ77 (pair version) naively.
//=============================================================================
template<typename text_offset_type>
void naive_leftmost_lz77(
    const std::uint8_t * const text,
    const std::uint64_t text_length,
    std::vector<std::pair<text_offset_type,
        text_offset_type> > &parsing) {

  // Parse.
  std::uint64_t parsed_length = 0;
  while (parsed_length < text_length) {
    std::uint64_t max_lcp = 0;
    std::uint64_t max_lcp_pos = 0;
    for (std::uint64_t j = 0; j < parsed_length; ++j) {
      std::uint64_t lcp = 0;
      while (parsed_length + lcp < text_length &&
          text[j + lcp] == text[parsed_length + lcp])
        ++lcp;
      if (lcp > max_lcp) {
        max_lcp = lcp;
        max_lcp_pos = j;
      }
    }

    if (max_lcp == 0) {
      parsing.push_back(std::make_pair(
            text[parsed_length],
            (text_offset_type)((std::uint64_t)0)));
      ++parsed_length;
    } else {
      parsing.push_back(std::make_pair(
            (text_offset_type)max_lcp_pos,
            (text_offset_type)max_lcp));
      parsed_length += max_lcp;
    }
  }
}

//=============================================================================
// Verify the correctness of the (possibly overlapping) parsing.
//=============================================================================
template<typename text_offset_type>
bool verify_parsing(
    const std::uint8_t *const text,
    const std::uint64_t text_length,
    const std::vector<std::pair<text_offset_type,
        text_offset_type> > &parsing) {

  // Allocate array for decoded text.
  std::uint8_t * const decoded_text =
    new std::uint8_t[text_length];
  std::uint64_t decoded_length = 0;

  // Decode parsing.
  for (std::uint64_t j = 0; j < parsing.size(); ++j) {
    const std::uint64_t len = parsing[j].second;
    if (len == 0) {

      // Return false if decoded text is too long.
      if (decoded_length == text_length)
        return false;

      // Decode next symbol.
      decoded_text[decoded_length++] =
        (std::uint8_t)parsing[j].first;
    } else {

      // Return false, if the next phrase
      // would create a text that is too long.
      if (decoded_length + len > text_length)
        return false;

      // Obtain the pointer to the source.
      const std::uint64_t pos = parsing[j].first;

      // Decode next phrase.
      for (std::uint64_t i = 0; i < len; ++i)
        decoded_text[decoded_length++] =
          decoded_text[pos + i];
    }
  }

  // Return false if the text was too short.
  if (decoded_length != text_length)
    return false;

  // Return false if decoded text equal from original.
  if (!std::equal(text, text + text_length, decoded_text))
    return false;

  // Clean up.
  delete[] decoded_text;

  // Return the result.
  return true;
}

//=============================================================================
// Verify the correctness of the non-overlapping parsing.
//=============================================================================
template<typename text_offset_type>
bool verify_parsing_nonov(
    const std::uint8_t *const text,
    const std::uint64_t text_length,
    const std::vector<std::pair<text_offset_type,
        text_offset_type> > &parsing) {

  // Allocate array for decoded text.
  std::uint8_t * const decoded_text =
    new std::uint8_t[text_length];
  std::uint64_t decoded_length = 0;

  // Decode parsing.
  for (std::uint64_t j = 0; j < parsing.size(); ++j) {
    const std::uint64_t len = parsing[j].second;
    if (len == 0) {

      // Return false if decoded text is too long.
      if (decoded_length == text_length)
        return false;

      // Decode next symbol.
      decoded_text[decoded_length++] =
        (std::uint8_t)parsing[j].first;
    } else {

      // Return false, if the next phrase
      // would create a text that is too long.
      if (decoded_length + len > text_length)
        return false;

      // Obtain the pointer to the source.
      const std::uint64_t pos = parsing[j].first;

      // Return false, if the source overlaps the phrase.
      if (pos + len > decoded_length)
        return false;

      // Decode next phrase.
      for (std::uint64_t i = 0; i < len; ++i)
        decoded_text[decoded_length++] =
          decoded_text[pos + i];
    }
  }

  // Return false if the text was too short.
  if (decoded_length != text_length)
    return false;

  // Return false if decoded text equal from original.
  if (!std::equal(text, text + text_length, decoded_text))
    return false;

  // Clean up.
  delete[] decoded_text;

  // Return the result.
  return true;
}

//=============================================================================
// Test the leftmost non-overlapping LZ77 (pair version).
//=============================================================================
template<typename text_offset_type>
void test_leftmost_nonov_pair(
    const std::uint8_t * const text,
    const std::uint64_t text_length) {

  // Compute parsing using the algorithm.
  std::vector<std::pair<text_offset_type, text_offset_type> >
    parsing_computed;
  compute_leftmost_lz77_nonov<text_offset_type>(
      text, text_length, parsing_computed);

  // Compute correct parsing.
  std::vector<std::pair<text_offset_type, text_offset_type> >
    parsing_correct;
  naive_leftmost_lz77_nonov<text_offset_type>(
      text, text_length, parsing_correct);

  // Compare results.
  if (parsing_correct.size() != parsing_computed.size() ||
      !verify_parsing_nonov<text_offset_type>(text,
        text_length, parsing_computed)) {

    // Print debug info.
    fprintf(stderr, "\nError:\n");
    fprintf(stderr, "  text = ");
    for (std::uint64_t j = 0; j < text_length; ++j)
      fprintf(stderr, "%c", text[j]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  text length = %lu\n", text_length);
    fprintf(stderr, "  correct parsing: ");
    for (std::uint64_t j = 0; j < parsing_correct.size(); ++j) {
      if ((std::uint64_t)parsing_correct[j].second == 0)
        fprintf(stderr, "(%c, %lu), ",
            (std::uint8_t)parsing_correct[j].first,
            (std::uint64_t)parsing_correct[j].second);
      else
        fprintf(stderr, "(%lu, %lu), ",
            (std::uint64_t)parsing_correct[j].first,
            (std::uint64_t)parsing_correct[j].second);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "  computed parsing: ");
    for (std::uint64_t j = 0; j < parsing_computed.size(); ++j) {
      if ((std::uint64_t)parsing_computed[j].second == 0)
        fprintf(stderr, "(%c, %lu), ",
            (std::uint8_t)parsing_computed[j].first,
            (std::uint64_t)parsing_computed[j].second);
      else
        fprintf(stderr, "(%lu, %lu), ",
            (std::uint64_t)parsing_computed[j].first,
            (std::uint64_t)parsing_computed[j].second);
    }
    fprintf(stderr, "\n");
    std::exit(EXIT_FAILURE);
  }
}

//=============================================================================
// Test the leftmost LZ77 (pair version).
//=============================================================================
template<typename text_offset_type>
void test_leftmost_pair(
    const std::uint8_t * const text,
    const std::uint64_t text_length) {

  // Compute parsing using the algorithm.
  std::vector<std::pair<text_offset_type, text_offset_type> >
    parsing_computed;
  compute_leftmost_lz77<text_offset_type>(
      text, text_length, parsing_computed);

  // Compute correct parsing.
  std::vector<std::pair<text_offset_type, text_offset_type> >
    parsing_correct;
  naive_leftmost_lz77<text_offset_type>(
      text, text_length, parsing_correct);

  // Compare results.
  if (parsing_correct.size() != parsing_computed.size() ||
      !verify_parsing<text_offset_type>(text,
        text_length, parsing_computed)) {

    // Print debug info.
    fprintf(stderr, "\nError:\n");
    fprintf(stderr, "  text = ");
    for (std::uint64_t j = 0; j < text_length; ++j)
      fprintf(stderr, "%c", text[j]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  text length = %lu\n", text_length);
    fprintf(stderr, "  correct parsing: ");
    for (std::uint64_t j = 0; j < parsing_correct.size(); ++j) {
      if ((std::uint64_t)parsing_correct[j].second == 0)
        fprintf(stderr, "(%c, %lu), ",
            (std::uint8_t)parsing_correct[j].first,
            (std::uint64_t)parsing_correct[j].second);
      else
        fprintf(stderr, "(%lu, %lu), ",
            (std::uint64_t)parsing_correct[j].first,
            (std::uint64_t)parsing_correct[j].second);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "  computed parsing: ");
    for (std::uint64_t j = 0; j < parsing_computed.size(); ++j) {
      if ((std::uint64_t)parsing_computed[j].second == 0)
        fprintf(stderr, "(%c, %lu), ",
            (std::uint8_t)parsing_computed[j].first,
            (std::uint64_t)parsing_computed[j].second);
      else
        fprintf(stderr, "(%lu, %lu), ",
            (std::uint64_t)parsing_computed[j].first,
            (std::uint64_t)parsing_computed[j].second);
    }
    fprintf(stderr, "\n");
    std::exit(EXIT_FAILURE);
  }
}

int main() {

  // Init rand.
  srand(time(0) + getpid());

  // Declare types.
  typedef std::uint8_t char_type;
  typedef std::uint32_t text_offset_type;

  {

    // Print initial message.
    fprintf(stderr, "Alphabet size = 2\n");

    // Allocate text.
#ifdef NDEBUG
    static const std::uint64_t max_text_length = 17;
#else
    static const std::uint64_t max_text_length = 11;
#endif
    char_type * const text = new char_type[max_text_length];

    // Run tests.
    for (std::uint64_t text_length = 1;
        text_length <= max_text_length; ++text_length) {

      // Print text length.
      fprintf(stderr, "Text length = %lu\n", text_length);

      // Run tests.
      const std::uint64_t n_texts = ((std::uint64_t)1 << text_length);
      for (std::uint64_t text_id = 0; text_id < n_texts; ++text_id) {
        for (std::uint64_t j = 0; j < text_length; ++j)
          if (text_id & ((std::uint64_t)1 << j)) text[j] = 'a';
          else text[j] = 'b';
        test_leftmost_nonov_pair<text_offset_type>(text, text_length);
        test_leftmost_pair<text_offset_type>(text, text_length);
      }
    }

    // Clean up.
    delete[] text;

    // Print summary.
    fprintf(stderr, "All tests passed.\n");
  }

  {
    // Print initial message.
    fprintf(stderr, "Alphabet size = 3\n");

    // Allocate text.
#ifdef NDEBUG
    static const std::uint64_t max_text_length = 11;
#else
    static const std::uint64_t max_text_length = 7;
#endif
    char_type * const text = new char_type[max_text_length];

    // Run tests.
    for (std::uint64_t text_length = 1;
        text_length <= max_text_length; ++text_length) {

      // Print text length.
      fprintf(stderr, "Text length = %lu\n", text_length);

      // Run tests.
      std::uint64_t n_texts = 1;
      for (std::uint64_t j = 0; j < text_length; ++j)
        n_texts *= 3;

      for (std::uint64_t text_id = 0; text_id < n_texts; ++text_id) {
        std::uint64_t temp_text_id = text_id;
        for (std::uint64_t j = 0; j < text_length; ++j) {
          const std::uint64_t rest = temp_text_id % 3;

          if (rest == 0) text[j] = 'a';
          else if (rest == 1) text[j] = 'b';
          else text[j] = 'c';
          temp_text_id /= 3;
        }

        test_leftmost_nonov_pair<text_offset_type>(text, text_length);
        test_leftmost_pair<text_offset_type>(text, text_length);
      }
    }

    // Clean up.
    delete[] text;

    // Print summary.
    fprintf(stderr, "All tests passed.\n");
  }

  {
    // Print initial message.
    fprintf(stderr, "Alphabet size = 4\n");

    // Allocate text.
#ifdef NDEBUG
    static const std::uint64_t max_text_length = 9;
#else
    static const std::uint64_t max_text_length = 5;
#endif
    char_type * const text = new char_type[max_text_length];

    // Run tests.
    for (std::uint64_t text_length = 1;
        text_length <= max_text_length; ++text_length) {

      // Print text length.
      fprintf(stderr, "Text length = %lu\n", text_length);

      // Run tests.
      std::uint64_t n_texts = 1;
      for (std::uint64_t j = 0; j < text_length; ++j)
        n_texts *= 4;

      for (std::uint64_t text_id = 0; text_id < n_texts; ++text_id) {
        std::uint64_t temp_text_id = text_id;
        for (std::uint64_t j = 0; j < text_length; ++j) {
          const std::uint64_t rest = temp_text_id % 4;

          if (rest == 0) text[j] = 'a';
          else if (rest == 1) text[j] = 'b';
          else if (rest == 2) text[j] = 'c';
          else text[j] = 'd';
          temp_text_id /= 4;
        }

        test_leftmost_nonov_pair<text_offset_type>(text, text_length);
        test_leftmost_pair<text_offset_type>(text, text_length);
      }
    }

    // Clean up.
    delete[] text;

    // Print summary.
    fprintf(stderr, "All tests passed.\n");
  }
}
