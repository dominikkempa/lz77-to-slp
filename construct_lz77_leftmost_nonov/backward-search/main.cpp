#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <unistd.h>

#include "utils.hpp"
#include "compute_sa.hpp"
#include "byte_rank.hpp"


// Main testing function.
template<typename char_type,
  typename text_offset_type>
void test(
    char_type *text,
    std::uint64_t text_length,
    std::uint64_t alph_size,
    std::uint64_t n_queries) {

  // Compute suffix array of text.
  text_offset_type *sa = new text_offset_type[text_length];
  compute_sa(text, text_length, sa);

  // Comptute isa0.
  std::uint64_t isa0 = 0;
  for (std::uint64_t j = 0; j < text_length; ++j)
    if (sa[j] == 0)
      isa0 = j;

  // Compute BWT from SA.
  char_type *bwt = new char_type[text_length];
  bwt[0] = text[text_length - 1];
  {
    std::uint64_t pos = 1;
    for (std::uint64_t j = 0; j < text_length; ++j)
      if (sa[j] > 0) bwt[pos++] = text[(std::uint64_t)sa[j] - 1];
  }

  // Create rank over bwt.
  typedef byte_rank<> rank_type;
  rank_type *rank = new rank_type(bwt, text_length);

  // Compute counts array.
  std::vector<std::uint64_t> count(256, (std::uint64_t)0);
  for (std::uint64_t j = 0; j < text_length; ++j)
    ++count[(std::uint64_t)bwt[j]];
  {
    std::uint64_t sum = 0;
    for (std::uint64_t j = 0; j < 256; ++j) {
      std::uint64_t temp = count[j];
      count[j] = sum;
      sum += temp;
    }
  }

  // Allocate patters.
  static const std::uint64_t max_pat_length = 2 * text_length;
  char_type *pat = new char_type[max_pat_length];
  char_type *pat2 = new char_type[max_pat_length + 1];

  // Run tests.
  for (std::uint64_t query_id = 0; query_id < n_queries; ++query_id) {
    std::uint64_t pat_length = utils::random_int64((std::int64_t)0,
        (std::uint64_t)max_pat_length);

    // Choose random pattern.
    for (std::uint64_t j = 0; j < pat_length; ++j)
      pat[j] = 'a' + utils::random_int64(
          (std::int64_t)0,
          (std::int64_t)alph_size - 1);

    char_type c = 'a' + utils::random_int64(
        (std::int64_t)0,
        (std::int64_t)alph_size - 1);

    // Prepend pattern with c to create pat2.
    pat2[0] = c;
    for (std::uint64_t j = 0; j < pat_length; ++j)
      pat2[1 + j] = pat[j];
    std::uint64_t pat2_length = pat_length + 1;

    // Compute the number of suffixes of text that are smaller than pat.
    std::uint64_t beg = 0;
    for (std::uint64_t j = 0; j < text_length; ++j) {
      std::uint64_t lcp = 0;
      while (j + lcp < text_length && lcp < pat_length &&
          text[j + lcp] == pat[lcp]) ++lcp;
      if ((lcp < pat_length && j + lcp < text_length &&
            text[j + lcp] < pat[lcp]) ||
          (lcp < pat_length && j + lcp == text_length))
        ++beg;
    }

    // Compute the number of suffixes that are either
    // smaller than pat or have pat as a prefix.
    std::uint64_t end = 0;
    {
      std::uint64_t prefixed_count = 0;
      for (std::uint64_t j = 0; j < text_length; ++j) {
        std::uint64_t lcp = 0;
        while (j + lcp < text_length && lcp < pat_length &&
            text[j + lcp] == pat[lcp]) ++lcp;
        if (lcp == pat_length)
          ++prefixed_count;
      }
      end = beg + prefixed_count;
    }

    // Compute the number of suffixes of text that are smaller than pat2.
    std::uint64_t newbeg_correct = 0;
    for (std::uint64_t j = 0; j < text_length; ++j) {
      std::uint64_t lcp = 0;
      while (j + lcp < text_length && lcp < pat2_length &&
          text[j + lcp] == pat2[lcp]) ++lcp;
      if ((lcp < pat2_length && j + lcp < text_length &&
            text[j + lcp] < pat2[lcp]) ||
          (lcp < pat2_length && j + lcp == text_length))
        ++newbeg_correct;
    }

    // Compute the number of suffixes of text that
    // are smaller than pat2 or have pat2 as a prefix.
    std::uint64_t newend_correct = 0;
    {
      std::uint64_t prefixed_count = 0;
      for (std::uint64_t j = 0; j < text_length; ++j) {
        std::uint64_t lcp = 0;
        while (j + lcp < text_length && lcp < pat2_length &&
            text[j + lcp] == pat2[lcp]) ++lcp;
        if (lcp == pat2_length)
          ++prefixed_count;
      }
      newend_correct = newbeg_correct + prefixed_count;
    }

    // Compute answer using the algorithm.
    std::uint64_t newbeg_computed = count[c] +
      rank->query(beg + (pat_length > 0) - (beg > isa0), c);
    std::uint64_t newend_computed = count[c] +
      rank->query(end + 1 - (end > isa0), c);

    // Print error message if necessary.
    if (newbeg_correct != newbeg_computed ||
        newend_correct != newend_computed) {
      fprintf(stderr, "\nError:\n");
      if (text_length < 100) {
        fprintf(stderr, "  text = ");
        for (std::uint64_t j = 0; j < text_length; ++j)
          fprintf(stderr, "%c", text[j]);
        fprintf(stderr, "\n");
      }
      fprintf(stderr, "  text_length = %lu\n", text_length);
      fprintf(stderr, "  pat_length = %lu\n", pat_length);
      fprintf(stderr, "  pat = ");
      for (std::uint64_t j = 0; j < pat_length; ++j)
        fprintf(stderr, "%c", pat[j]);
      fprintf(stderr, "\n");
      fprintf(stderr, "  c = %c\n", c);
      fprintf(stderr, "  isa0 = %lu\n", isa0);
      fprintf(stderr, "  bwt = ");
      for (std::uint64_t j = 0; j < text_length; ++j)
        fprintf(stderr, "%c", bwt[j]);
      fprintf(stderr, "\n");
      fprintf(stderr, "  beg = %lu\n", beg);
      fprintf(stderr, "  newbeg_correct = %lu\n", newbeg_correct);
      fprintf(stderr, "  newbeg_computed = %lu\n", newbeg_computed);
      fprintf(stderr, "  end = %lu\n", end);
      fprintf(stderr, "  newend_correct = %lu\n", newend_correct);
      fprintf(stderr, "  newend_computed = %lu\n", newend_computed);
      std::exit(EXIT_FAILURE);
    }
  }

  // Clean up.
  delete rank;
  delete[] pat;
  delete[] pat2;
  delete[] bwt;
  delete[] sa;
}

int main() {

  // Init rand.
  srand(time(0) + getpid());

  // Declare types.
  typedef std::uint8_t char_type;
  typedef std::uint32_t text_offset_type;

  // Set number of tests.
  static const std::uint64_t n_tests = 100;

  {

    // Print initial message.
    fprintf(stderr, "Alphabet size = 2\n");

    // Allocate text.
    static const std::uint64_t max_text_length = 20;
    char_type *text = new char_type[max_text_length];

    // Run tests.
    for (std::uint64_t text_length = 1;
        text_length <= max_text_length; ++text_length) {

      // Print text length.
      fprintf(stderr, "Text length = %lu\n", text_length);

      // Run tests.
      std::uint64_t n_texts = ((std::uint64_t)1 << text_length);
      for (std::uint64_t text_id = 0; text_id < n_texts; ++text_id) {
        for (std::uint64_t j = 0; j < text_length; ++j)
          if (text_id & ((std::uint64_t)1 << j)) text[j] = 'a';
          else text[j] = 'b';
        test<char_type, text_offset_type>(text, text_length,
            (std::uint64_t)2, n_tests);
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
    static const std::uint64_t max_text_length = 13;
    char_type *text = new char_type[max_text_length];

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
          std::uint64_t rest = temp_text_id % 3;

          if (rest == 0) text[j] = 'a';
          else if (rest == 1) text[j] = 'b';
          else text[j] = 'c';
          temp_text_id /= 3;
        }

        test<char_type, text_offset_type>(text, text_length,
            (std::uint64_t)3, n_tests);
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
    static const std::uint64_t max_text_length = 10;
    char_type *text = new char_type[max_text_length];

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
          std::uint64_t rest = temp_text_id % 4;

          if (rest == 0) text[j] = 'a';
          else if (rest == 1) text[j] = 'b';
          else if (rest == 2) text[j] = 'c';
          else text[j] = 'd';
          temp_text_id /= 4;
        }

        test<char_type, text_offset_type>(text, text_length,
            (std::uint64_t)4, n_tests);
      }
    }

    // Clean up.
    delete[] text;

    // Print summary.
    fprintf(stderr, "All tests passed.\n");
  }
}
