#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <ctime>
#include <unistd.h>

#include "utils.hpp"
#include "byte_rank.hpp"


std::uint64_t naive_rank(
    const std::uint8_t *text,
    std::uint64_t length,
    std::uint64_t i,
    std::uint8_t c) {

  i = std::min(i, length);
  std::uint64_t ret = 0;
  for (std::uint64_t j = 0; j < i; ++j)
    if (text[j] == c) ++ret;
  return ret;
}

template<std::uint64_t block_size_log,
  std::uint64_t sblock_size_log>
void test(
    const std::uint8_t *text,
    std::uint64_t length,
    std::uint64_t n_queries) {

  fprintf(stderr, "TEST(block_size_log = %lu, sblock_size_log = %lu, "
      "length = %lu, n_queries = %lu)\n", block_size_log, sblock_size_log,
      length, n_queries);

  typedef byte_rank<block_size_log, sblock_size_log> rank_type;
  rank_type *rank = new rank_type(text, length);

  for (std::uint64_t query_id = 0; query_id < n_queries; ++query_id) {
    std::uint64_t i = utils::random_int64(
        (std::uint64_t)0,
        (std::uint64_t)(2 * length));
    std::uint8_t c = utils::random_int64(
        (std::uint64_t)0,
        (std::uint64_t)255);

    std::uint64_t ret_correct = naive_rank(text, length, i, c);
    std::uint64_t ret_computed = rank->query(i, c);

    if (ret_computed != ret_correct) {
      fprintf(stderr, "\nError:\n");
      fprintf(stderr, "  i = %lu\n", i);
      fprintf(stderr, "  c = %lu\n", (std::uint64_t)c);
      fprintf(stderr, "  ret correct = %lu\n", ret_correct);
      fprintf(stderr, "  ret computed = %lu\n", ret_computed);
      std::exit(EXIT_FAILURE);
    }
  }

  delete rank;
}

void test(
    const std::uint8_t *text,
    std::uint64_t length,
    std::uint64_t n_queries) {
  test<0, 0>(text, length, n_queries);
  test<0, 1>(text, length, n_queries);
  test<0, 2>(text, length, n_queries);
  test<0, 3>(text, length, n_queries);
  test<0, 4>(text, length, n_queries);
  test<0, 5>(text, length, n_queries);
  test<0, 6>(text, length, n_queries);
  test<0, 7>(text, length, n_queries);
  test<0, 8>(text, length, n_queries);
  test<0, 9>(text, length, n_queries);
  test<1, 1>(text, length, n_queries);
  test<1, 2>(text, length, n_queries);
  test<1, 3>(text, length, n_queries);
  test<1, 4>(text, length, n_queries);
  test<1, 5>(text, length, n_queries);
  test<1, 6>(text, length, n_queries);
  test<1, 7>(text, length, n_queries);
  test<1, 8>(text, length, n_queries);
  test<1, 9>(text, length, n_queries);
  test<2, 2>(text, length, n_queries);
  test<2, 3>(text, length, n_queries);
  test<2, 4>(text, length, n_queries);
  test<2, 5>(text, length, n_queries);
  test<2, 6>(text, length, n_queries);
  test<2, 7>(text, length, n_queries);
  test<2, 8>(text, length, n_queries);
  test<2, 9>(text, length, n_queries);
  test<3, 3>(text, length, n_queries);
  test<3, 4>(text, length, n_queries);
  test<3, 5>(text, length, n_queries);
  test<3, 6>(text, length, n_queries);
  test<3, 7>(text, length, n_queries);
  test<3, 8>(text, length, n_queries);
  test<3, 9>(text, length, n_queries);
  test<4, 4>(text, length, n_queries);
  test<4, 5>(text, length, n_queries);
  test<4, 6>(text, length, n_queries);
  test<4, 7>(text, length, n_queries);
  test<4, 8>(text, length, n_queries);
  test<4, 9>(text, length, n_queries);
  test<5, 5>(text, length, n_queries);
  test<5, 6>(text, length, n_queries);
  test<5, 7>(text, length, n_queries);
  test<5, 8>(text, length, n_queries);
  test<5, 9>(text, length, n_queries);
  test<6, 6>(text, length, n_queries);
  test<6, 7>(text, length, n_queries);
  test<6, 8>(text, length, n_queries);
  test<6, 9>(text, length, n_queries);
  test<7, 7>(text, length, n_queries);
  test<7, 8>(text, length, n_queries);
  test<7, 9>(text, length, n_queries);
  test<8, 8>(text, length, n_queries);
  test<8, 9>(text, length, n_queries);
  test<9, 9>(text, length, n_queries);
}

void test(std::uint64_t max_length) {
  std::uint64_t length = utils::random_int64(
      (std::uint64_t)1,
      (std::uint64_t)max_length);
  std::uint8_t *text = new std::uint8_t[length];
  for (std::uint64_t i = 0; i < length; ++i)
    text[i] = utils::random_int64(
        (std::uint64_t)0,
        (std::uint64_t)255);
  test(text, length, 10000);
  delete[] text;
}

int main() {
  srand(time(0) + getpid());

  for (std::uint64_t max_length = 1;
      max_length <= (1 << 19); max_length <<= 1)
    test(max_length);

  fprintf(stderr, "All tests passed.\n");
}

