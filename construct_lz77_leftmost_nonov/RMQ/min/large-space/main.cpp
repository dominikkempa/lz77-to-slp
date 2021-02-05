#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <ctime>
#include <unistd.h>

#include "utils.hpp"
#include "rmq_tree.hpp"


// Return the position of min in the range [beg..end).
// If no valid position exists, return `size'.
template<typename value_type>
std::uint64_t naive_rmq(
    const value_type *tab,
    std::uint64_t size,
    std::uint64_t beg,
    std::uint64_t end) {
  std::uint64_t ret_pos = size;
  std::uint64_t ret_val = std::numeric_limits<value_type>::max();

  beg = std::min(beg, size);
  end = std::min(end, size);

  if (beg >= end)
    return size;

  for (std::uint64_t i = beg; i < end; ++i) {
    if ((std::uint64_t)tab[i] < ret_val) {
      ret_val = tab[i];
      ret_pos = i;
    }
  }

  return ret_pos;
}

template<typename value_type>
bool naive_less(
    const value_type *tab,
    std::uint64_t size,
    std::uint64_t beg,
    std::uint64_t end,
    std::uint64_t threshold) {

  beg = std::min(beg, size);
  end = std::min(end, size);

  if (beg >= end)
    return false;

  for (std::uint64_t i = beg; i < end; ++i) 
    if ((std::uint64_t)tab[i] < threshold)
      return true;

  return false;
}

template<typename value_type>
void test_less(
    const value_type *tab,
    std::uint64_t size,
    std::uint64_t n_queries) {

  fprintf(stderr, "LESSTEST(sizeof(value_type) = %lu, size = %lu, "
      "n_queries = %lu)\n", (std::uint64_t)sizeof(value_type), size, n_queries);

  typedef rmq_tree<value_type> rmq_type;
  rmq_type *rmq = new rmq_type(tab, size);

  for (std::uint64_t query_id = 0; query_id < n_queries; ++query_id) {
    std::uint64_t beg =
      utils::random_int64((std::int64_t)0, (std::int64_t)(2 * size));
    std::uint64_t end =
      utils::random_int64((std::int64_t)0, (std::int64_t)(2 * size));
    std::uint64_t threshold = 0;

    if (utils::random_int64((std::int64_t)0, (std::int64_t)1) == 0)
      threshold = tab[utils::random_int64(
          (std::uint64_t)0, (std::uint64_t)(size - 1))];
    else {
      if (sizeof(value_type) < 8)
        threshold = utils::random_int64(
            (std::uint64_t)0,
            (std::uint64_t)std::numeric_limits<value_type>::max());
      else
        threshold = utils::random_int64(
            (std::uint64_t)0,
            (std::uint64_t)std::numeric_limits<std::int64_t>::max());
    }

    bool ret_correct = naive_less<value_type>(tab, size,
        beg, end, threshold);
    bool ret_computed = rmq->less(beg, end, threshold);

    if (ret_computed != ret_correct) {
      fprintf(stderr, "\nError:\n");
      fprintf(stderr, "  beg = %lu\n", beg);
      fprintf(stderr, "  end = %lu\n", end);
      fprintf(stderr, "  size = %lu\n", size);
      fprintf(stderr, "  threshold = %lu\n", threshold);
      fprintf(stderr, "  ret_correct = %lu\n", (std::uint64_t)ret_correct);
      fprintf(stderr, "  ret_computed = %lu\n", (std::uint64_t)ret_computed);
      if (size < 1000) {
        fprintf(stderr, "  tab: ");
        for (std::uint64_t i = 0; i < size; ++i)
          fprintf(stderr, "%lu ", (std::uint64_t)tab[i]);
        fprintf(stderr, "\n");
      }
      std::exit(EXIT_FAILURE);
    }
  }

  delete rmq;
}


template<typename value_type>
void test_rmq(
    const value_type *tab,
    std::uint64_t size,
    std::uint64_t n_queries) {

  fprintf(stderr, "RMQTEST(sizeof(value_type) = %lu, size = %lu, "
      "n_queries = %lu)\n", (std::uint64_t)sizeof(value_type), size, n_queries);

  typedef rmq_tree<value_type> rmq_type;
  rmq_type *rmq = new rmq_type(tab, size);

  for (std::uint64_t query_id = 0; query_id < n_queries; ++query_id) {
    std::uint64_t beg =
      utils::random_int64((std::int64_t)0, (std::int64_t)(2 * size));
    std::uint64_t end =
      utils::random_int64((std::int64_t)0, (std::int64_t)(2 * size));

    std::uint64_t ret_correct = naive_rmq<value_type>(tab, size, beg, end);
    std::uint64_t ret_computed = rmq->rmq(beg, end);

    if (ret_computed != ret_correct &&
        tab[ret_computed] != tab[ret_correct]) {

      fprintf(stderr, "\nError:\n");
      fprintf(stderr, "  beg = %lu\n", beg);
      fprintf(stderr, "  end = %lu\n", end);
      fprintf(stderr, "  size = %lu\n", size);
      fprintf(stderr, "  ret_correct = %lu\n", ret_correct);
      fprintf(stderr, "  ret_computed = %lu\n", ret_computed);
      if (size < 1000) {
        fprintf(stderr, "  tab: ");
        for (std::uint64_t i = 0; i < size; ++i)
          fprintf(stderr, "%lu ", (std::uint64_t)tab[i]);
        fprintf(stderr, "\n");
      }
      std::exit(EXIT_FAILURE);
    }
  }

  delete rmq;
}

template<typename value_type>
void test_rmq(
    std::uint64_t max_size,
    std::uint64_t n_queries) {

  std::uint64_t size = utils::random_int64((std::int64_t)1, (std::int64_t)max_size);
  value_type *tab = new value_type[size];
  if (sizeof(value_type) < 8) {
    for (std::uint64_t i = 0; i < size; ++i)
      tab[i] = utils::random_int64(
          (std::int64_t)0,
          (std::int64_t)std::numeric_limits<value_type>::max());
  } else {
    for (std::uint64_t i = 0; i < size; ++i)
      tab[i] = utils::random_int64(
          (std::int64_t)0,
          (std::int64_t)std::numeric_limits<std::int64_t>::max());
  }
  test_rmq(tab, size, n_queries);
  delete[] tab;
}

template<typename value_type>
void test_rmq() {
  for (std::uint64_t size = 1; size <= (1 << 19); size <<= 1)
    test_rmq<value_type>(size, 100000);
}

template<typename value_type>
void test_less(
    std::uint64_t max_size,
    std::uint64_t n_queries) {

  std::uint64_t size = utils::random_int64((std::int64_t)1, (std::int64_t)max_size);
  value_type *tab = new value_type[size];
  if (sizeof(value_type) < 8) {
    for (std::uint64_t i = 0; i < size; ++i)
      tab[i] = utils::random_int64(
          (std::int64_t)0,
          (std::int64_t)std::numeric_limits<value_type>::max());
  } else {
    for (std::uint64_t i = 0; i < size; ++i)
      tab[i] = utils::random_int64(
          (std::int64_t)0,
          (std::int64_t)std::numeric_limits<std::int64_t>::max());
  }
  test_less(tab, size, n_queries);
  delete[] tab;
}

template<typename value_type>
void test_less() {
  for (std::uint64_t size = 1; size <= (1 << 19); size <<= 1)
    test_less<value_type>(size, 1000000);
}

int main() {
  srand(time(0) + getpid());
  test_rmq<std::uint16_t>();
  test_rmq<std::uint32_t>();
  test_rmq<std::uint64_t>();
  test_less<std::uint16_t>();
  test_less<std::uint32_t>();
  test_less<std::uint64_t>();
  fprintf(stderr, "All tests passed.\n");
}

