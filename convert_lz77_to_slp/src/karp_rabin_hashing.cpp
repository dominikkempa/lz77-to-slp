#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "../include/utils.hpp"
#include "../include/karp_rabin_hashing.hpp"


// Return (a * b) mod p, where p = (2^k) - 1.
// Requires a, b <= 2^k. Tested for k = 1, .., 63.
std::uint64_t mul_mod_meresenne(std::uint64_t a,
    std::uint64_t b, std::uint64_t k) {
  std::uint64_t p = ((std::uint64_t)1 << k) - 1;
  __extension__ unsigned __int128 ab =
    (unsigned __int128)a *
    (unsigned __int128)b;
  std::uint64_t lo = (std::uint64_t)ab;
  std::uint64_t hi = (ab >> 64);
  lo = (lo & p) + ((lo >> k) + (hi << (64 - k)));
  lo = (lo & p) + (lo >> k);
  return lo == p ? 0 : lo;
}

// Return a mod p, where p = (2^k) - 1.
// Works for any a in [0..2^64).
// Tested for k = 1, .., 63.
std::uint64_t mod_mersenne(std::uint64_t a, std::uint64_t k) {
  std::uint64_t p = ((std::uint64_t)1 << k) - 1;
  if (k < 32) {

    // We need to check if a <= 2^(2k).
    std::uint64_t threshold = ((std::uint64_t)1 << (k << 1));
    if (a <= threshold) {
      a = (a & p) + (a >> k);
      a = (a & p) + (a >> k);
      return a == p ? 0 : a;
    } else return a % p;
  } else {

    // We are guaranteed that a < 2^(2k)
    // because a < 2^64 <= 2^(2k).
    a = (a & p) + (a >> k);
    a = (a & p) + (a >> k);
    return a == p ? 0 : a;
  }
}

// Return random number x in [0..p), where p = (2^k) - 1.
inline std::uint64_t rand_mod_mersenne(std::uint64_t k) {
  std::uint64_t p = ((std::uint64_t)1 << k) - 1;
  return utils::random_int64((std::uint64_t)0, p - 1);
}
