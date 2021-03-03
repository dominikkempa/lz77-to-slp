#ifndef __KARP_RABIN_HASHING_HPP_INCLUDED
#define __KARP_RABIN_HASHING_HPP_INCLUDED


//=============================================================================
// Compute a^n (mod p) in O(log n) time.
//=============================================================================
template<typename int_type = std::uint64_t>
int_type mod_pow(
    const int_type a,
    int_type n,
    const int_type p) {
  int_type pow = a % p;
  int_type ret = 1 % p;
  while (n > 0) {
    if (n & 1)
      ret = (ret * pow) % p;
    pow = (pow * pow) % p;
    n >>= 1;
  }
  return ret;
}

std::uint64_t mul_mod_meresenne(const std::uint64_t a,
    const std::uint64_t b, const std::uint64_t k);
std::uint64_t mod_mersenne(std::uint64_t a, const std::uint64_t k);
std::uint64_t rand_mod_mersenne(const std::uint64_t k);
std::uint64_t pow_mod_mersenne(const std::uint64_t a,
    std::uint64_t n, const std::uint64_t k);

#endif  // __KARP_RABIN_HASHING_HPP_INCLUDED
