#ifndef __KARP_RABIN_HASHING_HPP_INCLUDED
#define __KARP_RABIN_HASHING_HPP_INCLUDED

#include <cstdint>


namespace karp_rabin_hashing {

extern std::uint64_t hash_variable;
extern std::uint64_t mersenne_prime_exponent;

std::uint64_t mul_mod_mersenne(const std::uint64_t a,
    const std::uint64_t b, const std::uint64_t k);
std::uint64_t mod_mersenne(std::uint64_t a, const std::uint64_t k);
std::uint64_t rand_mod_mersenne(const std::uint64_t k);
std::uint64_t pow_mod_mersenne(const std::uint64_t a,
    std::uint64_t n, const std::uint64_t k);
std::uint64_t concat(const std::uint64_t lh,
    const std::uint64_t rh, const std::uint64_t rl);

template<typename char_type>
std::uint64_t hash_char(char_type c) {
  return mod_mersenne((std::uint64_t)c, mersenne_prime_exponent);
}

void init();

}  // namespace karp_rabin_hashing

#endif  // __KARP_RABIN_HASHING_HPP_INCLUDED
