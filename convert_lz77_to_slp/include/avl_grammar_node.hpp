#ifndef __AVL_GRAMMAR_NODE_HPP_INCLUDED
#define __AVL_GRAMMAR_NODE_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <map>
#include <algorithm>

#include "hash_table.hpp"
#include "karp_rabin_hashing.hpp"


//=============================================================================
// A class used to represent the nonterminal in the AVL grammar.
//=============================================================================
template<typename char_type>
struct avl_grammar_node {
  typedef avl_grammar_node<char_type> node_type;

  // Class members.
  const char_type m_char;
  const std::uint8_t m_height;
  const std::uint64_t m_exp_len;
  const node_type * const m_left;
  const node_type * const m_right;

  // Default constructor.
  avl_grammar_node() :
      m_char((char_type)0),
      m_height(0),
      m_exp_len(1),
      m_left(NULL),
      m_right(NULL) {}

  // Constructor for a node expanding to a single symbol.
  avl_grammar_node(const char_type c) :
      m_char(c),
      m_height(0),
      m_exp_len(1),
      m_left(NULL),
      m_right(NULL) {}

  // Constructor for a standard nonterminal
  // (expanding to two other nonterminals).
  avl_grammar_node(
      const node_type * const left,
      const node_type * const right) :
        m_char((char_type)0),
        m_height(std::max(left->m_height, right->m_height) + 1),
        m_exp_len(left->m_exp_len + right->m_exp_len),
        m_left(left),
        m_right(right) {}

  // Print the string encoded by the grammar.
  void print_expansion() const {
    if (m_height == 0)
      fprintf(stderr, "%c", (char)m_char);
    else {
      m_left->print_expansion();
      m_right->print_expansion();
    }
  }

  // Write the expansion into the given array.
  void write_expansion(char_type * const text) const {
    if (m_height == 0)
      text[0] = m_char;
    else {
      m_left->write_expansion(text);
      m_right->write_expansion(text + m_left->m_exp_len);
    }
  }

  // Test the AVL propert of a subtree.
  bool test_avl_property() const {
    if (m_height == 0)
      return true;
    else {
      if ((!(m_left->test_avl_property())) ||
          (!(m_right->test_avl_property())) ||
          (m_right->m_height > m_left->m_height &&
           m_right->m_height - m_left->m_height > 1) ||
          (m_left->m_height > m_right->m_height &&
           m_left->m_height - m_right->m_height > 1))
        return false;
      else return true;
    }
  }

  // Collect Karp-Rabin hashes of all nonterminals.
  std::uint64_t collect_karp_rabin_hashes(
      std::vector<std::uint64_t> &hashes,
      const std::uint64_t a,
      const std::uint64_t p) const {
    if (m_height == 0) {
      const std::uint64_t h = m_char % p;
      hashes.push_back(h);
      return h;
    } else {
      const std::uint64_t left_hash =
        m_left->collect_karp_rabin_hashes(hashes, a, p);
      const std::uint64_t right_hash =
        m_right->collect_karp_rabin_hashes(hashes, a, p);
      std::uint64_t h =
        (left_hash * mod_pow<std::uint64_t>(a, m_right->m_exp_len, p)) % p;
      h = (h + right_hash) % p;
      hashes.push_back(h);
      return h;
    }
  }

  // Collect Mersenne Karp-Rabin hashes of all nonterminals.
  std::uint64_t collect_mersenne_karp_rabin_hashes(
      std::vector<std::uint64_t> &hashes,
      const std::uint64_t hash_variable,
      const std::uint64_t mersenne_prime_exponent) const {
    if (m_height == 0) {
      const std::uint64_t h = mod_mersenne(m_char, mersenne_prime_exponent);
      hashes.push_back(h);
      return h;
    } else {
      const std::uint64_t left_hash =
        m_left->collect_mersenne_karp_rabin_hashes(hashes,
            hash_variable, mersenne_prime_exponent);
      const std::uint64_t right_hash =
        m_right->collect_mersenne_karp_rabin_hashes(hashes,
            hash_variable, mersenne_prime_exponent);
      std::uint64_t h = 0;
      {
        const std::uint64_t pow = pow_mod_mersenne(hash_variable,
            m_right->m_exp_len, mersenne_prime_exponent);
        const std::uint64_t tmp = mul_mod_meresenne(left_hash,
            pow, mersenne_prime_exponent);
        h = mod_mersenne(tmp + right_hash, mersenne_prime_exponent);
      }
      hashes.push_back(h);
      return h;
    }
  }

  // Collect reachable nonterminals.
  void collect_nonterminal_pointers(
      std::vector<const node_type*> &pointers) const {
    pointers.push_back(this);
    if (m_height > 0) {
      m_left->collect_nonterminal_pointers(pointers);
      m_right->collect_nonterminal_pointers(pointers);
    }
  }

  // Collect Mersenne Karp-Rabin hashes of all nonterminals.
  std::uint64_t collect_mersenne_karp_rabin_hashes_2(
      hash_table<const node_type*, std::uint64_t> &hashes,
      const std::uint64_t hash_variable,
      const std::uint64_t mersenne_prime_exponent) const {
    if (m_height == 0) {
      const std::uint64_t h =
        mod_mersenne(m_char, mersenne_prime_exponent);
      hashes.insert(this, h);
      return h;
    } else {
      const std::uint64_t left_hash =
        m_left->collect_mersenne_karp_rabin_hashes_2(hashes,
            hash_variable, mersenne_prime_exponent);
      const std::uint64_t right_hash =
        m_right->collect_mersenne_karp_rabin_hashes_2(hashes,
            hash_variable, mersenne_prime_exponent);
      std::uint64_t h = 0;
      {
        const std::uint64_t pow = pow_mod_mersenne(hash_variable,
            m_right->m_exp_len, mersenne_prime_exponent);
        const std::uint64_t tmp = mul_mod_meresenne(left_hash,
            pow, mersenne_prime_exponent);
        h = mod_mersenne(tmp + right_hash, mersenne_prime_exponent);
      }
      hashes.insert(this, h);
      return h;
    }
  }

  // Compute the number of nodes in the pruned grammar.
  void count_nodes_in_pruned_grammar(
      hash_table<const node_type*, std::uint64_t> &hashes,
      hash_table<std::uint64_t, bool> &seen_hashes,
      std::uint64_t &current_count) const {
    if (m_height == 0) {
      const std::uint64_t *h = hashes.find(this);
      if (seen_hashes.find(*h) == NULL) {
        seen_hashes.insert(*h, true);
        ++current_count;
      }
    } else {
      const std::uint64_t *h = hashes.find(this);
      if (seen_hashes.find(*h) == NULL) {
        seen_hashes.insert(*h, true);
        ++current_count;
        m_left->count_nodes_in_pruned_grammar(hashes,
            seen_hashes, current_count);
        m_right->count_nodes_in_pruned_grammar(hashes,
            seen_hashes, current_count);
      }
    }
  }
};

//=============================================================================
// Hash functions of the appropriate type.
// Used in the hash table used to prune the grammar.
//=============================================================================
typedef const avl_grammar_node<std::uint8_t>* get_hash_ptr_type;

template<>
std::uint64_t get_hash(const get_hash_ptr_type &x) {
  return (std::uint64_t)x * (std::uint64_t)29996224275833;
}

template<>
std::uint64_t get_hash(const std::uint64_t &x) {
  return (std::uint64_t)x * (std::uint64_t)4972548694736365;
}

#endif  // __AVL_GRAMMAR_NODE_HPP_INCLUDED