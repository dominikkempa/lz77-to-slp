#ifndef __AVL_GRAMMAR_HPP_INCLUDED
#define __AVL_GRAMMAR_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "../utils/hash_table.hpp"
#include "../utils/karp_rabin_hashing.hpp"
#include "avl_grammar_node.hpp"
#include "avl_grammar_add_concat_nonterminal.hpp"
#include "avl_grammar_add_substring_nonterminal.hpp"


//=============================================================================
// A class storing AVL grammar.
//=============================================================================
template<typename char_type>
struct avl_grammar {
  typedef avl_grammar_node<char_type> node_type;

  // Class members.
  hash_table<std::uint64_t, const node_type*> m_hashes;
  std::vector<const node_type*> m_nonterminals;
  const node_type *m_root;
  const std::uint64_t m_hash_variable;
  const std::uint64_t m_mersenne_prime_exponent;

  // Constructor.
  avl_grammar(
      const std::uint64_t hash_variable,
      const std::uint64_t mersenne_prime_exponent) :
    m_root(NULL),
    m_hash_variable(hash_variable),
    m_mersenne_prime_exponent(mersenne_prime_exponent) {}

  // Print the string encoded by the grammar.
  void print_expansion() const {
    m_root->print_expansion();
  }

  // Return the number of nonterminals.
  std::uint64_t size() const {
    return m_nonterminals.size();
  }

  // Decode the text and write to a given array.
  void decode(
      char_type* &text,
      std::uint64_t &text_length) const {
    text_length = m_root->m_exp_len;
    text = new char_type[text_length];
    m_root->write_expansion(text);
  }

  // Test the AVL property of all nonterminals.
  bool test_avl_property() const {
    return m_root->test_avl_property();
  }

  // Collect Mersenne Karp-Rabin hashes in a vector.
  // Allows specifying variable and prime exponent.
  void collect_mersenne_karp_rabin_hashes(
      std::vector<std::uint64_t> &hashes) const {
    (void) m_root->collect_mersenne_karp_rabin_hashes(hashes);
  }

  // Collect Mersenne Karp-Rabin hashes in a hash table.
  // Allows specifying variable and prime exponent.
  void collect_mersenne_karp_rabin_hashes_2(
      hash_table<const node_type*, std::uint64_t> &hashes) const {
    (void) m_root->collect_mersenne_karp_rabin_hashes_2(hashes);
  }

  // Collect Mersenne Karp_Rabin hashes in a vector.
  // Relies on automatic choice of variable and exponent.
  void collect_mersenne_karp_rabin_hashes(
      std::vector<std::uint64_t> &hashes) const {
    collect_mersenne_karp_rabin_hashes(hashes);
  }

  // Collect Mersenne Karp_Rabin hashes in a hash table.
  // Relies on automatic choice of variable and exponent.
  void collect_mersenne_karp_rabin_hashes_2(
      hash_table<const node_type*, std::uint64_t> &hashes) const {
    collect_mersenne_karp_rabin_hashes_2(hashes);
  }

  // Count nodes in the pruned grammar.
  void count_nodes_in_pruned_grammar(
      hash_table<const node_type*, std::uint64_t> &hashes,
      hash_table<std::uint64_t, bool> &seen_hashes,
      std::uint64_t &current_count) const {
    m_root->count_nodes_in_pruned_grammar(hashes,
        seen_hashes, current_count);
  }

  // Collect pointers to all nonterminals reachable from the root.
  void collect_nonterminal_pointers(
      std::vector<const node_type*> &pointers) const {
    m_root->collect_nonterminal_pointers(pointers);
  }
};

#endif  // __AVL_GRAMMAR_HPP_INCLUDED
