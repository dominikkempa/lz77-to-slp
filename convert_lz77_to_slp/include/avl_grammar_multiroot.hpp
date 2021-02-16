#ifndef __AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
#define __AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <map>
#include <algorithm>

#include "hash_table.hpp"
#include "karp_rabin_hashing.hpp"
#include "avl_grammar_node.hpp"
#include "avl_grammar_add_concat_nonterminal.hpp"
#include "avl_grammar_add_substring_nonterminal.hpp"


//=============================================================================
// A class storing multiroot AVL grammar.
//=============================================================================
template<typename char_type>
struct avl_grammar_multiroot {
  typedef avl_grammar_node<char_type> node_type;
  typedef std::uint64_t key_type;
  typedef const node_type* value_type;
  typedef typename std::map<key_type, value_type> map_type;
  typedef typename map_type::const_iterator const_iter_type;
  typedef typename map_type::iterator iter_type;

  // Class members.
  std::vector<const node_type*> m_nonterminals;
  map_type m_roots;

  // Constructor.
  avl_grammar_multiroot() {}

  // Print the string encoded by the grammar.
  void print_expansion() const {
    for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
      it->second->print_expansion();
  }

  // Return the number of nonterminals.
  std::uint64_t size() const {
    return m_nonterminals.size();
  }

  // Decode the text and write to a given array.
  void decode(
      char_type* &text,
      std::uint64_t &text_length) const {
    text_length = 0;
    for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
      text_length += it->second->m_exp_len;
    text = new char_type[text_length];
    std::uint64_t ptr = 0;
    for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it) {
      it->second->write_expansion(text + ptr);
      ptr += it->second->m_exp_len;
    }
  }

  // Test the AVL property of all nonterminals.
  bool test_avl_property() const {
    for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
      if (it->second->test_avl_property() == false)
        return false;
    return true;
  }

  // Collect Karp-Rabin hashes in a vector.
  void collect_karp_rabin_hashes(
      std::vector<std::uint64_t> &hashes,
      const std::uint64_t a = (std::uint64_t)999285268,
      const std::uint64_t p = (std::uint64_t)1000000009) const {
    for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
      (void) it->second->collect_karp_rabin_hashes(hashes, a, p);
  }

  // Collect Mersenne Karp-Rabin hashes in a vector.
  // Allows specifying variable and prime exponent.
  void collect_mersenne_karp_rabin_hashes(
      std::vector<std::uint64_t> &hashes,
      const std::uint64_t hash_variable,
      const std::uint64_t mersenne_prime_exponent) const {
    for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
      (void) it->second->collect_mersenne_karp_rabin_hashes(
          hashes, hash_variable, mersenne_prime_exponent);
  }

  // Collect Mersenne Karp-Rabin hashes in a hash table.
  // Allows specifying variable and prime exponent.
  void collect_mersenne_karp_rabin_hashes_2(
      hash_table<const node_type*, std::uint64_t> &hashes,
      const std::uint64_t hash_variable,
      const std::uint64_t mersenne_prime_exponent) const {
    for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
      (void) it->second->collect_mersenne_karp_rabin_hashes_2(
          hashes, hash_variable, mersenne_prime_exponent);
  }

  // Collect Mersenne Karp_Rabin hashes in a vector.
  // Relies on automatic choice of variable and exponent.
  void collect_mersenne_karp_rabin_hashes(
      std::vector<std::uint64_t> &hashes) const {
    const std::uint64_t mersenne_prime_exponent = 61;
    const std::uint64_t hash_variable = 
      rand_mod_mersenne(mersenne_prime_exponent);
    collect_mersenne_karp_rabin_hashes(hashes,
        hash_variable, mersenne_prime_exponent);
  }

  // Collect Mersenne Karp_Rabin hashes in a hash table.
  // Relies on automatic choice of variable and exponent.
  void collect_mersenne_karp_rabin_hashes_2(
      hash_table<const node_type*, std::uint64_t> &hashes) const {
    const std::uint64_t mersenne_prime_exponent = 61;
    const std::uint64_t hash_variable = 
      rand_mod_mersenne(mersenne_prime_exponent);
    collect_mersenne_karp_rabin_hashes_2(hashes,
        hash_variable, mersenne_prime_exponent);
  }

  // Count nodes in the pruned grammar.
  void count_nodes_in_pruned_grammar(
      hash_table<const node_type*, std::uint64_t> &hashes,
      hash_table<std::uint64_t, bool> &seen_hashes,
      std::uint64_t &current_count) const {
    for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
      it->second->count_nodes_in_pruned_grammar(
          hashes, seen_hashes, current_count);
  }

  // Collect pointers to all nonterminals reachable from the root.
  void collect_nonterminal_pointers(
      std::vector<const node_type*> &pointers) const {
    for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
      it->second->collect_nonterminal_pointers(pointers);
  }

  // Merge nonterminals enclosed in the interval [begin..end).
  void merge_enclosed_nonterminals(
      const std::uint64_t begin,
      const std::uint64_t end) {

    // Check arguments.
    if (begin > end) {
      fprintf(stderr, "n\nError: merge_enclosed_nonterminals: "
          "begin > end!\n");
      std::exit(EXIT_FAILURE);
    }

    // Handle special case.
    if (m_roots.empty())
      return;

    // Find the range [it_begin..it_end)
    // of nonterminals to merge.
    iter_type it_begin = m_roots.end();
    std::uint64_t prev_nonterminals_total_length = 0;
    if (begin == 0) {
      it_begin = m_roots.begin();
      prev_nonterminals_total_length = 0;
    } else {
      it_begin = m_roots.lower_bound(begin);
      prev_nonterminals_total_length = it_begin->first;
      ++it_begin;
    }
    iter_type it_end = it_begin;
    while (it_end != m_roots.end() &&
        it_end->first <= end)
      ++it_end;

    // Merge nonterminals in [it_begin..it_end).
    // TODO: reduce complexity from O(n^2) to O(n log n).
    if (it_begin != it_end) {

      // Collect the sequence of pointers into a vector.
      std::vector<value_type> v;
      for (iter_type it = it_begin; it != it_end; ++it)
        v.push_back(it->second);

      // We can now erase the roots in the grammar.
      m_roots.erase(it_begin, it_end);

      // Merge neighbors, always starting with the shortest one
      // as long as there is >= 2 nonterminals left.
      while (v.size() > 1) {

        // Find the nonterminal with the smallest height.
        std::uint64_t smallest_height_id = 0;
        for (std::uint64_t i = 1; i < v.size(); ++i) {
          if (v[i]->m_height < v[smallest_height_id]->m_height)
            smallest_height_id = i;
        }

        // Merge the nonterminal with the smaller height with
        // one of its beighbors (whichever is shorter).
        if (smallest_height_id == 0 ||
            (smallest_height_id + 1 < v.size() &&
             v[smallest_height_id + 1]->m_height <=
             v[smallest_height_id - 1]->m_height)) {

          // Only right neighbor exists, or both exist
          // and the right one is not taller than the left
          // one. End result: merge with the right neighbor.
          value_type ret =
            add_concat_nonterminal<char_type>(m_nonterminals,
                v[smallest_height_id], v[smallest_height_id + 1]);
          v.erase(v.begin() + smallest_height_id);
          v[smallest_height_id] = ret;
        } else {

          // Only left neighbor exists, or both exists
          // and the left one is not taller than the
          // right one. End result: merge with left neighbor.
          value_type ret =
            add_concat_nonterminal<char_type>(m_nonterminals,
                v[smallest_height_id - 1], v[smallest_height_id]);
          v.erase(v.begin() + (smallest_height_id - 1));
          v[smallest_height_id - 1] = ret;
        }
      }

      // Add the resulting nonterminal into the list of roots.
      value_type val = v[0];
      key_type key = prev_nonterminals_total_length + val->m_exp_len;
      m_roots[key] = val;
    }
  }

  // Add the nonterminal expanding to the string T[begin..end).
  // There are many ways of doing this. For now we employ a simple
  // strategy, but in the future, this function will most likely
  // return the sequence of nonterminals, and not just one.
  const node_type *add_substring_nonterminal(
      const std::uint64_t begin,
      const std::uint64_t end) {

    // Check arguments.
    if (begin > end) {
      fprintf(stderr, "n\nError: merge_enclosed_nonterminals: "
          "begin > end!\n");
      std::exit(EXIT_FAILURE);
    }

    // First, ensure that [begin..end) overlaps
    // at most three roots of the grammar.
    merge_enclosed_nonterminals(begin, end);

    // Find the leftmost root whose expansion
    // overlaps or touches block T[begin..end).
    iter_type it = m_roots.lower_bound(begin);
    --it; // delete me.

    // TODO
    return NULL;
  }
};

#endif  // __AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
