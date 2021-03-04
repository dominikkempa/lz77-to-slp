#ifndef __AVL_GRAMMAR_NODE_HPP_INCLUDED
#define __AVL_GRAMMAR_NODE_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <map>
#include <algorithm>

#include "../utils/hash_table.hpp"
#include "../utils/karp_rabin_hashing.hpp"


//=============================================================================
// A class used to represent the nonterminal in the AVL grammar.
//=============================================================================
template<typename char_type>
struct avl_grammar_node {
  private:

    // Declare typedefs.
    typedef avl_grammar_node<char_type> node_type;

  public:

    // Class members.
    const char_type m_char;
    const std::uint8_t m_height;
    const std::uint64_t m_exp_len;
    const std::uint64_t m_kr_hash;
    const node_type * const m_left;
    const node_type * const m_right;

    // Default constructor.
    avl_grammar_node() :
      m_char((char_type)0),
      m_height(0),
      m_exp_len(1),
      m_kr_hash(0),
      m_left(NULL),
      m_right(NULL) {}

    // Constructor for a node expanding to a single symbol.
    avl_grammar_node(const char_type c) :
      m_char(c),
      m_height(0),
      m_exp_len(1),
      m_kr_hash(karp_rabin_hashing::hash_char(c)),
      m_left(NULL),
      m_right(NULL) {}

    // Constructor for non-single-symbol nonterminal.
    avl_grammar_node(
        const node_type * const left,
        const node_type * const right) :
          m_char((char_type)0),
          m_height(std::max(left->m_height, right->m_height) + 1),
          m_exp_len(left->m_exp_len + right->m_exp_len),
          m_kr_hash(
              karp_rabin_hashing::concat(
                left->m_kr_hash,
                right->m_kr_hash,
                right->m_exp_len)),
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
      if (m_height > 0 &&
          ((!(m_left->test_avl_property())) ||
           (!(m_right->test_avl_property())) ||
           (m_right->m_height > m_left->m_height &&
            m_right->m_height - m_left->m_height > 1) ||
           (m_left->m_height > m_right->m_height &&
            m_left->m_height - m_right->m_height > 1)))
        return false;
      return true;
    }

    // Collect Mersenne Karp-Rabin hashes of all nonterminals.
    std::uint64_t collect_mersenne_karp_rabin_hashes(
        std::vector<std::uint64_t> &hashes) const {
      if (m_height == 0) {
        const std::uint64_t h = karp_rabin_hashing::hash_char(m_char);
        hashes.push_back(h);
        return h;
      } else {
        const std::uint64_t left_hash =
          m_left->collect_mersenne_karp_rabin_hashes(hashes);
        const std::uint64_t right_hash =
          m_right->collect_mersenne_karp_rabin_hashes(hashes);
        const std::uint64_t right_len = m_right->m_exp_len;
        const std::uint64_t h = karp_rabin_hashing::concat(
            left_hash, right_hash, right_len);
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
        hash_table<const node_type*, std::uint64_t> &hashes) const {
      if (m_height == 0) {
        const std::uint64_t h = karp_rabin_hashing::hash_char(m_char);
        hashes.insert(this, h);
        return h;
      } else {
        const std::uint64_t left_hash =
          m_left->collect_mersenne_karp_rabin_hashes_2(hashes);
        const std::uint64_t right_hash =
          m_right->collect_mersenne_karp_rabin_hashes_2(hashes);
        const std::uint64_t right_len = m_right->m_exp_len;
        const std::uint64_t h = karp_rabin_hashing::concat(
            left_hash, right_hash, right_len);
        hashes.insert(this, h);
        return h;
      }
    }

    // Compute the number of nodes in the pruned grammar.
    void count_nodes_in_pruned_grammar(
        hash_table<const node_type*, std::uint64_t> &hashes,
        hash_table<std::uint64_t, bool> &seen_hashes,
        std::uint64_t &current_count) const {
      const std::uint64_t * const h = hashes.find(this);
      if (seen_hashes.find(*h) == NULL) {
        seen_hashes.insert(*h, true);
        ++current_count;
        if (m_height != 0) {
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

template<typename char_type>
std::uint64_t merge_hashes(
    const avl_grammar_node<char_type> * const left,
    const avl_grammar_node<char_type> * const right) {
  const std::uint64_t left_hash = left->m_kr_hash;
  const std::uint64_t right_hash = right->m_kr_hash;
  const std::uint64_t right_len = right->m_exp_len;
  const std::uint64_t h = karp_rabin_hashing::concat(
      left_hash, right_hash, right_len);
  return h;
}

template<typename char_type>
std::uint64_t append_hash(
    const std::uint64_t left_hash,
    const avl_grammar_node<char_type> * const right) {
  const std::uint64_t right_hash = right->m_kr_hash;
  const std::uint64_t right_len = right->m_exp_len;
  const std::uint64_t h = karp_rabin_hashing::concat(
      left_hash, right_hash, right_len);
  return h;
}

#endif  // __AVL_GRAMMAR_NODE_HPP_INCLUDED
