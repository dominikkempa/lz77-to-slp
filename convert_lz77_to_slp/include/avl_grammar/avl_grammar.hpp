#ifndef __AVL_GRAMMAR_HPP_INCLUDED
#define __AVL_GRAMMAR_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "../utils/hash_table.hpp"
#include "avl_grammar_node.hpp"
#include "avl_grammar_add_concat_nonterminal.hpp"


//=============================================================================
// A class storing AVL grammar.
//=============================================================================
template<typename char_type>
struct avl_grammar {
  private:

    // Declare typedefs.
    typedef avl_grammar_node<char_type> node_type;

  public:

    // Class members.
    std::vector<const node_type*> m_nonterminals;
    const node_type *m_root;

  public:

    // Constructor.
    avl_grammar() :
      m_root(NULL) {}

    // Destructor.
    ~avl_grammar() {
      for (std::uint64_t i = 0; i < m_nonterminals.size(); ++i)
        delete m_nonterminals[i];
    }

    // Print the string encoded by the grammar.
    void print_expansion() const {
      m_root->print_expansion();
    }

    // Return the number of nonterminals.
    std::uint64_t size() const {
      return m_nonterminals.size();
    }

    // Add a nonterminal.
    void add_nonterminal(const node_type* nonterm) {
      m_nonterminals.push_back(nonterm);
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
    void collect_mersenne_karp_rabin_hashes(
        std::vector<std::uint64_t> &hashes) const {
      (void) m_root->collect_mersenne_karp_rabin_hashes(hashes);
    }

    // Collect Mersenne Karp-Rabin hashes in a hash table.
    void collect_mersenne_karp_rabin_hashes_2(
        hash_table<const node_type*, std::uint64_t> &hashes) const {
      (void) m_root->collect_mersenne_karp_rabin_hashes_2(hashes);
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

    // Add a nonterminal expanding to a substring of a given nonterminal.
    const node_type* add_substring_nonterminal(
        const node_type *x,
        const std::uint64_t begin,
        const std::uint64_t end) {
      std::vector<const node_type*> v = x->decomposition(begin, end);
      return greedy_merge(v);
    }

    // Add a substring expanding to a substring of root.
    const node_type* add_substring_nonterminal(
        const std::uint64_t begin,
        const std::uint64_t end) {
      return add_substring_nonterminal(m_root, begin, end);
    }

  private:

    // Merge greedily (shortest first) sequence of nonterminals.
    const node_type* greedy_merge(std::vector<const node_type*> &seq) {
      while (seq.size() > 1) {

        // Find the nonterminal with the smallest height.
        std::uint64_t smallest_height_id = 0;
        for (std::uint64_t i = 1; i < seq.size(); ++i) {
          if (seq[i]->m_height < seq[smallest_height_id]->m_height)
            smallest_height_id = i;
        }

        // Merge the nonterminal with the smaller height with
        // one of its beighbors (whichever is shorter).
        if (smallest_height_id == 0 ||
            (smallest_height_id + 1 < seq.size() &&
             seq[smallest_height_id + 1]->m_height <=
             seq[smallest_height_id - 1]->m_height)) {

          // Only right neighbor exists, or both exist
          // and the right one is not taller than the left
          // one. End result: merge with the right neighbor.
          const node_type * const left = seq[smallest_height_id];
          const node_type * const right = seq[smallest_height_id + 1];
          seq.erase(seq.begin() + smallest_height_id);
          seq[smallest_height_id] =
            add_concat_nonterminal<char_type>(m_nonterminals, left, right);
        } else {

          // Only left neighbor exists, or both exists
          // and the left one is not taller than the
          // right one. End result: merge with left neighbor.
          const node_type * const left = seq[smallest_height_id - 1];
          const node_type * const right = seq[smallest_height_id];
          seq.erase(seq.begin() + (smallest_height_id - 1));
          seq[smallest_height_id - 1] =
            add_concat_nonterminal<char_type>(m_nonterminals, left, right);
        }
      }
      return seq[0];
    }
};

#endif  // __AVL_GRAMMAR_HPP_INCLUDED
