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
#include "avl_grammar_add_substring_nonterminal.hpp"


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
    hash_table<std::uint64_t, const node_type*> m_hashes;

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
};

#endif  // __AVL_GRAMMAR_HPP_INCLUDED
