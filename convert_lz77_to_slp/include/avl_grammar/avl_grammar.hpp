#ifndef __AVL_GRAMMAR_HPP_INCLUDED
#define __AVL_GRAMMAR_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "../utils/hash_table.hpp"
#include "nonterminal.hpp"


//=============================================================================
// A class storing AVL grammar.
//=============================================================================
template<typename char_type>
struct avl_grammar {
  private:

  //=========================================================================
  // Declare typedefs.
  //=========================================================================
  typedef nonterminal<char_type> nonterminal_type;

  private:

    //=========================================================================
    // Class members.
    //=========================================================================
    std::vector<const nonterminal_type*> m_nonterminals;
    const nonterminal_type *m_root;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    avl_grammar() :
      m_root(NULL) {}

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~avl_grammar() {
      for (std::uint64_t i = 0; i < m_nonterminals.size(); ++i)
        delete m_nonterminals[i];
    }

    //=========================================================================
    // Print the string encoded by the grammar.
    //=========================================================================
    void print_expansion() const {
      m_root->print_expansion();
    }

    //=========================================================================
    // Return the number of nonterminals.
    //=========================================================================
    std::uint64_t size() const {
      return m_nonterminals.size();
    }

    //=========================================================================
    // Set the root.
    //=========================================================================
    void set_root(
        const nonterminal_type * const newroot) {
      m_root = newroot;
    }

    //=========================================================================
    // Add a nonterminal.
    //=========================================================================
    void add_nonterminal(const nonterminal_type* nonterm) {
      m_nonterminals.push_back(nonterm);
    }

    //=========================================================================
    // Decode the text and write to a given array.
    //=========================================================================
    void decode(
        char_type* &text,
        std::uint64_t &text_length) const {
      text_length = m_root->m_exp_len;
      text = new char_type[text_length];
      m_root->write_expansion(text);
    }

    //=========================================================================
    // Test the AVL property of all nonterminals.
    //=========================================================================
    bool test_avl_property() const {
      return m_root->test_avl_property();
    }

    //=========================================================================
    // Collect Mersenne Karp-Rabin hashes in a vector.
    //=========================================================================
    void collect_mersenne_karp_rabin_hashes(
        std::vector<std::uint64_t> &hashes) const {
      (void) m_root->collect_mersenne_karp_rabin_hashes(hashes);
    }

    //=========================================================================
    // Collect Mersenne Karp-Rabin hashes in a hash table.
    //=========================================================================
    void collect_mersenne_karp_rabin_hashes_2(
        hash_table<const nonterminal_type*, std::uint64_t> &hashes) const {
      (void) m_root->collect_mersenne_karp_rabin_hashes_2(hashes);
    }

    //=========================================================================
    // Count nonterminals in the pruned grammar.
    //=========================================================================
    void count_nonterminals_in_pruned_grammar(
        hash_table<const nonterminal_type*, std::uint64_t> &hashes,
        hash_table<std::uint64_t, bool> &seen_hashes,
        std::uint64_t &current_count) const {
      m_root->count_nonterminals_in_pruned_grammar(hashes,
          seen_hashes, current_count);
    }

    //=========================================================================
    // Collect pointers to all nonterminals reachable from the root.
    //=========================================================================
    void collect_nonterminal_pointers(
        std::vector<const nonterminal_type*> &pointers) const {
      m_root->collect_nonterminal_pointers(pointers);
    }

    //=========================================================================
    // Add a nonterminal expanding to a substring of a given nonterminal.
    //=========================================================================
    const nonterminal_type* add_substring_nonterminal(
        const nonterminal_type *x,
        const std::uint64_t begin,
        const std::uint64_t end) {
      std::vector<const nonterminal_type*> v = x->decomposition(begin, end);
      return greedy_merge(v);
    }

    //=========================================================================
    // Add a substring expanding to a substring of grammar.
    //=========================================================================
    const nonterminal_type* add_substring_nonterminal(
        const std::uint64_t begin,
        const std::uint64_t end) {
      return add_substring_nonterminal(m_root, begin, end);
    }

    //=========================================================================
    // Given two nonterminals `left' and `right' expanding to X and Y, add
    // nonterminals that expands to XY, and return the pointer to it.
    //=========================================================================
    const nonterminal_type *add_concat_nonterminal(
        const nonterminal_type * const left,
        const nonterminal_type * const right) {
  
      // Consider two cases, depending on whether
      // left of right nonterminal is taller.
      if (left->m_height >= right->m_height) {
        if (left->m_height - right->m_height <= 1) {

          // Height are close. Just merge and return.
          const nonterminal_type * const newroot =
            new nonterminal_type(left, right);
          add_nonterminal(newroot);
          return newroot;
        } else {
          const nonterminal_type * const newright =
            add_concat_nonterminal(left->m_right, right);
          if (newright->m_height > left->m_left->m_height &&
              newright->m_height - left->m_left->m_height > 1) {

            // Rebalancing needed.
            if (newright->m_left->m_height > newright->m_right->m_height) {

              // Double (right-left) rotation.
              const nonterminal_type * const X = new nonterminal_type(
                  left->m_left, newright->m_left->m_left);
              const nonterminal_type * const Z = new nonterminal_type(
                  newright->m_left->m_right, newright->m_right);
              const nonterminal_type * const Y = new nonterminal_type(X, Z);
              add_nonterminal(X);
              add_nonterminal(Y);
              add_nonterminal(Z);
              return Y;
            } else {

              // Single (left) rotation.
              const nonterminal_type * const X =
                new nonterminal_type(left->m_left, newright->m_left);
              const nonterminal_type * const Y =
                new nonterminal_type(X, newright->m_right);
              add_nonterminal(X);
              add_nonterminal(Y);
              return Y;
            }
          } else {

            // No need to rebalance.
            const nonterminal_type * const newroot =
              new nonterminal_type(left->m_left, newright);
            add_nonterminal(newroot);
            return newroot;
          }
        }
      } else {
        if (right->m_height - left->m_height <= 1) {

          // Heights are close. Just merge and return.
          const nonterminal_type * const newroot =
            new nonterminal_type(left, right);
          add_nonterminal(newroot);
          return newroot;
        } else {
          const nonterminal_type * const newleft =
            add_concat_nonterminal(left, right->m_left);
          if (newleft->m_height > right->m_right->m_height &&
              newleft->m_height - right->m_right->m_height > 1) {

            // Rebalancing needed.
            if (newleft->m_right->m_height > newleft->m_left->m_height) {

              // Double (left-right) rotation.
              const nonterminal_type * const X = new nonterminal_type(
                  newleft->m_left, newleft->m_right->m_left);
              const nonterminal_type * const Z = new nonterminal_type(
                  newleft->m_right->m_right, right->m_right);
              const nonterminal_type * const Y = new nonterminal_type(X, Z);
              add_nonterminal(X);
              add_nonterminal(Y);
              add_nonterminal(Z);
              return Y;
            } else {

              // Single (right) rotation.
              const nonterminal_type * const Y =
                new nonterminal_type(newleft->m_right, right->m_right);
              const nonterminal_type * const X =
                new nonterminal_type(newleft->m_left, Y);
              add_nonterminal(X);
              add_nonterminal(Y);
              return X;
            }
          } else {

            // No need to rebalance.
            const nonterminal_type * const newroot =
              new nonterminal_type(newleft, right->m_right);
            add_nonterminal(newroot);
            return newroot;
          }
        }
      }
    }

  private:

    //=========================================================================
    // Merge greedily (shortest first) sequence of nonterminals.
    //=========================================================================
    const nonterminal_type* greedy_merge(
        std::vector<const nonterminal_type*> &seq) {
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
          const nonterminal_type * const left = seq[smallest_height_id];
          const nonterminal_type * const right = seq[smallest_height_id + 1];
          seq.erase(seq.begin() + smallest_height_id);
          seq[smallest_height_id] =
            add_concat_nonterminal(left, right);
        } else {

          // Only left neighbor exists, or both exists
          // and the left one is not taller than the
          // right one. End result: merge with left neighbor.
          const nonterminal_type * const left = seq[smallest_height_id - 1];
          const nonterminal_type * const right = seq[smallest_height_id];
          seq.erase(seq.begin() + (smallest_height_id - 1));
          seq[smallest_height_id - 1] =
            add_concat_nonterminal(left, right);
        }
      }
      return seq[0];
    }
};

#endif  // __AVL_GRAMMAR_HPP_INCLUDED
