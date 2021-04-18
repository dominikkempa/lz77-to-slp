#ifndef __AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
#define __AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <map>
#include <algorithm>

#include "../utils/hash_table.hpp"
#include "avl_grammar_node.hpp"
#include "avl_grammar_add_concat_nonterminal.hpp"
#include "dp_grouping_algorithm.hpp"


//=============================================================================
// A class storing multiroot AVL grammar.
//=============================================================================
template<typename char_type>
struct avl_grammar_multiroot {

  // Declare typedefs
  typedef avl_grammar_node<char_type> node_type;
  typedef std::uint64_t key_type;
  typedef const node_type* value_type;
  typedef typename std::map<key_type, value_type> map_type;
  typedef typename map_type::const_iterator const_iter_type;
  typedef typename map_type::iterator iter_type;

  private:

    // Class members.
    map_type m_roots;
    std::vector<const node_type*> m_nonterminals;
    hash_table<std::uint64_t, const node_type*> m_hashes;

  public:
  
    // Constructor.
    avl_grammar_multiroot() {
      m_roots.insert(std::make_pair(0, (value_type)NULL));
    }

    // Destructor.
    ~avl_grammar_multiroot() {
      for (std::uint64_t i = 0; i < m_nonterminals.size(); ++i)
        delete m_nonterminals[i];
    }

    // Print the string encoded by the grammar.
    void print_expansion() const {
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
        if ((std::uint64_t)it->first != 0)
          it->second->print_expansion();
    }

    // Return the number of nonterminals.
    std::uint64_t size() const {
      return m_nonterminals.size();
    }

    // Return the number of roots.
    std::uint64_t number_of_roots() const {
      return m_roots.size() - 1;
    }

    // Add a root.
    void add_root(
        const key_type pos,
        const value_type ptr) {
      m_roots[pos] = ptr;
    }

    // Add a nonterminal.
    void add_nonterminal(const node_type* nonterm) {
      m_nonterminals.push_back(nonterm);
    }

    // Decode the text and write to a given array.
    void decode(
        char_type* &text,
        std::uint64_t &text_length) const {
      text_length = 0;
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
        if ((std::uint64_t)it->first != 0)
          text_length += it->second->m_exp_len;
      text = new char_type[text_length];
      std::uint64_t ptr = 0;
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it) {
        if ((std::uint64_t)it->first != 0) {
          it->second->write_expansion(text + ptr);
          ptr += it->second->m_exp_len;
        }
      }
    }

    // Test the AVL property of all nonterminals.
    bool test_avl_property() const {
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
        if ((std::uint64_t)it->first != 0 &&
            it->second->test_avl_property() == false)
          return false;
      return true;
    }

    // Collect Mersenne Karp-Rabin hashes in a vector.
    void collect_mersenne_karp_rabin_hashes(
        std::vector<std::uint64_t> &hashes) const {
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
        if ((std::uint64_t)it->first != 0)
          (void) it->second->collect_mersenne_karp_rabin_hashes(hashes);
    }

    // Collect Mersenne Karp-Rabin hashes in a hash table.
    void collect_mersenne_karp_rabin_hashes_2(
        hash_table<const node_type*, std::uint64_t> &hashes) const {
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
        if ((std::uint64_t)it->first != 0)
          (void) it->second->collect_mersenne_karp_rabin_hashes_2(hashes);
    }

    // Count nodes in the pruned grammar.
    void count_nodes_in_pruned_grammar(
        hash_table<const node_type*, std::uint64_t> &hashes,
        hash_table<std::uint64_t, bool> &seen_hashes,
        std::uint64_t &current_count) const {
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
        if ((std::uint64_t)it->first != 0)
          it->second->count_nodes_in_pruned_grammar(
              hashes, seen_hashes, current_count);
    }

    // Collect pointers to all nonterminals reachable from the root.
    void collect_nonterminal_pointers(
        std::vector<const node_type*> &pointers) const {
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it)
        if ((std::uint64_t)it->first != 0)
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

      // Initialize length of previous nonterminals.
      std::uint64_t prev_nonterminals_total_length = 0;

      // Find range [it_begin..it_end) of roots to merge.
      iter_type it_begin = m_roots.end();
      it_begin = m_roots.lower_bound(begin);
      prev_nonterminals_total_length = it_begin->first;
      ++it_begin;

      // Now find it_end.
      iter_type it_end = it_begin;
      while (it_end != m_roots.end() && it_end->first <= end)
        ++it_end;

      // Merge nonterminals in [it_begin..it_end).
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
            const node_type * const left = v[smallest_height_id];
            const node_type * const right = v[smallest_height_id + 1];
            const std::uint64_t h = merge_hashes<char_type>(left, right);
            v.erase(v.begin() + smallest_height_id);
            if (m_hashes.find(h) != NULL)
              v[smallest_height_id] = *(m_hashes.find(h));
            else
              v[smallest_height_id] = add_concat_nonterminal<char_type>(
                  m_hashes, m_nonterminals, left, right);
          } else {

            // Only left neighbor exists, or both exists
            // and the left one is not taller than the
            // right one. End result: merge with left neighbor.
            const node_type * const left = v[smallest_height_id - 1];
            const node_type * const right = v[smallest_height_id];
            const std::uint64_t h = merge_hashes<char_type>(left, right);
            v.erase(v.begin() + (smallest_height_id - 1));
            if (m_hashes.find(h) != NULL)
              v[smallest_height_id - 1] = *(m_hashes.find(h));
            else
              v[smallest_height_id - 1] = add_concat_nonterminal<char_type>(
                  m_hashes, m_nonterminals, left, right);
          }
        }

        // Add the resulting nonterminal into the list of roots.
        value_type val = v[0];
        key_type key = prev_nonterminals_total_length + val->m_exp_len;
        m_roots[key] = val;
      }
    }

    // Get the sequence of nonterminals expanding to T[begin..end).
    std::vector<const node_type*> decomposition(
        std::uint64_t begin,
        std::uint64_t end) const {

      // Find leftmost root whose expansion overlaps/touches T[begin..end).
      const_iter_type it = m_roots.lower_bound(begin);

      // Proper substring or suffix of expansion of `it'.
      std::vector<const node_type*> ret;
      if (begin < it->first) {
        const std::uint64_t it_exp_size = it->second->m_exp_len;
        const std::uint64_t it_exp_beg = it->first - it_exp_size;
        const std::uint64_t local_beg = begin - it_exp_beg;
        const std::uint64_t local_end = std::min(it->first, end) - it_exp_beg;
        const std::uint64_t local_size = local_end - local_beg;
        std::vector<const node_type*> dec =
          it->second->decomposition(local_beg, local_end);
        ret.insert(ret.end(), dec.begin(), dec.end());
        begin += local_size;
      }

      // Full expansions of nonterminals.
      ++it;
      while (begin < end && it->first <= end) {
        ret.push_back(it->second);
        begin = it->first;
        ++it;
      }

      // Proper suffix of of expansion of `it'.
      if (begin < end) {
        const std::uint64_t it_exp_size = it->second->m_exp_len;
        const std::uint64_t it_exp_beg = it->first - it_exp_size;
        const std::uint64_t local_end = end - it_exp_beg;
        std::vector<const node_type*> dec =
          it->second->decomposition(0, local_end);
        ret.insert(ret.end(), dec.begin(), dec.end());
        begin = end;
      }

      // Return shorter equivalent sequence of nonterminals.
      std::vector<const node_type*> ret_opt =
        dp_grouping_algorithm<char_type>(m_hashes, ret);
      return ret_opt;
    }
};

#endif  // __AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
