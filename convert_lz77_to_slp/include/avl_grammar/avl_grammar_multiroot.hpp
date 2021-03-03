#ifndef __AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
#define __AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <map>
#include <algorithm>

#include "../utils/hash_table.hpp"
#include "../utils/karp_rabin_hashing.hpp"
#include "avl_grammar_node.hpp"
#include "avl_grammar_add_concat_nonterminal.hpp"
#include "avl_grammar_add_substring_nonterminal.hpp"
#include "avl_grammar_decomposition.hpp"
#include "dp_grouping_algorithm.hpp"


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
  hash_table<std::uint64_t, const node_type*> m_hashes;
  std::vector<const node_type*> m_nonterminals;
  map_type m_roots;
  const std::uint64_t m_hash_variable;
  const std::uint64_t m_mersenne_prime_exponent;

  // Constructor.
  avl_grammar_multiroot(
      const std::uint64_t hash_variable,
      const std::uint64_t mersenne_prime_exponent) :
    m_hash_variable(hash_variable),
    m_mersenne_prime_exponent(mersenne_prime_exponent) {}

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
          std::uint64_t h_left = v[smallest_height_id]->m_kr_hash;
          std::uint64_t h_right = v[smallest_height_id + 1]->m_kr_hash;
          std::uint64_t len_right = v[smallest_height_id + 1]->m_exp_len;
          std::uint64_t h =
            mod_mersenne(
                mul_mod_meresenne(
                  h_left,
                  pow_mod_mersenne(
                    m_hash_variable,
                    len_right,
                    m_mersenne_prime_exponent),
                  m_mersenne_prime_exponent) + h_right,
                m_mersenne_prime_exponent);
          if (m_hashes.find(h) != NULL) {
            v.erase(v.begin() + smallest_height_id);
            v[smallest_height_id] = *(m_hashes.find(h));
          } else {
            value_type ret =
              add_concat_nonterminal<char_type>(m_hashes, m_nonterminals,
                  v[smallest_height_id], v[smallest_height_id + 1],
                  m_hash_variable, m_mersenne_prime_exponent);
            v.erase(v.begin() + smallest_height_id);
            v[smallest_height_id] = ret;
          }
        } else {

          // Only left neighbor exists, or both exists
          // and the left one is not taller than the
          // right one. End result: merge with left neighbor.
          std::uint64_t h_left = v[smallest_height_id - 1]->m_kr_hash;
          std::uint64_t h_right = v[smallest_height_id]->m_kr_hash;
          std::uint64_t len_right = v[smallest_height_id]->m_exp_len;
          std::uint64_t h =
            mod_mersenne(
                mul_mod_meresenne(
                  h_left,
                  pow_mod_mersenne(
                    m_hash_variable,
                    len_right,
                    m_mersenne_prime_exponent),
                  m_mersenne_prime_exponent) + h_right,
                m_mersenne_prime_exponent);
          if (m_hashes.find(h) != NULL) {
            v.erase(v.begin() + (smallest_height_id - 1));
            v[smallest_height_id - 1] = *(m_hashes.find(h));
          } else {
            value_type ret =
              add_concat_nonterminal<char_type>(m_hashes, m_nonterminals,
                  v[smallest_height_id - 1], v[smallest_height_id],
                  m_hash_variable, m_mersenne_prime_exponent);
            v.erase(v.begin() + (smallest_height_id - 1));
            v[smallest_height_id - 1] = ret;
          }
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
  std::vector<const node_type*> add_substring_nonterminal(
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

    // Consider all possible cases in which the block
    // [begin..end) can overlap the roots of the grammar.
    std::vector<const node_type*> ret_vec;
    if (begin == 0 || it->first == begin) {

      // If the block [begin..end) starts at the beginning
      // of the expansion of the root in the grammar, there
      // are three cases to consider. First, obtain the
      // point to the grammar root whose expansion overlaps
      // block [begin..end).
      iter_type it_left = it;
      if (begin > 0)
        ++it_left;

      // Consider all three cases.
      if (it_left->first <= end) {
        if (it_left->first == end) {

          // Case I: the block [begin..end) is equal
          // to the expansion of the grammar root. No
          // need to merge anything.
          ret_vec.push_back(it_left->second);
        } else {

          // Case II: it_left->first > end. We need
          // to create a nonterminal expanding to some
          // propert prefix of the expansion of the right
          // neighbor of it_left->second, and then merge
          // it with it_left->second.
          const node_type * const left = it_left->second;
          ret_vec.push_back(left);
          const std::uint64_t rbegin = 0;
          const std::uint64_t rend = end - it_left->first;
          iter_type it_right = it_left;
          ++it_right;
          //{
            //const node_type * const right =
            //  ::add_substring_nonterminal<char_type>(
            //      m_hashes, m_nonterminals, it_right->second, rbegin, rend,
            //      m_hash_variable, m_mersenne_prime_exponent);
            //ret_vec.push_back(right);
          //}
          {
            std::vector<const node_type*> right_decomposition =
              decomposition<char_type>(it_right->second, rbegin, rend);
            std::vector<const node_type*> right_opt_decomposition =
              //dp_grouping_algorithm<char_type>(m_hashes, right_decomposition,
              //    m_hash_variable, m_mersenne_prime_exponent);
              right_decomposition;
            ret_vec.insert(ret_vec.end(),
                right_opt_decomposition.begin(),
                right_opt_decomposition.end());
          }
        }
      } else {

        // Case III: block [begin..end) is entirely inside
        // the expansion of it_left->second. It suffices to
        // call add_substring_nonterminal.
        const std::uint64_t lbegin = 0;
        const std::uint64_t lend = end - begin;
        //{
        //  const node_type * const ret =
        //    ::add_substring_nonterminal<char_type>(
        //        m_hashes, m_nonterminals, it_left->second, lbegin, lend,
        //        m_hash_variable, m_mersenne_prime_exponent);
        //  ret_vec.push_back(ret);
        //}
        {
          std::vector<const node_type*> ret_decomposition =
            decomposition<char_type>(it_left->second, lbegin, lend);
          std::vector<const node_type*> ret_opt_decomposition =
            //dp_grouping_algorithm<char_type>(m_hashes, ret_decomposition,
            //    m_hash_variable, m_mersenne_prime_exponent);
            ret_decomposition;
          ret_vec.insert(ret_vec.end(),
              ret_opt_decomposition.begin(),
              ret_opt_decomposition.end());
        }
      }
    } else {

      // It it->first > begin, then there are two subcases to
      // consider. Either the block [begin..end) it entirely
      // contain in the expansion of it->second, or not.
      if (end <= it->first) {

        // Case I: the block [begin..end) it entirely contained
        // in the expansion of nonterminal it->second. Since
        // we separately considered the case of begin = 0 before,
        // we are not therefore quaranteed that it is moreover
        // a proper substring.
        const std::uint64_t it_exp_size = it->second->m_exp_len;
        const std::uint64_t it_exp_beg = it->first - it_exp_size;
        const std::uint64_t lbegin = begin - it_exp_beg;
        const std::uint64_t lend = end - it_exp_beg;
        //{
        //  const node_type * const ret =
        //    ::add_substring_nonterminal<char_type>(
        //        m_hashes, m_nonterminals, it->second, lbegin, lend,
        //        m_hash_variable, m_mersenne_prime_exponent);
        //  ret_vec.push_back(ret);
        //}
        {
          std::vector<const node_type*> ret_decomposition =
            decomposition<char_type>(it->second, lbegin, lend);
          std::vector<const node_type*> ret_opt_decomposition =
            //dp_grouping_algorithm<char_type>(m_hashes, ret_decomposition,
            //    m_hash_variable, m_mersenne_prime_exponent);
            ret_decomposition;
          ret_vec.insert(ret_vec.end(),
              ret_opt_decomposition.begin(),
              ret_opt_decomposition.end());
        }
      } else {

        // Case II: the expansion of it->first and block [begin..end)
        // overlap each other in the proper way (i.e., neither is
        // contained inside the other one). Thus, we can already compute
        // the nonterminal expanding to the proper prefix of [begin..end).
        const std::uint64_t it_block_size = it->second->m_exp_len;
        const std::uint64_t lsize = it->first - begin;
        const std::uint64_t lbegin = it_block_size - lsize;
        const std::uint64_t lend = it_block_size;
        //{
        //  const node_type * const left_nonterminal =
        //    ::add_substring_nonterminal<char_type>(
        //        m_hashes, m_nonterminals, it->second, lbegin, lend,
        //        m_hash_variable, m_mersenne_prime_exponent);
        //  ret_vec.push_back(left_nonterminal);
        //}
        {
          std::vector<const node_type*> left_decomposition =
            decomposition<char_type>(it->second, lbegin, lend);
          std::vector<const node_type*> left_opt_decomposition =
            //dp_grouping_algorithm<char_type>(m_hashes, left_decomposition,
            //    m_hash_variable, m_mersenne_prime_exponent);
            left_decomposition;
          ret_vec.insert(ret_vec.end(),
              left_opt_decomposition.begin(),
              left_opt_decomposition.end());
        }

        // We now have three cases: either the block [begin..end)
        // overlaps expansions of two grammar roots (and there there
        // are two subcases, depending, on whether it perfectly ends
        // at the end of the right one or not), or three. If three,
        // then we are guaranteed that for the third one, we only
        // overlap by the proper prefix.
        iter_type it_mid = it;
        ++it_mid;
        if (end < it_mid->first) {

          // Case IIa: the block [begin..end) properly overlaps
          // the expansion of it_mid->second. We first create the
          // nonterminal expanding to that prefix, and them merge
          // it with left_nonterminal.
          const std::uint64_t mbegin = 0;
          const std::uint64_t mend = end - it->first;
          //{
          //  const node_type * const mid_nonterminal =
          //    ::add_substring_nonterminal<char_type>(
          //        m_hashes, m_nonterminals, it_mid->second, mbegin, mend,
          //        m_hash_variable, m_mersenne_prime_exponent);
          //  ret_vec.push_back(mid_nonterminal);
          //}
          {
            std::vector<const node_type*> mid_decomposition =
              decomposition<char_type>(it_mid->second, mbegin, mend);
            std::vector<const node_type*> mid_opt_decomposition =
              //dp_grouping_algorithm<char_type>(m_hashes, mid_decomposition,
              //    m_hash_variable, m_mersenne_prime_exponent);
              mid_decomposition;
            ret_vec.insert(ret_vec.end(),
                mid_opt_decomposition.begin(),
                mid_opt_decomposition.end());
          }
        } else if (end == it_mid->first) {

          // Case IIb: the block [begin..end) ends at the
          // expansion of it_mid->second. Thus, ut suffices
          // to merge left_nonterminal with it_mid->second.
          ret_vec.push_back(it_mid->second);
        } else {

          // Case IIc: the block [begin..end) overlaps (as
          // a proper prefix) the expansion of the grammar
          // root to the right of it_mid. First, create the
          // nonterminal expanding to that prefix.
          ret_vec.push_back(it_mid->second);
          iter_type it_right = it_mid;
          ++it_right;
          const std::uint64_t rbegin = 0;
          const std::uint64_t rend = end - it_mid->first;
          //{
          //  const node_type * const right_nonterminal =
          //    ::add_substring_nonterminal<char_type>(
          //        m_hashes, m_nonterminals, it_right->second, rbegin, rend,
          //        m_hash_variable, m_mersenne_prime_exponent);
          //  ret_vec.push_back(right_nonterminal);
          //}
          {
            std::vector<const node_type*> right_decomposition =
              decomposition<char_type>(it_right->second, rbegin, rend);
            std::vector<const node_type*> right_opt_decomposition =
              //dp_grouping_algorithm<char_type>(m_hashes, right_decomposition,
              //    m_hash_variable, m_mersenne_prime_exponent);
              right_decomposition;
            ret_vec.insert(ret_vec.end(),
                right_opt_decomposition.begin(),
                right_opt_decomposition.end());
          }
        }
      }
    }
    {
      std::vector<const node_type*> ret_vec_opt =
        dp_grouping_algorithm<char_type>(m_hashes, ret_vec,
            m_hash_variable, m_mersenne_prime_exponent);
      return ret_vec_opt;
    }
    //return ret_vec;
  }
};

#endif  // __AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
