/**
 * @file    lazy_avl_grammar.hpp
 * @section LICENCE
 *
 * This file is part of Lazy-AVLG v0.1.0
 * See: https://github.com/dominikkempa/lz77-to-slp
 *
 * Copyright (C) 2021
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __LAZY_AVL_GRAMMAR_HPP_INCLUDED
#define __LAZY_AVL_GRAMMAR_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "../utils/hash_table.hpp"
#include "../utils/karp_rabin_hashing.hpp"
#include "../utils/space_efficient_vector.hpp"
#include "../utils/cache.hpp"
#include "../utils/packed_pair.hpp"
#include "../utils/packed_triple.hpp"
#include "../io/async_stream_writer.hpp"


//=============================================================================
// Class used to represent the lazy AVL grammar. Forward declaration.
//=============================================================================
template<typename char_type, typename text_offset_type>
struct lazy_avl_grammar;

//=============================================================================
// Class used to represent nonterminal. Forward declaration.
//=============================================================================
template<typename char_type, typename text_offset_type>
struct nonterminal {
  private:

    //=========================================================================
    // Declare types.
    //=========================================================================
    typedef nonterminal<char_type, text_offset_type> nonterminal_type;
    typedef text_offset_type ptr_type;
    typedef lazy_avl_grammar<char_type, text_offset_type> grammar_type;
    typedef packed_pair<text_offset_type, text_offset_type> pair_type;

    //=========================================================================
    // Class members.
    //=========================================================================
    std::uint8_t m_height;
    std::uint8_t m_truncated_exp_len;
    ptr_type m_left_p;
    ptr_type m_right_p;

  public:

    //=========================================================================
    // Constructors.
    //=========================================================================
    nonterminal();
    nonterminal(const char_type);
    nonterminal(const std::uint8_t, const std::uint8_t,
      const ptr_type, const ptr_type);

    //=========================================================================
    // Access methods.
    //=========================================================================
    std::uint64_t get_height() const;
    std::uint64_t get_truncated_exp_len() const;
    char_type get_char() const;
    ptr_type get_left_p() const;
    ptr_type get_right_p() const;

    //=========================================================================
    // Key methods.
    //=========================================================================
    void decomposition(ptr_type, const std::uint64_t,
        const std::uint64_t, space_efficient_vector<pair_type> &,
        const grammar_type * const g) const;
    std::uint64_t write_expansion(const ptr_type, char_type * const,
        const grammar_type * const) const;

    //=========================================================================
    // Mostly unused.
    //=========================================================================
    void print_expansion(const ptr_type,
        const grammar_type * const) const;
    bool compare_expansion_to_text(const ptr_type,
        const char_type * const, const grammar_type * const) const;
    bool test_avl_property(const ptr_type,
        const grammar_type * const) const;
    std::uint64_t collect_mersenne_karp_rabin_hashes(const ptr_type,
        std::vector<std::uint64_t> &, const grammar_type * const) const;
    void collect_nonterminal_pointers(const ptr_type,
        std::vector<ptr_type> &, const grammar_type * const) const;
    std::uint64_t collect_mersenne_karp_rabin_hashes_2(const ptr_type,
        hash_table<ptr_type, std::uint64_t> &,
        const grammar_type * const) const;
    void count_nonterminals_in_pruned_grammar(const ptr_type,
        hash_table<ptr_type, std::uint64_t> &,
        hash_table<std::uint64_t, bool> &,
        std::uint64_t &, const grammar_type * const) const;
} __attribute__((packed));

//=============================================================================
// Hash functions of the appropriate type.
// Used in the hash table used to prune the grammar.
//=============================================================================
typedef const uint40 get_hash_ptr_type;

template<>
std::uint64_t get_hash(const get_hash_ptr_type &x) {
  return (std::uint64_t)x * (std::uint64_t)29996224275833;
}

template<>
std::uint64_t get_hash(const std::uint64_t &x) {
  return (std::uint64_t)x * (std::uint64_t)4972548694736365;
}

//=============================================================================
// Implementation of the lazy_avl_grammar class.
//=============================================================================
template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
struct lazy_avl_grammar {
  static_assert(sizeof(char_type) <= sizeof(text_offset_type),
      "Error: sizeof(char_type) > sizeof(text_offset_type)!");

  private:

    //=========================================================================
    // Declare typedefs.
    //=========================================================================
    typedef nonterminal<char_type, text_offset_type> nonterminal_type;
    typedef text_offset_type ptr_type;
    typedef packed_pair<text_offset_type, text_offset_type> pair_type;
    typedef packed_triple<text_offset_type, text_offset_type, std::uint64_t>
      triple_type;
    typedef packed_pair<text_offset_type, std::uint64_t> hash_pair_type;

    //=========================================================================
    // Class members.
    //=========================================================================
    space_efficient_vector<nonterminal_type> m_nonterminals;
    space_efficient_vector<pair_type> m_roots_vec;
    space_efficient_vector<pair_type> m_long_exp_len;
    space_efficient_vector<hash_pair_type> m_long_exp_hashes;
    hash_table<std::uint64_t, ptr_type, text_offset_type> m_hashes;
    cache<ptr_type, std::uint64_t> *m_kr_hash_cache;
    char_type *m_snippet;
    std::uint64_t m_empty_step_counter;
    bool m_enable_kr_hashing;
    long double m_kr_hashing_prob;

    //=========================================================================
    // Internal statistics.
    //=========================================================================
    std::uint64_t m_merge_count;
    std::uint64_t m_avoided_merges;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    lazy_avl_grammar(
        bool enable_kr_hashing,
        long double kr_hashing_prob = (long double)0.0) {
      m_merge_count = 0;
      m_avoided_merges = 0;
      const std::uint64_t cache_size = (1 << 14);

      // Initialize standard members.
      m_empty_step_counter = 0;
      m_enable_kr_hashing = enable_kr_hashing;
      m_kr_hashing_prob = kr_hashing_prob;
      m_roots_vec.push_back(pair_type(0, std::numeric_limits<ptr_type>::max()));

      // Initialize KR hashing (optional).
      if (m_enable_kr_hashing) {
        m_kr_hash_cache = new cache<ptr_type, std::uint64_t>(cache_size);
        m_snippet = utils::allocate_array<char_type>(256);
      } else {
        m_kr_hash_cache = NULL;
        m_snippet = NULL;
      }
    }

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~lazy_avl_grammar() {
      if (m_enable_kr_hashing) {
        utils::deallocate(m_snippet);
        delete m_kr_hash_cache;
      }
    }

    //=========================================================================
    // Print the string encoded by the grammar.
    //=========================================================================
    void print_expansion() {
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i)) {
        const std::uint64_t preflen = m_roots_vec[i].first;
        const ptr_type id = m_roots_vec[i].second;
        if (preflen != 0) {
          const nonterminal_type &nonterm = get_nonterminal(id);
          nonterm.print_expansion(id, this);
        }
      }
    }

    //=========================================================================
    // Return the number of nonterminals.
    //=========================================================================
    std::uint64_t size() const {
      return m_nonterminals.size();
    }

    //=========================================================================
    // Return the number of roots.
    //=========================================================================
    std::uint64_t number_of_roots() {
      std::uint64_t ret = 0;
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i))
        ++ret;
      return ret - 1;
    }

    //=========================================================================
    // Add a new root at the end of the sequence.
    //=========================================================================
    void push_root(
        const std::uint64_t pos,
        const ptr_type id) {
      m_roots_vec.push_back(pair_type(pos, id));
    }

    //=========================================================================
    // Returns true if the given nonterminal is a liteval (i.e., exp len = 1).
    //========================================================================
    bool is_literal(const ptr_type p) const {
      const nonterminal_type &nonterm = get_nonterminal(p);
      const std::uint64_t height = nonterm.get_height();
      return (height == 0);
    }

    //=========================================================================
    // Return the total length of right-hand sides of all expansions.
    //=========================================================================
    std::uint64_t total_rhs_length() const {
      std::uint64_t ret = 0;
      for (std::uint64_t i = 0; i < m_nonterminals.size(); ++i) {
        if (is_literal(i))
          ret += 1;
        else ret += 2;
      }
      return ret;
    }

    //=========================================================================
    // Gives access to a given nonterminal.
    //=========================================================================
    const nonterminal_type& get_nonterminal(const ptr_type id) const {
      return m_nonterminals[(std::uint64_t)id];
    }

    //=========================================================================
    // Return the expansion length of a given nonterminal.
    //=========================================================================
    std::uint64_t get_exp_len(const ptr_type id) const {
      const nonterminal_type &nonterm = get_nonterminal(id);
      std::uint64_t exp_len = 0;
      const std::uint64_t truncated_exp_len = nonterm.get_truncated_exp_len();
      if (truncated_exp_len == 255) {

        // Binary search in m_long_exp_len.
        std::uint64_t beg = 0;
        std::uint64_t end = m_long_exp_len.size();
        while (beg + 1 < end) {
          const std::uint64_t mid = (beg + end) / 2;
          const std::uint64_t mid_id = m_long_exp_len[mid].first;
          if (mid_id <= (std::uint64_t)id)
            beg = mid;
          else end = mid;
        }
        exp_len = m_long_exp_len[beg].second;
      } else exp_len = truncated_exp_len;
      return exp_len;
    }

    //=========================================================================
    // Return the Karp-Rabin hash of a given nonterminal.
    //=========================================================================
    std::uint64_t get_kr_hash(const ptr_type id) const {

      // Check, if the value is in cache.
      const std::uint64_t *cache_ret =
        m_kr_hash_cache->lookup(id);
      if (cache_ret != NULL)
        return *cache_ret;

      // Obtain/compute the hash.
      const nonterminal_type &nonterm = get_nonterminal(id);
      const std::uint64_t truncated_exp_len =
        nonterm.get_truncated_exp_len();
      std::uint64_t kr_hash = 0;
      if (truncated_exp_len < 255) {

        // Recompute the hash from scratch.
        (void) nonterm.write_expansion(id, m_snippet, this);
        kr_hash = karp_rabin_hashing::hash_string<char_type>(
            m_snippet, truncated_exp_len);
      } else {

        // Binary search in m_long_exp_hashes.
        std::uint64_t beg = 0;
        std::uint64_t end = m_long_exp_hashes.size();
        while (beg + 1 < end) {
          const std::uint64_t mid = (beg + end) / 2;
          const std::uint64_t mid_id = m_long_exp_hashes[mid].first;
          if (mid_id <= id)
            beg = mid;
          else end = mid;
        }
        kr_hash = m_long_exp_hashes[beg].second;
      }

      // Update cache.
      m_kr_hash_cache->insert(id, kr_hash);

      // Return the result.
      return kr_hash;
    }

    //=========================================================================
    // Add a nonterminal.
    //=========================================================================
    ptr_type add_nonterminal(const nonterminal_type &nonterm) {
      const ptr_type new_nonterm_p = m_nonterminals.size();
      m_nonterminals.push_back(nonterm);

      // With some probability add to hash table.
      if (m_enable_kr_hashing &&
          utils::random_real() < m_kr_hashing_prob) {
        const std::uint64_t kr_hash = get_kr_hash(new_nonterm_p);
        ptr_type * const ret = m_hashes.find(kr_hash);
        if (ret == NULL)
          m_hashes.insert(kr_hash, new_nonterm_p);
        else *ret = new_nonterm_p;
      }

      // Return the id of the nonterminal.
      return new_nonterm_p;
    }

    //=========================================================================
    // Add a new binary nonterminal.
    //=========================================================================
    ptr_type add_nonterminal(
        const ptr_type left_p,
        const ptr_type right_p) {
      const nonterminal_type &left = get_nonterminal(left_p);
      const nonterminal_type &right = get_nonterminal(right_p);

      // Compute values for the new nonterminal.
      const std::uint64_t left_exp_len = get_exp_len(left_p);
      const std::uint64_t right_exp_len = get_exp_len(right_p);
      const std::uint64_t exp_len = left_exp_len + right_exp_len;
      const std::uint64_t truncated_exp_len =
        std::min((std::uint64_t)255, exp_len);
      const std::uint8_t left_height = left.get_height();
      const std::uint8_t right_height = right.get_height();
      const std::uint8_t height = std::max(left_height, right_height) + 1;

      // Create and add new nonterminal.
      const ptr_type new_nonterm_p = m_nonterminals.size();
      nonterminal_type new_nonterm(height, truncated_exp_len, left_p, right_p);
      m_nonterminals.push_back(new_nonterm);

      // With some probability add to hash table.
      std::uint64_t kr_hash = 0;
      bool hash_computed = false;
      if (m_enable_kr_hashing &&
          utils::random_real() < m_kr_hashing_prob) {
        const std::uint64_t left_hash = get_kr_hash(left_p);
        const std::uint64_t right_hash = get_kr_hash(right_p);
        kr_hash = karp_rabin_hashing::concat(
            left_hash, right_hash, right_exp_len);
        hash_computed = true;
        ptr_type * const ret = m_hashes.find(kr_hash);
        if (ret == NULL)
          m_hashes.insert(kr_hash, new_nonterm_p);
        else *ret = new_nonterm_p;
      }

      // Update list of long nonterminals.
      if (exp_len >= 255)
        m_long_exp_len.push_back(pair_type(new_nonterm_p, exp_len));

      // Update list of hashes for long nonterminals.
      if (m_enable_kr_hashing && exp_len >= 255) {
        if (!hash_computed) {
          const std::uint64_t left_hash = get_kr_hash(left_p);
          const std::uint64_t right_hash = get_kr_hash(right_p);
          kr_hash = karp_rabin_hashing::concat(
              left_hash, right_hash, right_exp_len);
          hash_computed = true;
        }
        m_long_exp_hashes.push_back(
            hash_pair_type(new_nonterm_p, kr_hash));
      }

      // Return the ptr to the new nonterminal.
      return new_nonterm_p;
    }

    //=========================================================================
    // Decode the text and write to a given array.
    //=========================================================================
    void decode(
        char_type* &text,
        std::uint64_t &text_length) {
      text_length = 0;
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i)) {
        const std::uint64_t preflen = m_roots_vec[i].first;
        const ptr_type id = m_roots_vec[i].second;
        if (preflen != 0)
          text_length += get_exp_len(id);
      }
      text = new char_type[text_length];
      std::uint64_t ptr = 0;
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i)) {
        const std::uint64_t preflen = m_roots_vec[i].first;
        const ptr_type id = m_roots_vec[i].second;
        if (preflen != 0) {
          const nonterminal_type &nonterm = get_nonterminal(id);
          const std::uint64_t exp_len =
            nonterm.write_expansion(id, text + ptr, this);
          ptr += exp_len;
        }
      }
    }

    //=========================================================================
    // Check if the grammar expands to the given string.
    //=========================================================================
    bool compare_expansion_to_text(
        const char_type * const text,
        const std::uint64_t text_length) {

      // Compute length of expansion.
      std::uint64_t exp_total_len = 0;
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i)) {
        const std::uint64_t preflen = m_roots_vec[i].first;
        const ptr_type id = m_roots_vec[i].second;
        if (preflen != 0)
          exp_total_len += get_exp_len(id);
      }

      // If they differ, return false.
      if (text_length != exp_total_len)
        return false;

      // Otherwise, compare the generated string with `text'.
      std::uint64_t ptr = 0;
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i)) {
        const std::uint64_t preflen = m_roots_vec[i].first;
        const ptr_type id = m_roots_vec[i].second;
        if (preflen != 0) {
          const nonterminal_type &nonterm = get_nonterminal(id);
          if (!nonterm.compare_expansion_to_text(id, text + ptr, this))
            return false;
          ptr += get_exp_len(id);
        }
      }
      return true;
    }

    //=========================================================================
    // Test the AVL property of all nonterminals.
    //=========================================================================
    bool test_avl_property() {
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i)) {
        const std::uint64_t preflen = m_roots_vec[i].first;
        const std::uint64_t id = m_roots_vec[i].second;
        if (preflen != 0) {
          const nonterminal_type &nonterm = get_nonterminal(id);
          if (nonterm.test_avl_property(id, this) == false)
            return false;
        }
      }
      return true;
    }

    //=========================================================================
    // Collect Mersenne Karp-Rabin hashes in a vector.
    //=========================================================================
    void collect_mersenne_karp_rabin_hashes(
        std::vector<std::uint64_t> &hashes) {
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i)) {
        const std::uint64_t preflen = m_roots_vec[i].first;
        const std::uint64_t id = m_roots_vec[i].second;
        if (preflen != 0) {
          const nonterminal_type &nonterm = get_nonterminal(id);
          (void) nonterm.collect_mersenne_karp_rabin_hashes(id, hashes, this);
        }
      }
    }

    //=========================================================================
    // Collect Mersenne Karp-Rabin hashes in a hash table.
    //=========================================================================
    void collect_mersenne_karp_rabin_hashes_2(
        hash_table<ptr_type, std::uint64_t> &hashes) {
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i)) {
        const std::uint64_t preflen = m_roots_vec[i].first;
        const std::uint64_t id = m_roots_vec[i].second;
        if (preflen != 0) {
          const nonterminal_type &nonterm = get_nonterminal(id);
          (void) nonterm.collect_mersenne_karp_rabin_hashes_2(
              id, hashes, this);
        }
      }
    }

    //=========================================================================
    // Count nonterminals in the pruned grammar.
    //=========================================================================
    void count_nonterminals_in_pruned_grammar(
        hash_table<ptr_type, std::uint64_t> &hashes,
        hash_table<std::uint64_t, bool> &seen_hashes,
        std::uint64_t &current_count) {
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i)) {
        const std::uint64_t preflen = m_roots_vec[i].first;
        const std::uint64_t id = m_roots_vec[i].second;
        if (preflen != 0) {
          const nonterminal_type &nonterm = get_nonterminal(id);
          nonterm.count_nonterminals_in_pruned_grammar(
            id, hashes, seen_hashes, current_count, this);
        }
      }
    }

    //=========================================================================
    // Collect pointers to all nonterminals reachable from the root.
    //=========================================================================
    void collect_nonterminal_pointers(
        std::vector<ptr_type> &pointers) {
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i)) {
        const std::uint64_t preflen = m_roots_vec[i].first;
        const std::uint64_t id = m_roots_vec[i].second;
        if (preflen != 0) {
          const nonterminal_type &nonterm = get_nonterminal(id);
          nonterm.collect_nonterminal_pointers(id, pointers, this);
        }
      }
    }

    //=========================================================================
    // Given two nonterminals `left' and `right' expanding to X and Y, add
    // nonterminals that expands to XY, and return the pointer to it.
    //=========================================================================
    ptr_type add_concat_nonterminal(
        const ptr_type left_p,
        const ptr_type right_p) {

      // Consider two cases, depending on whether
      // left or right nonterminal is taller.
      const nonterminal_type &left = get_nonterminal(left_p);
      const nonterminal_type &right = get_nonterminal(right_p);
      if (left.get_height() >= right.get_height()) {
        if (left.get_height() - right.get_height() <= 1) {

          // Height are close. Just merge and return.
          const ptr_type newroot_p = add_nonterminal(left_p, right_p);
          return newroot_p;
        } else {
          const ptr_type leftleft_p = left.get_left_p();
          const ptr_type leftright_p = left.get_right_p();

          // Careful here: add_concat_nonterminal may reallocate
          // m_nonterminals, which may invalidate refs to m_nonterminals.
          // This was a rather nasty bug, that took a while to find.
          const ptr_type newright_p = add_concat_nonterminal(leftright_p, right_p);
          const nonterminal_type &leftleft = get_nonterminal(leftleft_p);
          const nonterminal_type &newright = get_nonterminal(newright_p);
          if (newright.get_height() > leftleft.get_height() &&
              newright.get_height() - leftleft.get_height() > 1) {

            // Rebalancing needed.
            const ptr_type newright_left_p = newright.get_left_p();
            const ptr_type newright_right_p = newright.get_right_p();
            const nonterminal_type &newright_left = get_nonterminal(newright_left_p);
            const nonterminal_type &newright_right = get_nonterminal(newright_right_p);
            if (newright_left.get_height() > newright_right.get_height()) {

              // Double (right-left) rotation. Be careful also here:
              // add_nonterminal can also invalidate refs to m_nonterminals.
              const ptr_type newright_leftleft_p = newright_left.get_left_p();
              const ptr_type newright_leftright_p = newright_left.get_right_p();
              const ptr_type X_p = add_nonterminal(leftleft_p, newright_leftleft_p);
              const ptr_type Z_p = add_nonterminal(newright_leftright_p, newright_right_p);
              const ptr_type Y_p = add_nonterminal(X_p, Z_p);
              return Y_p;
            } else {

              // Single (left) rotation.
              const ptr_type X_p = add_nonterminal(leftleft_p, newright_left_p);
              const ptr_type Y_p = add_nonterminal(X_p, newright_right_p);
              return Y_p;
            }
          } else {

            // No need to rebalance.
            const ptr_type newroot_p = add_nonterminal(leftleft_p, newright_p);
            return newroot_p;
          }
        }
      } else {
        if (right.get_height() - left.get_height() <= 1) {

          // Heights are close. Just merge and return.
          const ptr_type newroot_p = add_nonterminal(left_p, right_p);
          return newroot_p;
        } else {
          const ptr_type rightleft_p = right.get_left_p();
          const ptr_type rightright_p = right.get_right_p();
          const ptr_type newleft_p = add_concat_nonterminal(left_p, rightleft_p);
          const nonterminal_type &rightright = get_nonterminal(rightright_p);
          const nonterminal_type &newleft = get_nonterminal(newleft_p);
          if (newleft.get_height() > rightright.get_height() &&
              newleft.get_height() - rightright.get_height() > 1) {

            // Rebalancing needed.
            const ptr_type newleft_left_p = newleft.get_left_p();
            const ptr_type newleft_right_p = newleft.get_right_p();
            const nonterminal_type &newleft_left = get_nonterminal(newleft_left_p);
            const nonterminal_type &newleft_right = get_nonterminal(newleft_right_p);
            if (newleft_right.get_height() > newleft_left.get_height()) {

              // Double (left-right) rotation.
              const ptr_type newleft_rightleft_p = newleft_right.get_left_p();
              const ptr_type newleft_rightright_p = newleft_right.get_right_p();
              const ptr_type X_p = add_nonterminal(newleft_left_p, newleft_rightleft_p);
              const ptr_type Z_p = add_nonterminal(newleft_rightright_p, rightright_p);
              const ptr_type Y_p = add_nonterminal(X_p, Z_p);
              return Y_p;
            } else {

              // Single (right) rotation.
              const ptr_type Y_p = add_nonterminal(newleft_right_p, rightright_p);
              const ptr_type X_p = add_nonterminal(newleft_left_p, Y_p);
              return X_p;
            }
          } else {

            // No need to rebalance.
            const ptr_type newroot_p = add_nonterminal(newleft_p, rightright_p);
            return newroot_p;
          }
        }
      }
    }

    //=========================================================================
    // Return RAM use.
    //=========================================================================
    std::uint64_t ram_use() const {
      const std::uint64_t m_roots_vec_ram_use = m_roots_vec.ram_use();
      const std::uint64_t m_nonterminals_ram_use = m_nonterminals.ram_use();
      const std::uint64_t m_long_exp_len_ram_use = m_long_exp_len.ram_use();
      std::uint64_t m_long_exp_hashes_ram_use = 0;
      std::uint64_t m_hashes_ram_use = 0;
      std::uint64_t m_kr_hash_cache_ram_use = 0;;
      if (m_enable_kr_hashing) {
        m_long_exp_hashes_ram_use = m_long_exp_hashes.ram_use();
        m_hashes_ram_use = m_hashes.ram_use();
        m_kr_hash_cache_ram_use = m_kr_hash_cache->ram_use();
      }
      const std::uint64_t total =
        m_roots_vec_ram_use +
        m_nonterminals_ram_use + 
        m_long_exp_hashes_ram_use +
        m_long_exp_len_ram_use +
        m_hashes_ram_use +
        m_kr_hash_cache_ram_use;
      return total;
    }

    //=========================================================================
    // Print statistics.
    //=========================================================================
    void print_stats() const {

      // Print RAM usage breakdown.
      const std::uint64_t m_roots_vec_ram_use = m_roots_vec.ram_use();
      const std::uint64_t m_nonterminals_ram_use = m_nonterminals.ram_use();
      const std::uint64_t m_long_exp_len_ram_use = m_long_exp_len.ram_use();
      std::uint64_t m_long_exp_hashes_ram_use = 0;
      std::uint64_t m_hashes_ram_use = 0;
      std::uint64_t m_kr_hash_cache_ram_use = 0;;
      if (m_enable_kr_hashing) {
        m_long_exp_hashes_ram_use = m_long_exp_hashes.ram_use();
        m_hashes_ram_use = m_hashes.ram_use();
        m_kr_hash_cache_ram_use = m_kr_hash_cache->ram_use();
      }
      const std::uint64_t total =
        m_roots_vec_ram_use +
        m_nonterminals_ram_use + 
        m_long_exp_hashes_ram_use +
        m_long_exp_len_ram_use +
        m_hashes_ram_use +
        m_kr_hash_cache_ram_use;
      fprintf(stderr, "RAM use:\n");
      fprintf(stderr, "  m_roots_vec: %.2LfMiB (%.2Lf%%)\n",
          (1.L * m_roots_vec_ram_use) / (1 << 20),
          (100.L * m_roots_vec_ram_use) / total);
      fprintf(stderr, "  m_nonterminals: %.2LfMiB (%.2Lf%%)\n",
          (1.L * m_nonterminals_ram_use) / (1 << 20),
          (100.L * m_nonterminals_ram_use) / total);
      fprintf(stderr, "  m_long_exp_len: %.2LfMiB (%.2Lf%%)\n",
          (1.L * m_long_exp_len_ram_use) / (1 << 20),
          (100.L * m_long_exp_len_ram_use) / total);
      fprintf(stderr, "  m_long_exp_hashes: %.2LfMiB (%.2Lf%%)\n",
          (1.L * m_long_exp_hashes_ram_use) / (1 << 20),
          (100.L * m_long_exp_hashes_ram_use) / total);
      fprintf(stderr, "  m_hashes: %.2LfMiB (%.2Lf%%)\n",
          (1.L * m_hashes_ram_use) / (1 << 20),
          (100.L * m_hashes_ram_use) / total);
      fprintf(stderr, "  m_kr_hash_cache: %.2LfMiB (%.2Lf%%)\n",
          (1.L * m_kr_hash_cache_ram_use) / (1 << 20),
          (100.L * m_kr_hash_cache_ram_use) / total);
      fprintf(stderr, "Total: %.2LfMiB (%.2Lf%%)\n",
          (1.L * total) / (1 << 20),
          (100.L * total) / total);

      // Show cache miss rate.
      if (m_enable_kr_hashing) {
        const std::uint64_t cache_misses = m_kr_hash_cache->get_cache_misses();
        const std::uint64_t query_counter = m_kr_hash_cache->get_query_counter();
        fprintf(stderr, "KR hash cache misses: %lu (%.2Lf%%)\n",
            cache_misses, (100.L * cache_misses) / query_counter);
      }
    }

  private:

    //=========================================================================
    // Find the leftmost nondeleted root >= key.
    //=========================================================================
    std::uint64_t roots_lower_bound(const std::uint64_t key) {
      std::uint64_t beg = 0;
      std::uint64_t end = roots_end();
      while (beg + 1 < end) {
        const std::uint64_t mid = (beg + end - 1) / 2;
        const std::uint64_t mid_end = m_roots_vec[mid].first;
        if (mid_end >= key)
          end = mid + 1;
        else beg = mid + 1;
      }

      // Skip deleted elements.
      while (m_roots_vec[beg].second ==
          std::numeric_limits<ptr_type>::max()) {
        ++beg;
        ++m_empty_step_counter;
      }

      // Return the result.
      return beg;
    }

    //=========================================================================
    // Find the first undeleted root.
    //=========================================================================
    inline std::uint64_t roots_begin() const {

      // We use the fact that there is a sentinel at the beginning.
      return 0;
    }

    //=========================================================================
    // Return the past-the-end position in the roots array.
    //=========================================================================
    inline std::uint64_t roots_end() const {
      return m_roots_vec.size();
    }


    //=========================================================================
    // Find the next undeleted root.
    //=========================================================================
    std::uint64_t roots_next(std::uint64_t pos) {
      ++pos;

      // Skip deleted elements.
      while (pos != roots_end() &&
          m_roots_vec[pos].second ==
          std::numeric_limits<ptr_type>::max()) {
        ++pos;
        ++m_empty_step_counter;
      }

      // Return the result.
      return pos;
    }

    //=========================================================================
    // Bring the elements in m_roots_vec together.
    //=========================================================================
    void roots_garbage_collector() {
      std::uint64_t filled = 0;
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i)) {
        if (i != filled)
          m_roots_vec[filled] = m_roots_vec[i];
        ++filled;
      }
      while (m_roots_vec.size() > filled)
        m_roots_vec.pop_back();
    }

  public:

    //=========================================================================
    // Merge roots enclosed in [begin..end).
    //=========================================================================
    void merge_enclosed_roots(
        const std::uint64_t begin,
        const std::uint64_t end) {

      // Compute the iterators of elements in m_roots_vec to merge.
      space_efficient_vector<triple_type> v;
      std::uint64_t range_beg = roots_lower_bound(begin);
      std::uint64_t prev_end = m_roots_vec[range_beg].first;
      range_beg = roots_next(range_beg);
      std::uint64_t range_end = range_beg;
      std::uint64_t newend = 0;
      while (range_end != roots_end() &&
          (std::uint64_t)m_roots_vec[range_end].first <= end) {
        const std::uint64_t cur_end = m_roots_vec[range_end].first;
        const std::uint64_t cur_exp_size = cur_end - prev_end;
        const std::uint64_t nonterm_p = m_roots_vec[range_end].second;
        std::uint64_t kr_hash = 0;
        if (m_enable_kr_hashing)
          kr_hash = get_kr_hash(nonterm_p);
        v.push_back(triple_type(nonterm_p, cur_exp_size, kr_hash));
        m_roots_vec[range_end].second = std::numeric_limits<ptr_type>::max();
        prev_end = cur_end;
        newend = range_end;
        range_end = roots_next(range_end);
      }

      // Merge the roots in the range.
      if (!v.empty()) {
        const std::uint64_t newroot_id = greedy_merge(v);
        m_roots_vec[newend].second = newroot_id;
      }
    }

    //=========================================================================
    // Run garbage collector if needed.
    //=========================================================================
    void check_gargage_collector() {
      if (m_empty_step_counter > m_roots_vec.size()) {
        roots_garbage_collector();
        m_empty_step_counter = 0;
      }
    }

    //=========================================================================
    // Get the sequence of nonterminals expanding to T[begin..end).
    //=========================================================================
    void decomposition(
        std::uint64_t begin,
        std::uint64_t end,
        space_efficient_vector<pair_type> &ret) {

      // Find leftmost root whose expansion overlaps/touches T[begin..end).
      std::uint64_t pos = roots_lower_bound(begin);

      // Proper substring or suffix of expansion of `it'.
      if (begin < (std::uint64_t)m_roots_vec[pos].first) {
        const std::uint64_t cur_end = m_roots_vec[pos].first;
        const std::uint64_t id = m_roots_vec[pos].second;
        const std::uint64_t exp_len = get_exp_len(id);
        const std::uint64_t it_exp_beg = cur_end - exp_len;
        const std::uint64_t local_beg = begin - it_exp_beg;
        const std::uint64_t local_end = std::min(cur_end, end) - it_exp_beg;
        const std::uint64_t local_size = local_end - local_beg;
        const nonterminal_type &nonterm = get_nonterminal(id);
        nonterm.decomposition(id, local_beg, local_end, ret, this);
        begin += local_size;
      }

      // Full expansions of nonterminals.
      std::uint64_t prev_end = m_roots_vec[pos].first;
      pos = roots_next(pos);
      while (begin < end && (std::uint64_t)m_roots_vec[pos].first <= end) {
        const std::uint64_t cur_end = m_roots_vec[pos].first;
        const std::uint64_t id = m_roots_vec[pos].second;
        const std::uint64_t exp_len = cur_end - prev_end;
        ret.push_back(pair_type(id, exp_len));
        begin = cur_end;
        prev_end = cur_end;
        pos = roots_next(pos);
      }

      // Proper suffix of expansion of `it'.
      if (begin < end) {
        const std::uint64_t cur_end = m_roots_vec[pos].first;
        const std::uint64_t id = m_roots_vec[pos].second;
        const std::uint64_t exp_len = cur_end - prev_end;
        const std::uint64_t it_exp_beg = cur_end - exp_len;
        const std::uint64_t local_end = end - it_exp_beg;
        const nonterminal_type &nonterm = get_nonterminal(id);
        space_efficient_vector<pair_type> dec;
        nonterm.decomposition(id, 0, local_end, dec, this);
        for (std::uint64_t i = 0; i < dec.size(); ++i)
          ret.push_back(dec[i]);
        begin = end;
      }
    }

    //=========================================================================
    // Return the sequence of nonterminals with the same expansion as seq.
    //=========================================================================
    void find_equivalent_seq(
      space_efficient_vector<pair_type> &seq) const {

      // Handle special case.
      if (!m_enable_kr_hashing || seq.empty())
        return;

      // Create the vector to hold the solution.
      space_efficient_vector<pair_type> ret;

      // Allocate the arrays used in the dynamic programming.
      const std::uint64_t length = seq.size();
      std::uint64_t *kr_hashes = utils::allocate_array<std::uint64_t>(length);
      std::uint64_t **dp = utils::allocate_array<std::uint64_t*>(length);
      std::uint64_t **dp_sol = utils::allocate_array<std::uint64_t*>(length);
      std::uint64_t **dp_explen = utils::allocate_array<std::uint64_t*>(length);
      ptr_type **dp_nonterm = utils::allocate_array<ptr_type*>(length);
      for (std::uint64_t i = 0; i < length; ++i) {
        dp[i] = utils::allocate_array<std::uint64_t>(length);
        dp_sol[i] = utils::allocate_array<std::uint64_t>(length);
        dp_explen[i] = utils::allocate_array<std::uint64_t>(length);
        dp_nonterm[i] = utils::allocate_array<ptr_type>(length);
      }

      // Fill in the array for len = 1.
      for (std::uint64_t i = 0; i < length; ++i) {
        const std::uint64_t id = seq[i].first;
        const std::uint64_t exp_len = seq[i].second;
        dp[i][i] = 1;
        dp_sol[i][i] = 1;
        dp_nonterm[i][i] = id;
        dp_explen[i][i] = exp_len;
        kr_hashes[i] = get_kr_hash(id);
      }

      // Solve for subarray of length > 1.
      for (std::uint64_t len = 2; len <= length; ++len) {
        for (std::uint64_t beg = 0; beg <= length - len; ++beg) {
          const std::uint64_t end = beg + len - 1;
          std::uint64_t exp_len = seq[beg].second;
          std::uint64_t h = kr_hashes[beg];

          // Initialize to solution to initial choice.
          dp[beg][end] = 1 + dp[beg + 1][end];
          dp_sol[beg][end] = 1;
          dp_nonterm[beg][end] = seq[beg].first;
          dp_explen[beg][end] = exp_len;

          // Try all other possible choices.
          for (std::uint64_t leftlen = 2; leftlen <= len; ++leftlen) {
            const std::uint64_t last = beg + leftlen - 1;
            const std::uint64_t right_hash = kr_hashes[last];
            const std::uint64_t right_len = seq[last].second;
            exp_len += right_len;
            h = karp_rabin_hashing::concat(h, right_hash, right_len);
            const ptr_type *nonterm_id_ptr = m_hashes.find(h);
            if (nonterm_id_ptr != NULL) {
              std::uint64_t sol_cost = 1;
              if (leftlen < len) sol_cost += dp[last + 1][end];
              if (sol_cost < dp[beg][end]) {
                dp[beg][end] = sol_cost;
                dp_sol[beg][end] = leftlen;
                dp_nonterm[beg][end] = *nonterm_id_ptr;
                dp_explen[beg][end] = exp_len;
              }
            }
          }
        }
      }

      // Restore the optimal solution.
      std::uint64_t prefix_length = 0;
      while (prefix_length < length) {
        const std::uint64_t id = dp_nonterm[prefix_length][length - 1];
        const std::uint64_t exp_len = dp_explen[prefix_length][length - 1];
        ret.push_back(pair_type(id, exp_len));
        prefix_length += dp_sol[prefix_length][length - 1];;
      }

      // Clean up.
      for (std::uint64_t len = length; len > 0; --len) {
        const std::uint64_t i = len - 1;
        utils::deallocate(dp_nonterm[i]);
        utils::deallocate(dp_explen[i]);
        utils::deallocate(dp_sol[i]);
        utils::deallocate(dp[i]);
      }
      utils::deallocate(dp_nonterm);
      utils::deallocate(dp_explen);
      utils::deallocate(dp_sol);
      utils::deallocate(dp);
      utils::deallocate(kr_hashes);

      // Store the result in seq.
      seq.set_empty();
      for (std::uint64_t i = 0; i < ret.size(); ++i)
        seq.push_back(ret[i]);
    }

    //=========================================================================
    // Return the fraction avoided merges.
    //=========================================================================
    long double get_avoided_merges() const {
      return (long double)m_avoided_merges / (long double)m_merge_count;
    }

    //=========================================================================
    // Write the grammar to a file.
    //=========================================================================
    void write_to_file(const std::string filename) {

      // Create the output writer.
      const std::uint64_t bufsize = (1 << 19);
      typedef async_stream_writer<text_offset_type> writer_type;
      writer_type *writer = new writer_type(filename, bufsize, 4);

      // Run the roots garbage collector.
      roots_garbage_collector();

      // Write its length, and then the roots sequence.
      writer->write((text_offset_type)(m_roots_vec.size() - 1));
      for (std::uint64_t i = 1; i < m_roots_vec.size(); ++i) {
        const std::uint64_t p = (std::uint64_t)m_roots_vec[i].second + 1;
        writer->write((text_offset_type)p);
      }

      // Write expansions of non-root nonterminals.
      // Note: some may be unused in the grammar!
      for (std::uint64_t i = 0; i < m_nonterminals.size(); ++i) {
        const nonterminal_type &nonterm = get_nonterminal(i);
        const std::uint64_t height = nonterm.get_height();
        if (height == 0) {
          const char_type c = nonterm.get_char();
          writer->write((text_offset_type)0);
          writer->write((text_offset_type)c);
        } else {
          const std::uint64_t left = (std::uint64_t)nonterm.get_left_p() + 1;
          const std::uint64_t right = (std::uint64_t)nonterm.get_right_p() + 1;
          writer->write((text_offset_type)left);
          writer->write((text_offset_type)right);
        }
      }

      // Clean up.
      delete writer;
    }

  private:

    //=========================================================================
    // Heap down routine.
    //=========================================================================
    void heap_down(
        std::uint64_t i,
        const space_efficient_vector<triple_type> &seq,
        text_offset_type * const heap,
        const std::uint64_t heap_size) const {
      ++i;
      std::uint64_t min_pos = i;
      std::uint64_t height = 0;
      const std::uint64_t val = heap[min_pos - 1];
      {
        const ptr_type nonterm_p = seq[val].first;
        const nonterminal_type &nonterm = get_nonterminal(nonterm_p);
        height = nonterm.get_height();
      }
      while (true) {
        std::uint64_t min_height = height;
        std::uint64_t min_val = val;
        if ((i << 1) <= heap_size) {
          const std::uint64_t left_val = heap[(i << 1) - 1];
          const ptr_type left_p = seq[left_val].first;
          const nonterminal_type &left = get_nonterminal(left_p);
          const std::uint64_t left_height = left.get_height();
          if (left_height < min_height ||
              (left_height == min_height && left_val < min_val)) {
            min_pos = (i << 1);
            min_height = left_height;
            min_val = left_val;
          }
        }
        if ((i << 1) + 1 <= heap_size) {
          const std::uint64_t right_val = heap[i << 1];
          const ptr_type right_p = seq[right_val].first;
          const nonterminal_type &right = get_nonterminal(right_p);
          const std::uint64_t right_height = right.get_height();
          if (right_height < min_height ||
              (right_height == min_height && right_val < min_val)) {
            min_pos = (i << 1) + 1;
            min_height = right_height;
            min_val = right_val;
          }
        }
        if (min_pos != i) {
          std::swap(heap[i - 1], heap[min_pos - 1]);
          i = min_pos;
        } else break;
      }
    }

    //=========================================================================
    // Extract min routine.
    //=========================================================================
    std::uint64_t extract_min(
        const space_efficient_vector<triple_type> &seq,
        text_offset_type * const heap,
        std::uint64_t &heap_size) const {
      std::uint64_t ret = heap[0];
      heap[0] = heap[--heap_size];
      if (heap_size != 0)
        heap_down(0, seq, heap, heap_size);
      return ret;
    }

    //=========================================================================
    // Make heap routine.
    //=========================================================================
    void make_heap(
        const space_efficient_vector<triple_type> &seq,
        text_offset_type * const heap,
        const std::uint64_t heap_size) const {
      for (std::uint64_t i = heap_size / 2; i > 0; --i)
        heap_down(i - 1, seq, heap, heap_size);
    }

    //=========================================================================
    // Merge greedily (shortest first) sequence of nonterminals.
    // Uses binary heap to achieve O(m log m) time.
    //=========================================================================
    ptr_type greedy_merge(space_efficient_vector<triple_type> &seq) {

      // Create the priority queue.
      const std::uint64_t num = seq.size();
      text_offset_type * const heap =
        utils::allocate_array<text_offset_type>(num);
      std::uint64_t heap_size = 0;

      // Allocate the arrays used to doubly-link remaining nonterminals.
      text_offset_type * const next =
        utils::allocate_array<text_offset_type>(num + 1);
      text_offset_type * const prev =
        utils::allocate_array<text_offset_type>(num + 1);
      std::uint8_t * const deleted =
        utils::allocate_array<std::uint8_t>(num);

      // Set initial linking and insert elements into heap.
      const std::uint64_t sentinel = num;
      for (std::uint64_t i = 0; i < num; ++i) {
        next[i] = i + 1;
        prev[i + 1] = i;
        heap[heap_size++] = i;
        deleted[i] = false;
      }
      next[sentinel] = 0;
      prev[0] = sentinel;
      make_heap(seq, heap, heap_size);

      // The main algorithm.
      ptr_type ret = std::numeric_limits<ptr_type>::max();
      while (true) {
        const std::uint64_t min_elem = heap[0];

        // If the element was already deleted, skip it.
        if (deleted[min_elem]) {
          (void) extract_min(seq, heap, heap_size);
          continue;
        }

        // Compute the height of the previous and next elements.
        std::uint64_t prev_height = 0;
        std::uint64_t next_height = 0;
        if ((std::uint64_t)prev[min_elem] != sentinel) {
          const std::uint64_t idx = prev[min_elem];
          const ptr_type prev_nonterm_p = seq[idx].first;
          const nonterminal_type &prev_nonterm =
            get_nonterminal(prev_nonterm_p);
          prev_height = prev_nonterm.get_height();
        }
        if ((std::uint64_t)next[min_elem] != sentinel) {
          const std::uint64_t idx = next[min_elem];
          const ptr_type next_nonterm_p = seq[idx].first;
          const nonterminal_type &next_nonterm =
            get_nonterminal(next_nonterm_p);
          next_height = next_nonterm.get_height();
        }

        // Merge min_elem with one of its
        // neighbors (whichever is shorter).
        if ((std::uint64_t)prev[min_elem] == sentinel &&
            (std::uint64_t)next[min_elem] == sentinel) {

          // We have reached the only nonterminal.
          // Store it as output and exit the loop.
          ret = seq[min_elem].first;
          break;
        } else if ((std::uint64_t)prev[min_elem] == sentinel ||
            ((std::uint64_t)next[min_elem] != sentinel &&
             next_height <= prev_height)) {

          // Only right neighbor exists, or both exist
          // and the right one is not taller than the left
          // one. End result: merge with the right neighbor.
          const std::uint64_t right_elem = next[min_elem];
          const ptr_type left_p = seq[min_elem].first;
          const ptr_type right_p = seq[right_elem].first;
          const std::uint64_t left_len = seq[min_elem].second;
          const std::uint64_t right_len = seq[right_elem].second;
          const std::uint64_t merged_len = left_len + right_len;
          std::uint64_t h = 0;
          ptr_type id_merged = 0;
          ++m_merge_count;
          if (m_enable_kr_hashing) {
            const std::uint64_t left_hash = seq[min_elem].third;
            const std::uint64_t right_hash = seq[right_elem].third;
            h = karp_rabin_hashing::concat(left_hash, right_hash, right_len);
            const ptr_type * const hash_ret = m_hashes.find(h);
            if (hash_ret != NULL) {
              id_merged = *hash_ret;
              ++m_avoided_merges;
            } else id_merged = add_concat_nonterminal(left_p, right_p);
          } else id_merged = add_concat_nonterminal(left_p, right_p);
          seq[min_elem].first = id_merged;
          seq[min_elem].second = merged_len;
          seq[min_elem].third = h;
          deleted[right_elem] = true;
          next[min_elem] = next[right_elem];
          prev[next[min_elem]] = min_elem;
          heap_down(0, seq, heap, heap_size);
        } else {

          // Only left neighbor exists, or both exists
          // and the left one is not taller than the
          // right one. End result: merge with left neighbor.
          const std::uint64_t left_elem = prev[min_elem];
          const ptr_type left_p = seq[left_elem].first;
          const ptr_type right_p = seq[min_elem].first;
          const std::uint64_t left_len = seq[left_elem].second;
          const std::uint64_t right_len = seq[min_elem].second;
          const std::uint64_t merged_len = left_len + right_len;
          std::uint64_t h = 0;
          ptr_type id_merged = 0;
          ++m_merge_count;
          if (m_enable_kr_hashing) {
            const std::uint64_t left_hash = seq[left_elem].third;
            const std::uint64_t right_hash = seq[min_elem].third;
            h = karp_rabin_hashing::concat(left_hash, right_hash, right_len);
            const ptr_type * const hash_ret = m_hashes.find(h);
            if (hash_ret != NULL) {
              id_merged = *hash_ret;
              ++m_avoided_merges;
            } else id_merged = add_concat_nonterminal(left_p, right_p);
          } else id_merged = add_concat_nonterminal(left_p, right_p);
          seq[min_elem].first = id_merged;
          seq[min_elem].second = merged_len;
          seq[min_elem].third = h;
          deleted[left_elem] = true;
          prev[min_elem] = prev[left_elem];
          next[prev[min_elem]] = min_elem;
          heap_down(0, seq, heap, heap_size);
        }
      }

      // Clean up.
      utils::deallocate(deleted);
      utils::deallocate(prev);
      utils::deallocate(next);
      utils::deallocate(heap);

      // Return the result.
      return ret;
    }
};

//=============================================================================
// Default constructor.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal()
  : m_height(0),
    m_truncated_exp_len(1),
    m_left_p(std::numeric_limits<text_offset_type>::max()),
    m_right_p(std::numeric_limits<text_offset_type>::max()) {}

//=============================================================================
// Constructor for a nonterminal expanding to a single symbol.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(const char_type c)
  : m_height(0),
    m_truncated_exp_len(1),
    m_left_p((text_offset_type)((std::uint64_t)c)),
    m_right_p(std::numeric_limits<text_offset_type>::max()) {}

//=============================================================================
// Constructor for non-single-symbol nonterminal.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(
      const std::uint8_t height,
      const std::uint8_t truncated_exp_len,
      const ptr_type left_p,
      const ptr_type right_p)
  : m_height(height),
    m_truncated_exp_len(truncated_exp_len),
    m_left_p(left_p),
    m_right_p(right_p) {}

//=============================================================================
// Get nonterminal height.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>::get_height() const {
  return (std::uint64_t)m_height;
}

//=============================================================================
// Get nonterminal expansion length.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>
::get_truncated_exp_len() const {
  return (std::uint64_t)m_truncated_exp_len;
}

//=============================================================================
// Return the char stored in a nonterminal.
//=============================================================================
template<typename char_type, typename text_offset_type>
char_type nonterminal<char_type, text_offset_type>::get_char() const {
  return (char_type)((std::uint64_t)m_left_p);
}

//=============================================================================
// Get nonterminal left ptr.
//=============================================================================
template<typename char_type, typename text_offset_type>
text_offset_type nonterminal<char_type, text_offset_type>::get_left_p() const {
  return m_left_p;
}

//=============================================================================
// Get nonterminal left ptr.
//=============================================================================
template<typename char_type, typename text_offset_type>
text_offset_type nonterminal<char_type, text_offset_type>::get_right_p() const {
  return m_right_p;
}

//=============================================================================
// Print expansion of a given nonterminal.
//=============================================================================
template<typename char_type, typename text_offset_type>
void nonterminal<char_type, text_offset_type>::print_expansion(
    const ptr_type x_p,
    const grammar_type * const g) const {
  const nonterminal_type &x = g->get_nonterminal(x_p);
  const std::uint64_t height = x.get_height();
  if (height == 0) {
    const char_type my_char = x.get_char();
    fprintf(stderr, "%c", (char)my_char);
  } else {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = g->get_nonterminal(x_left_p);
    const nonterminal_type &x_right = g->get_nonterminal(x_right_p);
    x_left.print_expansion(x_left_p, g);
    x_right.print_expansion(x_right_p, g);
  }
}

//=============================================================================
// Write the expansion into the given array.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>::write_expansion(
    const ptr_type x_p,
    char_type * const text,
    const grammar_type * const g) const {
  const nonterminal_type &x = g->get_nonterminal(x_p);
  const std::uint64_t x_height = x.get_height();
  if (x_height == 0) {
    const char_type my_char = x.get_char();
    text[0] = my_char;
    return 1;
  } else {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = g->get_nonterminal(x_left_p);
    const nonterminal_type &x_right = g->get_nonterminal(x_right_p);
    const std::uint64_t x_left_exp_len =
      x_left.write_expansion(x_left_p, text, g);
    const std::uint64_t x_right_exp_len =
      x_right.write_expansion(x_right_p, text + x_left_exp_len, g);
    return x_left_exp_len + x_right_exp_len;
  }
}

//=============================================================================
// Compare the expansion of the nonterminal to the given text.
//=============================================================================
template<typename char_type, typename text_offset_type>
bool nonterminal<char_type, text_offset_type>::compare_expansion_to_text(
    const ptr_type x_p,
    const char_type * const text,
    const grammar_type * const g) const {
  const nonterminal_type &x = g->get_nonterminal(x_p);
  const std::uint64_t x_height = x.get_height();
  if (x_height == 0) {
    const char_type my_char = x.get_char();
    return (text[0] == my_char);
  } else {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = g->get_nonterminal(x_left_p);
    const nonterminal_type &x_right = g->get_nonterminal(x_right_p);
    const std::uint64_t x_left_exp_len = g->get_exp_len(x_left_p);
    if (!x_left.compare_expansion_to_text(x_left_p, text, g))
      return false;
    if (!x_right.compare_expansion_to_text(x_right_p, text + x_left_exp_len, g))
      return false;
    return true;
  }
}

//=============================================================================
// Test the AVL property of a subtree.
//=============================================================================
template<typename char_type, typename text_offset_type>
bool nonterminal<char_type, text_offset_type>::test_avl_property(
    const ptr_type x_p,
    const grammar_type * const g) const {
  const nonterminal_type &x = g->get_nonterminal(x_p);
  const std::uint64_t x_height = x.get_height();

  // Handle special case.
  if (x_height == 0)
    return true;

  // Recursively check AVL property for children.
  const ptr_type x_left_p = x.get_left_p();
  const ptr_type x_right_p = x.get_right_p();
  const nonterminal_type &x_left = g->get_nonterminal(x_left_p);
  const nonterminal_type &x_right = g->get_nonterminal(x_right_p);
  if (!x_left.test_avl_property(x_left_p, g) ||
      !x_right.test_avl_property(x_right_p, g))
    return false;

  // Check the AVL property for the node.
  const std::uint64_t x_left_height = x_left.get_height();
  const std::uint64_t x_right_height = x_right.get_height();
  if ((x_right_height > x_left_height && x_right_height - x_left_height > 1) ||
      (x_right_height < x_left_height && x_left_height - x_right_height > 1))
    return false;

  // The subtree satisfies the AVL property.
  return true;
}

//=============================================================================
// Collect Mersenne Karp-Rabin hashes of all nonterminals.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>
::collect_mersenne_karp_rabin_hashes(
    const ptr_type x_p,
    std::vector<std::uint64_t> &hashes,
    const grammar_type * const g) const {
  const nonterminal_type &x = g->get_nonterminal(x_p);
  const std::uint64_t x_height = x.get_height();
  if (x_height == 0) {
    const char_type my_char = x.get_char();
    const std::uint64_t h = karp_rabin_hashing::hash_char(my_char);
    hashes.push_back(h);
    return h;
  } else {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = g->get_nonterminal(x_left_p);
    const nonterminal_type &x_right = g->get_nonterminal(x_right_p);
    const std::uint64_t x_right_len = g->get_exp_len(x_right_p);
    const std::uint64_t x_left_hash =
      x_left.collect_mersenne_karp_rabin_hashes(x_left_p, hashes, g);
    const std::uint64_t x_right_hash =
      x_right.collect_mersenne_karp_rabin_hashes(x_right_p, hashes, g);
    const std::uint64_t h = karp_rabin_hashing::concat(
        x_left_hash, x_right_hash, x_right_len);
    hashes.push_back(h);
    return h;
  }
}

//=============================================================================
// Collect reachable nonterminals.
//=============================================================================
template<typename char_type, typename text_offset_type>
void nonterminal<char_type, text_offset_type>::collect_nonterminal_pointers(
    const ptr_type x_p,
    std::vector<ptr_type> &pointers,
    const grammar_type * const g) const {
  const nonterminal_type &x = g->get_nonterminal(x_p);
  const std::uint64_t x_height = x.get_height();
  pointers.push_back(x_p);
  if (x_height > 0) {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = g->get_nonterminal(x_left_p);
    const nonterminal_type &x_right = g->get_nonterminal(x_right_p);
    x_left.collect_nonterminal_pointers(x_left_p, pointers, g);
    x_right.collect_nonterminal_pointers(x_right_p, pointers, g);
  }
}

//=============================================================================
// Collect Mersenne Karp-Rabin hashes of all nonterminals.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>
::collect_mersenne_karp_rabin_hashes_2(
    const ptr_type x_p,
    hash_table<ptr_type, std::uint64_t> &hashes,
    const grammar_type * const g) const {
  const nonterminal_type &x = g->get_nonterminal(x_p);
  const std::uint64_t x_height = x.get_height();
  if (x_height == 0) {
    const char_type my_char = x.get_char();
    const std::uint64_t h = karp_rabin_hashing::hash_char(my_char);
    hashes.insert(x_p, h);
    return h;
  } else {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = g->get_nonterminal(x_left_p);
    const nonterminal_type &x_right = g->get_nonterminal(x_right_p);
    const std::uint64_t x_right_len = g->get_exp_len(x_right_p);
    const std::uint64_t x_left_hash =
      x_left.collect_mersenne_karp_rabin_hashes_2(x_left_p, hashes, g);
    const std::uint64_t x_right_hash =
      x_right.collect_mersenne_karp_rabin_hashes_2(x_right_p, hashes, g);
    const std::uint64_t h = karp_rabin_hashing::concat(
        x_left_hash, x_right_hash, x_right_len);
    hashes.insert(x_p, h);
    return h;
  }
}

//=============================================================================
// Compute the number of nonterminals in the pruned grammar.
//=============================================================================
template<typename char_type, typename text_offset_type>
void nonterminal<char_type, text_offset_type>
::count_nonterminals_in_pruned_grammar(
    const ptr_type x_p,
    hash_table<ptr_type, std::uint64_t> &hashes,
    hash_table<std::uint64_t, bool> &seen_hashes,
    std::uint64_t &current_count,
    const grammar_type * const g) const {
  const std::uint64_t * const h = hashes.find(x_p);
  if (seen_hashes.find(*h) == NULL) {
    seen_hashes.insert(*h, true);
    ++current_count;
    const nonterminal_type &x = g->get_nonterminal(x_p);
    const std::uint64_t x_height = x.get_height();
    if (x_height != 0) {
      const ptr_type x_left_p = x.get_left_p();
      const ptr_type x_right_p = x.get_right_p();
      const nonterminal_type &x_left = g->get_nonterminal(x_left_p);
      const nonterminal_type &x_right = g->get_nonterminal(x_right_p);
      x_left.count_nonterminals_in_pruned_grammar(
          x_left_p, hashes, seen_hashes, current_count, g);
      x_right.count_nonterminals_in_pruned_grammar(
          x_right_p, hashes, seen_hashes, current_count, g);
    }
  }
}

//=============================================================================
// Assuming S is the expansions of the nonterminal, return the
// sequence of nonterminals expanding to S[begin..end).
//=============================================================================
template<typename char_type, typename text_offset_type>
void nonterminal<char_type, text_offset_type>::decomposition(
    ptr_type x_p,
    const std::uint64_t begin,
    const std::uint64_t end,
    space_efficient_vector<pair_type> &ret,
    const grammar_type * const g) const {

  // Handle boundary case.
  if (begin == end)
    return;

  // Compute height and expansion length for x.
  std::uint64_t x_height = 0;
  std::uint64_t x_exp_len = 0;
  {
    const nonterminal_type &x = g->get_nonterminal(x_p);
    x_height = x.get_height();
    x_exp_len = g->get_exp_len(x_p);
  }

  // Find the deepest nonterminal in the parse tree
  // containing the range [begin..end).
  std::uint64_t cur_range_beg = 0;
  std::uint64_t cur_range_end = x_exp_len;
  while (x_height > 0) {
    const nonterminal_type &x = g->get_nonterminal(x_p);
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = g->get_nonterminal(x_left_p);
    const nonterminal_type &x_right = g->get_nonterminal(x_right_p);
    const std::uint64_t x_left_exp_len = g->get_exp_len(x_left_p);
    const std::uint64_t cur_range_mid = cur_range_beg + x_left_exp_len;
    if (end <= cur_range_mid) {
      cur_range_end = cur_range_mid;
      x_p = x_left_p;
      x_height = x_left.get_height();
    } else if (begin >= cur_range_mid) {
      cur_range_beg = cur_range_mid;
      x_p = x_right_p;
      x_height = x_right.get_height();
    } else break;
  }

  // Check if the range of x is exactly [begin..end).
  if (cur_range_beg == begin && cur_range_end == end) {

    // If yes, return x as the answer.
    const std::uint64_t exp_len = end - begin;
    ret.push_back(pair_type(x_p, exp_len));
  } else {

    // Otherwise, we perform two traversals in the tree.
    {
      const nonterminal_type &x = g->get_nonterminal(x_p);
      const ptr_type x_left_p = x.get_left_p();
      const std::uint64_t x_left_exp_len = g->get_exp_len(x_left_p);
      const std::uint64_t left_range_end = cur_range_beg + x_left_exp_len;
      std::uint64_t suffix_length = left_range_end - begin;
      ptr_type y_p = x_left_p;
      while (suffix_length > 0) {
        const nonterminal_type &y = g->get_nonterminal(y_p);
        const std::uint64_t y_exp_len = g->get_exp_len(y_p);
        const ptr_type y_left_p = y.get_left_p();
        const ptr_type y_right_p = y.get_right_p();
        if (y_exp_len == suffix_length) {
          ret.push_back(pair_type(y_p, y_exp_len));
          suffix_length -= y_exp_len;
        } else {
          const std::uint64_t y_right_exp_len = g->get_exp_len(y_right_p);
          if (suffix_length > y_right_exp_len) {
            ret.push_back(pair_type(y_right_p, y_right_exp_len));
            suffix_length -= y_right_exp_len;
            y_p = y_left_p;
          } else y_p = y_right_p;
        }
      }
    }

    // Reverse the first sequence of nonterminals
    // collected during the left downward traversal.
    ret.reverse();

    // Perform the analogous operation for the right side.
    {
      const nonterminal_type &x = g->get_nonterminal(x_p);
      const ptr_type x_left_p = x.get_left_p();
      const ptr_type x_right_p = x.get_right_p();
      const std::uint64_t x_left_exp_len = g->get_exp_len(x_left_p);
      const std::uint64_t right_range_beg = cur_range_beg + x_left_exp_len;
      std::uint64_t prefix_length = end - right_range_beg;
      ptr_type y_p = x_right_p;
      while (prefix_length > 0) {
        const nonterminal_type &y = g->get_nonterminal(y_p);
        const std::uint64_t y_exp_len = g->get_exp_len(y_p);
        const ptr_type y_left_p = y.get_left_p();
        const ptr_type y_right_p = y.get_right_p();
        if (y_exp_len == prefix_length) {
          ret.push_back(pair_type(y_p, y_exp_len));
          prefix_length -= y_exp_len;
        } else {
          const std::uint64_t y_left_exp_len = g->get_exp_len(y_left_p);
          if (prefix_length > y_left_exp_len) {
            ret.push_back(pair_type(y_left_p, y_left_exp_len));
            prefix_length -= y_left_exp_len;
            y_p = y_right_p;
          } else y_p = y_left_p;
        }
      }
    }
  }
}

#endif  // __LAZY_AVL_GRAMMAR_HPP_INCLUDED
