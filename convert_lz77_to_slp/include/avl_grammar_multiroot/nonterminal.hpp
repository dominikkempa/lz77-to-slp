#ifndef __NONTERMINAL_HPP_INCLUDED
#define __NONTERMINAL_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <map>
#include <algorithm>

#include "../utils/hash_table.hpp"
#include "../utils/karp_rabin_hashing.hpp"


//=============================================================================
// A class used to represent the nonterminal.
//=============================================================================
template<
  typename char_type,
  typename text_offset_type>
struct nonterminal {
  private:

    //=========================================================================
    // Declare typedefs.
    //=========================================================================
    typedef nonterminal<char_type, text_offset_type> nonterminal_type;

  public:

    //=========================================================================
    // Class members.
    //=========================================================================
    const char_type m_char;
    const std::uint8_t m_height;
    const std::uint64_t m_exp_len;
    const std::uint64_t m_kr_hash;
    const text_offset_type m_left;
    const text_offset_type m_right;

    //=========================================================================
    // Default constructor.
    //=========================================================================
    nonterminal() :
      m_char((char_type)0),
      m_height(0),
      m_exp_len(1),
      m_kr_hash(0),
      m_left(std::numeric_limits<text_offset_type>::max()),
      m_right(std::numeric_limits<text_offset_type>::max()) {}

    //=========================================================================
    // Constructor for a nonterminal expanding to a single symbol.
    //=========================================================================
    nonterminal(const char_type c) :
      m_char(c),
      m_height(0),
      m_exp_len(1),
      m_kr_hash(karp_rabin_hashing::hash_char(c)),
      m_left(std::numeric_limits<text_offset_type>::max()),
      m_right(std::numeric_limits<text_offset_type>::max()) {}

    //=========================================================================
    // Constructor for non-single-symbol nonterminal.
    //=========================================================================
    nonterminal(
        const std::uint64_t left,
        const std::uint64_t right,
        const std::vector<nonterminal_type> &nonterminals) :
          m_char((char_type)0),
          m_height(std::max(
                nonterminals[(std::uint64_t)left].m_height,
                nonterminals[(std::uint64_t)right].m_height) + 1),
          m_exp_len(
              nonterminals[(std::uint64_t)left].m_exp_len +
              nonterminals[(std::uint64_t)right].m_exp_len),
          m_kr_hash(
              karp_rabin_hashing::concat(
                nonterminals[(std::uint64_t)left].m_kr_hash,
                nonterminals[(std::uint64_t)right].m_kr_hash,
                nonterminals[(std::uint64_t)right].m_exp_len)),
          m_left(left),
          m_right(right) {}

    //=========================================================================
    // Print the string encoded by the grammar.
    //=========================================================================
    void print_expansion(
        const std::uint64_t id,
        const std::vector<nonterminal_type> &nonterminals) const {
      const std::uint64_t height = nonterminals[id].m_height;
      if (height == 0) {
        const char_type my_char = nonterminals[id].m_char;
        fprintf(stderr, "%c", (char)my_char);
      } else {
        const std::uint64_t left = nonterminals[id].m_left;
        const std::uint64_t right = nonterminals[id].m_right;
        nonterminals[left].print_expansion(left, nonterminals);
        nonterminals[right].print_expansion(right, nonterminals);
      }
    }

    //=========================================================================
    // Write the expansion into the given array.
    //=========================================================================
    void write_expansion(
        const std::uint64_t id,
        char_type * const text,
        const std::vector<nonterminal_type> &nonterminals) const {
      const std::uint64_t height = nonterminals[id].m_height;
      if (height == 0) {
        const char_type my_char = nonterminals[id].m_char;
        text[0] = my_char;
      } else {
        const std::uint64_t left = nonterminals[id].m_left;
        const std::uint64_t right = nonterminals[id].m_right;
        const std::uint64_t left_exp_len = nonterminals[left].m_exp_len;
        nonterminals[left].write_expansion(left, text, nonterminals);
        nonterminals[right].write_expansion(right, text + left_exp_len, nonterminals);
      }
    }

    //=========================================================================
    // Test the AVL propert of a subtree.
    //=========================================================================
    bool test_avl_property(
        const std::uint64_t id,
        const std::vector<nonterminal_type> &nonterminals) const {
      const std::uint64_t height = nonterminals[id].m_height;
      if (height == 0)
        return true;

      const std::uint64_t left = nonterminals[id].m_left;
      const std::uint64_t right = nonterminals[id].m_right;
      if (!nonterminals[left].test_avl_property(left, nonterminals) ||
          !nonterminals[right].test_avl_property(right, nonterminals))
        return false;

      const std::uint64_t left_height = nonterminals[left].m_height;
      const std::uint64_t right_height = nonterminals[right].m_height;
      if ((right_height > left_height && right_height - left_height > 1) ||
          (right_height < left_height && left_height - right_height > 1))
        return false;

      return true;
    }

    //=========================================================================
    // Collect Mersenne Karp-Rabin hashes of all nonterminals.
    //=========================================================================
    std::uint64_t collect_mersenne_karp_rabin_hashes(
        const std::uint64_t id,
        std::vector<std::uint64_t> &hashes,
        const std::vector<nonterminal_type> &nonterminals) const {
      const std::uint64_t height = nonterminals[id].m_height;
      if (height == 0) {
        const char_type my_char = nonterminals[id].m_char;
        const std::uint64_t h = karp_rabin_hashing::hash_char(my_char);
        hashes.push_back(h);
        return h;
      } else {
        const std::uint64_t left = nonterminals[id].m_left;
        const std::uint64_t right = nonterminals[id].m_right;
        const std::uint64_t left_hash =
          nonterminals[left].collect_mersenne_karp_rabin_hashes(
              left, hashes, nonterminals);
        const std::uint64_t right_hash =
          nonterminals[right].collect_mersenne_karp_rabin_hashes(
              right, hashes, nonterminals);
        const std::uint64_t right_len = nonterminals[right].m_exp_len;
        const std::uint64_t h = karp_rabin_hashing::concat(
            left_hash, right_hash, right_len);
        hashes.push_back(h);
        return h;
      }
    }

    //=========================================================================
    // Collect reachable nonterminals.
    //=========================================================================
    void collect_nonterminal_pointers(
        const std::uint64_t id,
        std::vector<text_offset_type> &pointers,
        const std::vector<nonterminal_type> &nonterminals) const {
      const std::uint64_t height = nonterminals[id].m_height;
      pointers.push_back(id);
      if (height > 0) {
        const std::uint64_t left = nonterminals[id].m_left;
        const std::uint64_t right = nonterminals[id].m_right;
        nonterminals[left].collect_nonterminal_pointers(
            left, pointers, nonterminals);
        nonterminals[right].collect_nonterminal_pointers(
            right, pointers, nonterminals);
      }
    }

    //=========================================================================
    // Collect Mersenne Karp-Rabin hashes of all nonterminals.
    //=========================================================================
    std::uint64_t collect_mersenne_karp_rabin_hashes_2(
        const std::uint64_t id,
        hash_table<text_offset_type, std::uint64_t> &hashes,
        const std::vector<nonterminal_type> &nonterminals) const {
      const std::uint64_t height = nonterminals[id].m_height;
      if (height == 0) {
        const char_type my_char = nonterminals[id].m_char;
        const std::uint64_t h = karp_rabin_hashing::hash_char(my_char);
        hashes.insert(id, h);
        return h;
      } else {
        const std::uint64_t left = nonterminals[id].m_left;
        const std::uint64_t right = nonterminals[id].m_right;
        const std::uint64_t left_hash =
          nonterminals[left].collect_mersenne_karp_rabin_hashes_2(
              left, hashes, nonterminals);
        const std::uint64_t right_hash =
          nonterminals[right].collect_mersenne_karp_rabin_hashes_2(
              right, hashes, nonterminals);
        const std::uint64_t right_len = nonterminals[right].m_exp_len;
        const std::uint64_t h = karp_rabin_hashing::concat(
            left_hash, right_hash, right_len);
        hashes.insert(id, h);
        return h;
      }
    }

    //=========================================================================
    // Compute the number of nonterminals in the pruned grammar.
    //=========================================================================
    void count_nonterminals_in_pruned_grammar(
        const std::uint64_t id,
        hash_table<text_offset_type, std::uint64_t> &hashes,
        hash_table<std::uint64_t, bool> &seen_hashes,
        std::uint64_t &current_count,
        const std::vector<nonterminal_type> &nonterminals) const {
      const std::uint64_t * const h = hashes.find((text_offset_type)id);
      if (seen_hashes.find(*h) == NULL) {
        seen_hashes.insert(*h, true);
        ++current_count;
        const std::uint64_t height = nonterminals[id].m_height;
        if (height != 0) {
          const std::uint64_t left = nonterminals[id].m_left;
          const std::uint64_t right = nonterminals[id].m_right;
          nonterminals[left].count_nonterminals_in_pruned_grammar(
              left, hashes, seen_hashes, current_count, nonterminals);
          nonterminals[right].count_nonterminals_in_pruned_grammar(
              right, hashes, seen_hashes, current_count, nonterminals);
        }
      }
    }

    //=========================================================================
    // Assuming S is the expansions of the nontermnal, return the
    // sequence of nonterminals expanding to S[begin..end).
    //=========================================================================
    std::vector<text_offset_type> decomposition(
        const std::uint64_t id,
        const std::uint64_t begin,
        const std::uint64_t end,
        const std::vector<nonterminal_type> &nonterminals) const {

      // Declare the vector storing the result.
      std::vector<text_offset_type> ret;

      // Handle boundary case.
      if (begin == end)
        return ret;

      // Find the deepest nonterminal in the parse tree containing the range
      // [begin..end).
      std::uint64_t x = id;
      std::uint64_t x_height = nonterminals[x].m_height;
      std::uint64_t cur_range_beg = 0;
      std::uint64_t cur_range_end = nonterminals[x].m_exp_len;
      while (x_height > 0) {
        const std::uint64_t x_left = nonterminals[x].m_left;
        const std::uint64_t x_right = nonterminals[x].m_right;
        const std::uint64_t x_left_exp_len = nonterminals[x_left].m_exp_len;
        const std::uint64_t cur_range_mid = cur_range_beg + x_left_exp_len;
        if (end <= cur_range_mid) {
          cur_range_end = cur_range_mid;
          x = x_left;
          x_height = nonterminals[x].m_height;
        } else if (begin >= cur_range_mid) {
          cur_range_beg = cur_range_mid;
          x = x_right;
          x_height = nonterminals[x].m_height;
        } else break;
      }

      // Check if the range of x is exactly [begin..end).
      if (cur_range_beg == begin && cur_range_end == end) {

        // If yes, return x as the answer.
        ret.push_back((text_offset_type)x);
      } else {

        // Otherwise, we perform two traversals in the tree.
        {
          const std::uint64_t x_left = nonterminals[x].m_left;
          const std::uint64_t x_left_exp_len = nonterminals[x_left].m_exp_len;
          const std::uint64_t left_range_end = cur_range_beg + x_left_exp_len;
          std::uint64_t suffix_length = left_range_end - begin;
          std::uint64_t y = x_left;
          while (suffix_length > 0) {
            const std::uint64_t y_exp_len = nonterminals[y].m_exp_len;
            const std::uint64_t y_left = nonterminals[y].m_left;
            const std::uint64_t y_right = nonterminals[y].m_right;
            if (y_exp_len == suffix_length) {
              ret.push_back(y);
              suffix_length -= y_exp_len;
            } else {
              const std::uint64_t y_right_exp_len =
                nonterminals[y_right].m_exp_len;
              if (suffix_length > y_right_exp_len) {
                ret.push_back(y_right);
                suffix_length -= y_right_exp_len;
                y = y_left;
              } else y = y_right;
            }
          }
        }

        // Reverse the first sequence of nonterminals
        // collected during the left downward traversal.
        std::reverse(ret.begin(), ret.end());

        // Perform the analogous operation for the right side.
        {
          const std::uint64_t x_left = nonterminals[x].m_left;
          const std::uint64_t x_right = nonterminals[x].m_right;
          const std::uint64_t x_left_exp_len = nonterminals[x_left].m_exp_len;
          const std::uint64_t right_range_beg = cur_range_beg + x_left_exp_len;
          std::uint64_t prefix_length = end - right_range_beg;
          std::uint64_t y = x_right;
          while (prefix_length > 0) {
            const std::uint64_t y_exp_len = nonterminals[y].m_exp_len;
            const std::uint64_t y_left = nonterminals[y].m_left;
            const std::uint64_t y_right = nonterminals[y].m_right;
            if (y_exp_len == prefix_length) {
              ret.push_back(y);
              prefix_length -= y_exp_len;
            } else {
              const std::uint64_t y_left_exp_len =
                nonterminals[y_left].m_exp_len;
              if (prefix_length > y_left_exp_len) {
                ret.push_back(y_left);
                prefix_length -= y_left_exp_len;
                y = y_right;
              } else y = y_left;
            }
          }
        }
      }

      // Return the result.
      return ret;
    }
};

//=============================================================================
// Hash functions of the appropriate type.
// Used in the hash table used to prune the grammar.
//=============================================================================
template<>
std::uint64_t get_hash(const uint40 &x) {
  return (std::uint64_t)x * (std::uint64_t)29996224275833;
}

template<>
std::uint64_t get_hash(const std::uint64_t &x) {
  return (std::uint64_t)x * (std::uint64_t)4972548694736365;
}

template<
  typename char_type,
  typename text_offset_type>
std::uint64_t merge_hashes(
    const std::uint64_t left,
    const std::uint64_t right,
    const std::vector<nonterminal<char_type, text_offset_type> >&nonterminals) {
  const std::uint64_t left_hash = nonterminals[left].m_kr_hash;
  const std::uint64_t right_hash = nonterminals[right].m_kr_hash;
  const std::uint64_t right_len = nonterminals[right].m_exp_len;
  const std::uint64_t h = karp_rabin_hashing::concat(
      left_hash, right_hash, right_len);
  return h;
}

template<
  typename char_type,
  typename text_offset_type>
std::uint64_t append_hash(
    const std::uint64_t left_hash,
    const std::uint64_t right,
    const std::vector<nonterminal<char_type, text_offset_type> > &nonterminals) {
  const std::uint64_t right_hash = nonterminals[right].m_kr_hash;
  const std::uint64_t right_len = nonterminals[right].m_exp_len;
  const std::uint64_t h = karp_rabin_hashing::concat(
      left_hash, right_hash, right_len);
  return h;
}

#endif  // __NONTERMINAL_HPP_INCLUDED
