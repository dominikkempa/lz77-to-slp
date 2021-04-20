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


//=============================================================================
// Class used to present multiroot AVL grammar. Forward declaration.
//=============================================================================
template<typename char_type, typename text_offset_type>
struct avl_grammar_multiroot;

//=============================================================================
// Class used to represent nonterminal. Forward declaration.
//=============================================================================
template<typename char_type, typename text_offset_type>
struct nonterminal {
  private:
    typedef nonterminal<char_type, text_offset_type> nonterminal_type;
    typedef avl_grammar_multiroot<char_type, text_offset_type> grammar_type;

  public:
    const char_type m_char;
    const std::uint8_t m_height;
    const std::uint64_t m_exp_len;
    const std::uint64_t m_kr_hash;
    const text_offset_type m_left;
    const text_offset_type m_right;

    nonterminal();
    nonterminal(const char_type);
    nonterminal(const std::uint64_t, const std::uint64_t,
        const grammar_type * const g);
    void print_expansion( const std::uint64_t,
        const grammar_type * const) const;
    void write_expansion(const std::uint64_t, char_type * const,
        const grammar_type * const) const;
    bool test_avl_property(const std::uint64_t,
        const grammar_type * const) const;
    std::uint64_t collect_mersenne_karp_rabin_hashes(const std::uint64_t,
        std::vector<std::uint64_t> &, const grammar_type * const) const;
    void collect_nonterminal_pointers(const std::uint64_t,
        std::vector<text_offset_type> &, const grammar_type * const) const;
    std::uint64_t collect_mersenne_karp_rabin_hashes_2(const std::uint64_t,
        hash_table<text_offset_type, std::uint64_t> &,
        const grammar_type * const) const;
    void count_nonterminals_in_pruned_grammar(const std::uint64_t,
        hash_table<text_offset_type, std::uint64_t> &,
        hash_table<std::uint64_t, bool> &,
        std::uint64_t &, const grammar_type * const) const;
    std::vector<text_offset_type> decomposition(
        const std::uint64_t, const std::uint64_t, const std::uint64_t,
        const grammar_type * const g) const;
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

//=============================================================================
// Implementation of the avl_grammar_multiroot class.
//=============================================================================
template<
  typename char_type,
  typename text_offset_type>
struct avl_grammar_multiroot {

  //===========================================================================
  // Declare typedefs
  //===========================================================================
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  typedef typename std::map<std::uint64_t, text_offset_type> map_type;
  typedef typename map_type::const_iterator const_iter_type;
  typedef typename map_type::iterator iter_type;

  private:

    //=========================================================================
    // Class members.
    //=========================================================================
    map_type m_roots;
    std::vector<nonterminal_type> m_nonterminals;
    hash_table<std::uint64_t, text_offset_type> m_hashes;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    avl_grammar_multiroot() {
      m_roots.insert(
          std::make_pair(0, std::numeric_limits<text_offset_type>::max()));
    }

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~avl_grammar_multiroot() {}

    //=========================================================================
    // Print the string encoded by the grammar.
    //=========================================================================
    void print_expansion() const {
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it) {
        std::uint64_t preflen = it->first;
        std::uint64_t id = it->second;
        if (preflen != 0)
          m_nonterminals[id].print_expansion(id, this);
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
    std::uint64_t number_of_roots() const {
      return m_roots.size() - 1;
    }

    //=========================================================================
    // Add a root.
    //=========================================================================
    void add_root(
        const std::uint64_t pos,
        const text_offset_type id) {
      m_roots[pos] = id;
    }

    //=========================================================================
    // Gives access to a given nonterminal.
    //=========================================================================
    const nonterminal_type& get_nonterminal(const std::uint64_t id) const {
      return m_nonterminals[id];
    }

    //=========================================================================
    // Return the height of a given nonterminal.
    //=========================================================================
    std::uint64_t get_height(const std::uint64_t id) const {
      return m_nonterminals[id].m_height;
    }

    //=========================================================================
    // Return the char associated with a given nonterminal.
    //=========================================================================
    char_type get_char(const std::uint64_t id) const {
      return m_nonterminals[id].m_char;
    }

    //=========================================================================
    // Return the expansion length of a given nonterminal.
    //=========================================================================
    std::uint64_t get_exp_len(const std::uint64_t id) const {
      return m_nonterminals[id].m_exp_len;
    }

    //=========================================================================
    // Return the Karp-Rabin hash a given nonterminal.
    //=========================================================================
    std::uint64_t get_kr_hash(const std::uint64_t id) const {
      return m_nonterminals[id].m_kr_hash;
    }

    //=========================================================================
    // Return the ID of the left child a given nonterminal.
    //=========================================================================
    std::uint64_t get_left_id(const std::uint64_t id) const {
      return m_nonterminals[id].m_left;
    }

    //=========================================================================
    // Return the ID of a right child a given nonterminal.
    //=========================================================================
    std::uint64_t get_right_id(const std::uint64_t id) const {
      return m_nonterminals[id].m_right;
    }

    //=========================================================================
    // Add a nonterminal.
    //=========================================================================
    std::uint64_t add_nonterminal(const nonterminal_type nonterm) {
      const std::uint64_t id = m_nonterminals.size();
      m_nonterminals.push_back(nonterm);

      // With probability 1/16 add to hash table.
      if (utils::random_int<std::uint64_t>(
            (std::uint64_t)0,
            (std::uint64_t)15) == 0)
        m_hashes.insert(get_kr_hash(id), id);

      // Return the id of the nonterminal.
      return id;
    }

    //=========================================================================
    // Decode the text and write to a given array.
    //=========================================================================
    void decode(
        char_type* &text,
        std::uint64_t &text_length) const {
      text_length = 0;
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it) {
        const std::uint64_t preflen = it->first;
        const std::uint64_t id = it->second;
        if (preflen != 0)
          text_length += get_exp_len(id);
      }
      text = new char_type[text_length];
      std::uint64_t ptr = 0;
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it) {
        const std::uint64_t preflen = it->first;
        const std::uint64_t id = it->second;
        if (preflen != 0) {
          m_nonterminals[id].write_expansion(
              id, text + ptr, this);
          ptr += get_exp_len(id);
        }
      }
    }

    //=========================================================================
    // Test the AVL property of all nonterminals.
    //=========================================================================
    bool test_avl_property() const {
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it) {
        const std::uint64_t preflen = it->first;
        const std::uint64_t id = it->second;
        if (preflen != 0 &&
            !m_nonterminals[id].test_avl_property(id, this))
          return false;
      }
      return true;
    }

    //=========================================================================
    // Collect Mersenne Karp-Rabin hashes in a vector.
    //=========================================================================
    void collect_mersenne_karp_rabin_hashes(
        std::vector<std::uint64_t> &hashes) const {
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it) {
        const std::uint64_t preflen = it->first;
        const std::uint64_t id = it->second;
        if (preflen != 0)
          (void) m_nonterminals[id].collect_mersenne_karp_rabin_hashes(
              id, hashes, this);
      }
    }

    //=========================================================================
    // Collect Mersenne Karp-Rabin hashes in a hash table.
    //=========================================================================
    void collect_mersenne_karp_rabin_hashes_2(
        hash_table<text_offset_type, std::uint64_t> &hashes) const {
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it) {
        const std::uint64_t preflen = it->first;
        const std::uint64_t id = it->second;
        if (preflen != 0)
          (void) m_nonterminals[id].collect_mersenne_karp_rabin_hashes_2(
             id, hashes, this);
      }
    }

    //=========================================================================
    // Count nonterminals in the pruned grammar.
    //=========================================================================
    void count_nonterminals_in_pruned_grammar(
        hash_table<text_offset_type, std::uint64_t> &hashes,
        hash_table<std::uint64_t, bool> &seen_hashes,
        std::uint64_t &current_count) const {
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it) {
        const std::uint64_t preflen = it->first;
        const std::uint64_t id = it->second;
        if (preflen != 0)
          m_nonterminals[id].count_nonterminals_in_pruned_grammar(
            id, hashes, seen_hashes, current_count, this);
      }
    }

    //=========================================================================
    // Collect pointers to all nonterminals reachable from the root.
    //=========================================================================
    void collect_nonterminal_pointers(
        std::vector<text_offset_type> &pointers) const {
      for (const_iter_type it = m_roots.begin(); it != m_roots.end(); ++it) {
        const std::uint64_t preflen = it->first;
        const std::uint64_t id = it->second;
        if (preflen != 0)
          m_nonterminals[id].collect_nonterminal_pointers(
              id, pointers, this);
      }
    }

    //=========================================================================
    // Merge roots enclosed in [begin..end).
    //=========================================================================
    void merge_enclosed_roots(
        const std::uint64_t begin,
        const std::uint64_t end) {

      // Find range [it_begin..it_end) to merge.
      iter_type it_begin = m_roots.end();
      it_begin = m_roots.lower_bound(begin);
      std::uint64_t begin2 = it_begin->first;
      ++it_begin;
      iter_type it_end = it_begin;
      while (it_end != m_roots.end() && it_end->first <= end)
        ++it_end;

      // Merge roots in [it_begin..it_end).
      if (it_begin != it_end) {
        std::vector<text_offset_type> v;
        for (iter_type it = it_begin; it != it_end; ++it)
          v.push_back(it->second);
        const std::uint64_t newroot_id = greedy_merge(v);
        const std::uint64_t end2 = begin2 + get_exp_len(newroot_id);

        // Update roots.
        m_roots.erase(it_begin, it_end);
        m_roots[end2] = newroot_id;
      }
    }

    //=========================================================================
    // Given two nonterminals `left' and `right' expanding to X and Y, add
    // nonterminals that expands to XY, and return the pointer to it.
    //=========================================================================
    std::uint64_t add_concat_nonterminal(
        const std::uint64_t left_id,
        const std::uint64_t right_id) {

      // Consider two cases, depending on whether
      // left of right nonterminal is taller.
      if (get_height(left_id) >= get_height(right_id)) {
        if (get_height(left_id) - get_height(right_id) <= 1) {

          // Height are close. Just merge and return.
          const nonterminal_type newroot(left_id, right_id, this);
          const std::uint64_t newroot_id = add_nonterminal(newroot);
          return newroot_id;
      } else {
        const std::uint64_t newright_id = 
          add_concat_nonterminal(get_right_id(left_id), right_id);
        if (get_height(newright_id) > get_height(get_left_id(left_id)) &&
            get_height(newright_id) - get_height(get_left_id(left_id)) > 1) {

            // Rebalancing needed.
            if (get_height(get_left_id(newright_id)) >
                get_height(get_right_id(newright_id))) {

              // Double (right-left) rotation.
              const nonterminal_type X(
                  get_left_id(left_id),
                  get_left_id(get_left_id(newright_id)), this);
              const nonterminal_type Z(
                  get_right_id(get_left_id(newright_id)),
                  get_right_id(newright_id), this);
              const std::uint64_t X_id = add_nonterminal(X);
              const std::uint64_t Z_id = add_nonterminal(Z);
              const nonterminal_type Y(X_id, Z_id, this);
              const std::uint64_t Y_id = add_nonterminal(Y);
              return Y_id;
            } else {

              // Single (left) rotation.
              const nonterminal_type X(
                  get_left_id(left_id), get_left_id(newright_id), this);
              const std::uint64_t X_id = add_nonterminal(X);
              const nonterminal_type Y(X_id, get_right_id(newright_id), this);
              const std::uint64_t Y_id = add_nonterminal(Y);
              return Y_id;
            }
          } else {

            // No need to rebalance.
            const nonterminal_type newroot(
                get_left_id(left_id), newright_id, this);
            const std::uint64_t newroot_id = add_nonterminal(newroot);
            return newroot_id;
          }
        }
      } else {
        if (get_height(right_id) - get_height(left_id) <= 1) {

          // Heights are close. Just merge and return.
          const nonterminal_type newroot(left_id, right_id, this);
          const std::uint64_t newroot_id = add_nonterminal(newroot);
          return newroot_id;
        } else {
          const std::uint64_t newleft_id =
            add_concat_nonterminal(left_id, get_left_id(right_id));
          if (get_height(newleft_id) > get_height(get_right_id(right_id)) &&
              get_height(newleft_id) - get_height(get_right_id(right_id)) > 1) {

            // Rebalancing needed.
            if (get_height(get_right_id(newleft_id)) >
                get_height(get_left_id(newleft_id))) {

              // Double (left-right) rotation.
              const nonterminal_type X(
                  get_left_id(newleft_id),
                  get_left_id(get_right_id(newleft_id)), this);
              const nonterminal_type Z(
                  get_right_id(get_right_id(newleft_id)),
                  get_right_id(right_id), this);
              const std::uint64_t X_id = add_nonterminal(X);
              const std::uint64_t Z_id = add_nonterminal(Z);
              const nonterminal_type Y(X_id, Z_id, this);
              const std::uint64_t Y_id = add_nonterminal(Y);
              return Y_id;
            } else {

              // Single (right) rotation.
              const nonterminal_type Y(
                  get_right_id(newleft_id), get_right_id(right_id), this);
              const std::uint64_t Y_id = add_nonterminal(Y);
              const nonterminal_type X(get_left_id(newleft_id), Y_id, this);
              const std::uint64_t X_id = add_nonterminal(X);
              return X_id;
            }
          } else {

            // No need to rebalance.
            const nonterminal_type newroot(
                newleft_id, get_right_id(right_id), this);
            const std::uint64_t newroot_id = add_nonterminal(newroot);
            return newroot_id;
          }
        }
      }
    }

    //=========================================================================
    // Get the sequence of nonterminals expanding to T[begin..end).
    //=========================================================================
    std::vector<text_offset_type> decomposition(
        std::uint64_t begin,
        std::uint64_t end) const {

      // Find leftmost root whose expansion overlaps/touches T[begin..end).
      const_iter_type it = m_roots.lower_bound(begin);

      // Proper substring or suffix of expansion of `it'.
      std::vector<text_offset_type> ret;
      if (begin < (std::uint64_t)(it->first)) {
        const std::uint64_t preflen = it->first;
        const std::uint64_t id = it->second;
        const std::uint64_t it_exp_size = get_exp_len(id);
        const std::uint64_t it_exp_beg = preflen - it_exp_size;
        const std::uint64_t local_beg = begin - it_exp_beg;
        const std::uint64_t local_end = std::min(preflen, end) - it_exp_beg;
        const std::uint64_t local_size = local_end - local_beg;
        std::vector<text_offset_type> dec =
          m_nonterminals[id].decomposition(id, local_beg, local_end, this);
        ret.insert(ret.end(), dec.begin(), dec.end());
        begin += local_size;
      }

      // Full expansions of nonterminals.
      ++it;
      while (begin < end && (std::uint64_t)(it->first) <= end) {
        ret.push_back(it->second);
        begin = (std::uint64_t)(it->first);
        ++it;
      }

      // Proper suffix of expansion of `it'.
      if (begin < end) {
        const std::uint64_t preflen = it->first;
        const std::uint64_t id = it->second;
        const std::uint64_t it_exp_size = get_exp_len(id);
        const std::uint64_t it_exp_beg = preflen - it_exp_size;
        const std::uint64_t local_end = end - it_exp_beg;
        std::vector<text_offset_type> dec =
          m_nonterminals[id].decomposition(id, 0, local_end, this);
        ret.insert(ret.end(), dec.begin(), dec.end());
        begin = end;
      }

      // Return result.
      return ret;
    }

    //=========================================================================
    // Return the sequence of nonterminals with the same expansion as seq.
    //=========================================================================
    std::vector<text_offset_type> find_equivalent_seq(
      const std::vector<text_offset_type> &seq) const {

      // Create the vector to hold the solution.
      std::vector<text_offset_type> ret;

      // Handle special case.
      if (seq.empty())
        return ret;

      // Allocate the arrays used in the dynamic programming.
      std::uint64_t length = seq.size();
      std::uint64_t *kr_hashes = utils::allocate_array<std::uint64_t>(length);
      std::uint64_t *exp_lengths = utils::allocate_array<std::uint64_t>(length);
      std::uint64_t **dp = utils::allocate_array<std::uint64_t*>(length);
      std::uint64_t **dp_sol = utils::allocate_array<std::uint64_t*>(length);
      text_offset_type **dp_nonterm = utils::allocate_array<text_offset_type*>(length);
      for (std::uint64_t i = 0; i < length; ++i) {
        dp[i] = utils::allocate_array<std::uint64_t>(length);
        dp_sol[i] = utils::allocate_array<std::uint64_t>(length);
        dp_nonterm[i] = utils::allocate_array<text_offset_type>(length);
      }

      // Fill in the array for len = 1.
      for (std::uint64_t i = 0; i < length; ++i) {
        dp[i][i] = 1;
        dp_sol[i][i] = 1;
        dp_nonterm[i][i] = seq[i];
        kr_hashes[i] = get_kr_hash((std::uint64_t)seq[i]);
        exp_lengths[i] = get_exp_len((std::uint64_t)seq[i]);
      }

      // Solve for subarray of length > 1.
      for (std::uint64_t len = 2; len <= length; ++len) {
        for (std::uint64_t beg = 0; beg <= length - len; ++beg) {
          const std::uint64_t end = beg + len - 1;

          // Initialize to solution to initial choice.
          dp[beg][end] = 1 + dp[beg + 1][end];
          dp_sol[beg][end] = 1;
          dp_nonterm[beg][end] = seq[beg];

          // Try all other possible choices.
          std::uint64_t h = kr_hashes[beg];
          for (std::uint64_t leftlen = 2;
              leftlen <= len; ++leftlen) {
            const std::uint64_t last = beg + leftlen - 1;
            const std::uint64_t right_hash = kr_hashes[last];
            const std::uint64_t right_len = exp_lengths[last];
            h = karp_rabin_hashing::concat(h, right_hash, right_len);
            const text_offset_type *nonterm_id_ptr = m_hashes.find(h);
            if (nonterm_id_ptr != NULL) {
              std::uint64_t sol_cost = 1;
              if (leftlen < len) sol_cost += dp[last + 1][end];
              if (sol_cost < dp[beg][end]) {
                dp[beg][end] = sol_cost;
                dp_sol[beg][end] = leftlen;
                dp_nonterm[beg][end] = *nonterm_id_ptr;
              }
            }
          }
        }
      }

      // Restore the optimal solution.
      std::uint64_t prefix_length = 0;
      while (prefix_length < length) {
        ret.push_back(dp_nonterm[prefix_length][length - 1]);
        prefix_length += dp_sol[prefix_length][length - 1];
      }

      // Clean up.
      for (std::uint64_t i = 0; i < length; ++i) {
        utils::deallocate(dp[i]);
        utils::deallocate(dp_sol[i]);
        utils::deallocate(dp_nonterm[i]);
      }
      utils::deallocate(dp);
      utils::deallocate(dp_sol);
      utils::deallocate(dp_nonterm);
      utils::deallocate(exp_lengths);
      utils::deallocate(kr_hashes);

      // Return the result.
      return ret;
    }

  private:

    //=========================================================================
    // Merge greedily (shortest first) sequence of nonterminals.
    //=========================================================================
    std::uint64_t greedy_merge(
        std::vector<text_offset_type> &seq) {
      while (seq.size() > 1) {

        // Find the nonterminal with the smallest height.
        std::uint64_t smallest_height_id = 0;
        for (std::uint64_t i = 1; i < seq.size(); ++i) {
          if (get_height((std::uint64_t)seq[i]) <
              get_height((std::uint64_t)seq[smallest_height_id]))
            smallest_height_id = i;
        }

        // Merge the nonterminal with the smaller height with
        // one of its beighbors (whichever is shorter).
        if (smallest_height_id == 0 ||
            (smallest_height_id + 1 < seq.size() &&
             get_height((std::uint64_t)seq[smallest_height_id + 1]) <=
             get_height((std::uint64_t)seq[smallest_height_id - 1]))) {

          // Only right neighbor exists, or both exist
          // and the right one is not taller than the left
          // one. End result: merge with the right neighbor.
          const std::uint64_t left_id = seq[smallest_height_id];
          const std::uint64_t right_id = seq[smallest_height_id + 1];
          const std::uint64_t left_hash = get_kr_hash(left_id);
          const std::uint64_t right_hash = get_kr_hash(right_id);
          const std::uint64_t right_len = get_exp_len(right_id);
          const std::uint64_t h =
            karp_rabin_hashing::concat(left_hash, right_hash, right_len);
          seq.erase(seq.begin() + smallest_height_id);
          if (m_hashes.find(h) != NULL)
            seq[smallest_height_id] = *(m_hashes.find(h));
          else
            seq[smallest_height_id] =
              add_concat_nonterminal(left_id, right_id);
        } else {

          // Only left neighbor exists, or both exists
          // and the left one is not taller than the
          // right one. End result: merge with left neighbor.
          const std::uint64_t left_id = seq[smallest_height_id - 1];
          const std::uint64_t right_id = seq[smallest_height_id];
          const std::uint64_t left_hash = get_kr_hash(left_id);
          const std::uint64_t right_hash = get_kr_hash(right_id);
          const std::uint64_t right_len = get_exp_len(right_id);
          const std::uint64_t h =
            karp_rabin_hashing::concat(left_hash, right_hash, right_len);
          seq.erase(seq.begin() + (smallest_height_id - 1));
          if (m_hashes.find(h) != NULL)
            seq[smallest_height_id - 1] = *(m_hashes.find(h));
          else
            seq[smallest_height_id - 1] =
              add_concat_nonterminal(left_id, right_id);
        }
      }
      return seq[0];
    }
};

//=============================================================================
// Default constructor.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal()
  : m_char((char_type)0),
    m_height(0),
    m_exp_len(1),
    m_kr_hash(0),
    m_left(std::numeric_limits<text_offset_type>::max()),
    m_right(std::numeric_limits<text_offset_type>::max()) {}

//=============================================================================
// Constructor for a nonterminal expanding to a single symbol.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(const char_type c)
  : m_char(c),
    m_height(0),
    m_exp_len(1),
    m_kr_hash(karp_rabin_hashing::hash_char(c)),
    m_left(std::numeric_limits<text_offset_type>::max()),
    m_right(std::numeric_limits<text_offset_type>::max()) {}

//=============================================================================
// Constructor for non-single-symbol nonterminal.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(
      const std::uint64_t left,
      const std::uint64_t right,
      const avl_grammar_multiroot<char_type, text_offset_type> * const g)
  : m_char((char_type)0),
    m_height(std::max(g->get_height(left), g->get_height(right)) + 1),
    m_exp_len(g->get_exp_len(left) + g->get_exp_len(right)),
    m_kr_hash(karp_rabin_hashing::concat(
        g->get_kr_hash(left),
        g->get_kr_hash(right),
        g->get_exp_len(right))),
    m_left(left),
    m_right(right) {}

//=============================================================================
// Print expansion of a given nonterminal.
//=============================================================================
template<typename char_type, typename text_offset_type>
void nonterminal<char_type, text_offset_type>::print_expansion(
    const std::uint64_t id,
    const avl_grammar_multiroot<char_type, text_offset_type> * const g) const {
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  const std::uint64_t height = g->get_height(id);
  if (height == 0) {
    const char_type my_char = g->get_char(id);
    fprintf(stderr, "%c", (char)my_char);
  } else {
    const std::uint64_t left_id = g->get_left_id(id);
    const std::uint64_t right_id = g->get_right_id(id);
    const nonterminal_type &left = g->get_nonterminal(left_id);
    const nonterminal_type &right = g->get_nonterminal(right_id);
    left.print_expansion(left_id, g);
    right.print_expansion(right_id, g);
  }
}

//=============================================================================
// Write the expansion into the given array.
//=============================================================================
template<typename char_type, typename text_offset_type>
void nonterminal<char_type, text_offset_type>::write_expansion(
    const std::uint64_t id,
    char_type * const text,
    const avl_grammar_multiroot<char_type, text_offset_type> * const g) const {
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  const std::uint64_t height = g->get_height(id);
  if (height == 0) {
    const char_type my_char = g->get_char(id);
    text[0] = my_char;
  } else {
    const std::uint64_t left_id = g->get_left_id(id);
    const std::uint64_t right_id = g->get_right_id(id);
    const std::uint64_t left_exp_len = g->get_exp_len(left_id);
    const nonterminal_type &left = g->get_nonterminal(left_id);
    const nonterminal_type &right = g->get_nonterminal(right_id);
    left.write_expansion(left_id, text, g);
    right.write_expansion(right_id, text + left_exp_len, g);
  }
}

//=============================================================================
// Test the AVL propert of a subtree.
//=============================================================================
template<typename char_type, typename text_offset_type>
bool nonterminal<char_type, text_offset_type>::test_avl_property(
    const std::uint64_t id,
    const avl_grammar_multiroot<char_type, text_offset_type> * const g) const {
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  const std::uint64_t height = g->get_height(id);
  if (height == 0)
    return true;

  const std::uint64_t left_id = g->get_left_id(id);
  const std::uint64_t right_id = g->get_right_id(id);
  const nonterminal_type &left = g->get_nonterminal(left_id);
  const nonterminal_type &right = g->get_nonterminal(right_id);
  if (!left.test_avl_property(left_id, g) ||
      !right.test_avl_property(right_id, g))
    return false;

  const std::uint64_t left_height = g->get_height(left_id);
  const std::uint64_t right_height = g->get_height(right_id);
  if ((right_height > left_height && right_height - left_height > 1) ||
      (right_height < left_height && left_height - right_height > 1))
    return false;

  return true;
}

//=============================================================================
// Collect Mersenne Karp-Rabin hashes of all nonterminals.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>
::collect_mersenne_karp_rabin_hashes(
    const std::uint64_t id,
    std::vector<std::uint64_t> &hashes,
    const avl_grammar_multiroot<char_type, text_offset_type> * const g) const {
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  const std::uint64_t height = g->get_height(id);
  if (height == 0) {
    const char_type my_char = g->get_char(id);
    const std::uint64_t h = karp_rabin_hashing::hash_char(my_char);
    hashes.push_back(h);
    return h;
  } else {
    const std::uint64_t left_id = g->get_left_id(id);
    const std::uint64_t right_id = g->get_right_id(id);
    const nonterminal_type &left = g->get_nonterminal(left_id);
    const nonterminal_type &right = g->get_nonterminal(right_id);
    const std::uint64_t left_hash =
      left.collect_mersenne_karp_rabin_hashes(left_id, hashes, g);
    const std::uint64_t right_hash =
      right.collect_mersenne_karp_rabin_hashes(right_id, hashes, g);
    const std::uint64_t right_len = g->get_exp_len(right_id);
    const std::uint64_t h = karp_rabin_hashing::concat(
        left_hash, right_hash, right_len);
    hashes.push_back(h);
    return h;
  }
}

//=============================================================================
// Collect reachable nonterminals.
//=============================================================================
template<typename char_type, typename text_offset_type>
void nonterminal<char_type, text_offset_type>::collect_nonterminal_pointers(
    const std::uint64_t id,
    std::vector<text_offset_type> &pointers,
    const avl_grammar_multiroot<char_type, text_offset_type> * const g) const {
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  const std::uint64_t height = g->get_height(id);
  pointers.push_back(id);
  if (height > 0) {
    const std::uint64_t left_id = g->get_left_id(id);
    const std::uint64_t right_id = g->get_right_id(id);
    const nonterminal_type &left = g->get_nonterminal(left_id);
    const nonterminal_type &right = g->get_nonterminal(right_id);
    left.collect_nonterminal_pointers(left_id, pointers, g);
    right.collect_nonterminal_pointers(right_id, pointers, g);
 }
}

//=============================================================================
// Collect Mersenne Karp-Rabin hashes of all nonterminals.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>
::collect_mersenne_karp_rabin_hashes_2(
    const std::uint64_t id,
    hash_table<text_offset_type, std::uint64_t> &hashes,
    const avl_grammar_multiroot<char_type, text_offset_type> * const g) const {
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  const std::uint64_t height = g->get_height(id);
  if (height == 0) {
    const char_type my_char = g->get_char(id);
    const std::uint64_t h = karp_rabin_hashing::hash_char(my_char);
    hashes.insert(id, h);
    return h;
  } else {
    const std::uint64_t left_id = g->get_left_id(id);
    const std::uint64_t right_id = g->get_right_id(id);
    const nonterminal_type &left = g->get_nonterminal(left_id);
    const nonterminal_type &right = g->get_nonterminal(right_id);
    const std::uint64_t left_hash =
      left.collect_mersenne_karp_rabin_hashes_2(left_id, hashes, g);
    const std::uint64_t right_hash =
      right.collect_mersenne_karp_rabin_hashes_2(right_id, hashes, g);
    const std::uint64_t right_len = g->get_exp_len(right_id);
    const std::uint64_t h = karp_rabin_hashing::concat(
        left_hash, right_hash, right_len);
    hashes.insert(id, h);
    return h;
  }
}

//=============================================================================
// Compute the number of nonterminals in the pruned grammar.
//=============================================================================
template<typename char_type, typename text_offset_type>
void nonterminal<char_type, text_offset_type>
::count_nonterminals_in_pruned_grammar(
    const std::uint64_t id,
    hash_table<text_offset_type, std::uint64_t> &hashes,
    hash_table<std::uint64_t, bool> &seen_hashes,
    std::uint64_t &current_count,
    const avl_grammar_multiroot<char_type, text_offset_type> * const g) const {
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  const std::uint64_t * const h = hashes.find((text_offset_type)id);
  if (seen_hashes.find(*h) == NULL) {
    seen_hashes.insert(*h, true);
    ++current_count;
    const std::uint64_t height = g->get_height(id);
    if (height != 0) {
      const std::uint64_t left_id = g->get_left_id(id);
      const std::uint64_t right_id = g->get_right_id(id);
      const nonterminal_type &left = g->get_nonterminal(left_id);
      const nonterminal_type &right = g->get_nonterminal(right_id);
      left.count_nonterminals_in_pruned_grammar(
          left_id, hashes, seen_hashes, current_count, g);
      right.count_nonterminals_in_pruned_grammar(
          right_id, hashes, seen_hashes, current_count, g);
    }
  }
}

//=============================================================================
// Assuming S is the expansions of the nontermnal, return the
// sequence of nonterminals expanding to S[begin..end).
//=============================================================================
template<typename char_type, typename text_offset_type>
std::vector<text_offset_type>
nonterminal<char_type, text_offset_type>::decomposition(
    const std::uint64_t id,
    const std::uint64_t begin,
    const std::uint64_t end,
    const avl_grammar_multiroot<char_type, text_offset_type> * const g) const {

  // Declare the vector storing the result.
  std::vector<text_offset_type> ret;

  // Handle boundary case.
  if (begin == end)
    return ret;

  // Find the deepest nonterminal in the parse tree containing the range
  // [begin..end).
  std::uint64_t x = id;
  std::uint64_t x_height = g->get_height(x);
  std::uint64_t cur_range_beg = 0;
  std::uint64_t cur_range_end = g->get_exp_len(x);
  while (x_height > 0) {
    const std::uint64_t x_left = g->get_left_id(x);
    const std::uint64_t x_right = g->get_right_id(x);
    const std::uint64_t x_left_exp_len = g->get_exp_len(x_left);
    const std::uint64_t cur_range_mid = cur_range_beg + x_left_exp_len;
    if (end <= cur_range_mid) {
      cur_range_end = cur_range_mid;
      x = x_left;
      x_height = g->get_height(x);
    } else if (begin >= cur_range_mid) {
      cur_range_beg = cur_range_mid;
      x = x_right;
      x_height = g->get_height(x);
    } else break;
  }

  // Check if the range of x is exactly [begin..end).
  if (cur_range_beg == begin && cur_range_end == end) {

    // If yes, return x as the answer.
    ret.push_back((text_offset_type)x);
  } else {

    // Otherwise, we perform two traversals in the tree.
    {
      const std::uint64_t x_left = g->get_left_id(x);
      const std::uint64_t x_left_exp_len = g->get_exp_len(x_left);
      const std::uint64_t left_range_end = cur_range_beg + x_left_exp_len;
      std::uint64_t suffix_length = left_range_end - begin;
      std::uint64_t y = x_left;
      while (suffix_length > 0) {
        const std::uint64_t y_exp_len = g->get_exp_len(y);
        const std::uint64_t y_left = g->get_left_id(y);
        const std::uint64_t y_right = g->get_right_id(y);
        if (y_exp_len == suffix_length) {
          ret.push_back(y);
          suffix_length -= y_exp_len;
        } else {
          const std::uint64_t y_right_exp_len = g->get_exp_len(y_right);
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
      const std::uint64_t x_left = g->get_left_id(x);
      const std::uint64_t x_right = g->get_right_id(x);
      const std::uint64_t x_left_exp_len = g->get_exp_len(x_left);
      const std::uint64_t right_range_beg = cur_range_beg + x_left_exp_len;
      std::uint64_t prefix_length = end - right_range_beg;
      std::uint64_t y = x_right;
      while (prefix_length > 0) {
        const std::uint64_t y_exp_len = g->get_exp_len(y);
        const std::uint64_t y_left = g->get_left_id(y);
        const std::uint64_t y_right = g->get_right_id(y);
        if (y_exp_len == prefix_length) {
          ret.push_back(y);
          prefix_length -= y_exp_len;
        } else {
          const std::uint64_t y_left_exp_len = g->get_exp_len(y_left);
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

#endif  // __AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
