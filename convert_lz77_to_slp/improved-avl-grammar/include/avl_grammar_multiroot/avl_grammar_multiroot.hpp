#ifndef __AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
#define __AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED

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


//=============================================================================
// Class used to represent multiroot AVL grammar. Forward declaration.
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
    typedef packed_pair<text_offset_type, text_offset_type> pair_type;

  public:
    std::uint8_t m_height;
    std::uint8_t m_exp_len;
    text_offset_type m_left;
    text_offset_type m_right;

    nonterminal();
    nonterminal(const char_type);
    nonterminal(const std::uint64_t, const std::uint64_t,
        const grammar_type * const g);
    nonterminal(const nonterminal_type &);
    void print_expansion( const std::uint64_t,
        const grammar_type * const) const;
    std::uint64_t write_expansion(const std::uint64_t, char_type * const,
        const grammar_type * const) const;
    bool compare_expansion_to_text(const std::uint64_t,
        const char_type * const, const grammar_type * const) const;
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
    void decomposition(const std::uint64_t, const std::uint64_t,
        const std::uint64_t, space_efficient_vector<pair_type> &,
        const grammar_type * const g) const;
} __attribute__((packed));

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
  static_assert(sizeof(char_type) <= sizeof(text_offset_type),
      "Error: sizeof(char_type) > sizeof(text_offset_type)!");

  //===========================================================================
  // Declare typedefs
  //===========================================================================
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  typedef packed_pair<text_offset_type, text_offset_type> pair_type;
  typedef packed_triple<text_offset_type, text_offset_type, std::uint64_t>
    triple_type;
  typedef packed_pair<text_offset_type, std::uint64_t> hash_pair_type;

  private:

    //=========================================================================
    // Class members.
    //=========================================================================
    space_efficient_vector<pair_type> m_roots_vec;
    space_efficient_vector<nonterminal_type> m_nonterminals;
    space_efficient_vector<pair_type> m_long_exp_len;
    space_efficient_vector<hash_pair_type> m_long_exp_hashes;
    hash_table<std::uint64_t, text_offset_type, text_offset_type> m_hashes;
    cache<text_offset_type, std::uint64_t> *m_kr_hash_cache;
    char_type *m_snippet;
    std::uint64_t empty_step_counter;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    avl_grammar_multiroot() {
      const std::uint64_t cache_size = (1 << 17);
      empty_step_counter = 0;
      m_roots_vec.push_back(pair_type(
            (text_offset_type)0,
            (text_offset_type)std::numeric_limits<text_offset_type>::max()));
      m_snippet = utils::allocate_array<char_type>(256);
      m_kr_hash_cache =
        new cache<text_offset_type, std::uint64_t>(cache_size);
    }

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~avl_grammar_multiroot() {
      utils::deallocate(m_snippet);
      delete m_kr_hash_cache;
    }

    void print_stats() const {

      // Print RAM usage breakdown.
      const std::uint64_t m_roots_vec_ram_use = m_roots_vec.ram_use();
      const std::uint64_t m_nonterminals_ram_use = m_nonterminals.ram_use();
      const std::uint64_t m_long_exp_hashes_ram_use = m_long_exp_hashes.ram_use();
      const std::uint64_t m_long_exp_len_ram_use = m_long_exp_len.ram_use();
      const std::uint64_t m_hashes_ram_use = m_hashes.ram_use();
      const std::uint64_t m_kr_hash_cache_ram_use = m_kr_hash_cache->ram_use();
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
      fprintf(stderr, "  m_long_exp_hashes: %.2LfMiB (%.2Lf%%)\n",
          (1.L * m_long_exp_hashes_ram_use) / (1 << 20),
          (100.L * m_long_exp_hashes_ram_use) / total);
      fprintf(stderr, "  m_long_exp_len: %.2LfMiB (%.2Lf%%)\n",
          (1.L * m_long_exp_len_ram_use) / (1 << 20),
          (100.L * m_long_exp_len_ram_use) / total);
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
      m_kr_hash_cache->print_cache_miss_rate();
    }

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
          std::numeric_limits<text_offset_type>::max()) {
        ++beg;
        ++empty_step_counter;
      }

      // Return the result.
      return beg;
    }

    //=========================================================================
    // Find the first undeleted root.
    //=========================================================================
    inline std::uint64_t roots_begin() const {

      // We use the fact that there is a sentinel at the begnning.
      return 0;
    }

    //=========================================================================
    // Return the past-the-end position in the roots array.
    //=========================================================================
    inline std::uint64_t roots_end() const {

      // We use the fact that there is a sentinel at the begnning.
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
          std::numeric_limits<text_offset_type>::max()) {
        ++pos;
        ++empty_step_counter;
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

    //=========================================================================
    // Run garbage collector if needed.
    //=========================================================================
    void check_gargage_collector() {
      if (empty_step_counter > m_roots_vec.size()) {
        roots_garbage_collector();
        empty_step_counter = 0;
      }
    }

    //=========================================================================
    // Print the string encoded by the grammar.
    //=========================================================================
    void print_expansion() {
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i)) {
        const std::uint64_t preflen = m_roots_vec[i].first;
        const std::uint64_t id = m_roots_vec[i].second;
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
        const std::uint64_t id) {
      m_roots_vec.push_back(pair_type(pos, id));
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
      const nonterminal_type &nonterm = get_nonterminal(id);
      return nonterm.m_height;
    }

    //=========================================================================
    // Return the char associated with a given nonterminal.
    //=========================================================================
    char_type get_char(const std::uint64_t id) const {
      const nonterminal_type &nonterm = get_nonterminal(id);
      if (nonterm.m_height == 0) {
        char_type c = (char_type)nonterm.m_left;
        return c;
      } else return (char_type)0;
    }

    //=========================================================================
    // Return the expansion length of a given nonterminal.
    //=========================================================================
    std::uint64_t get_exp_len(const std::uint64_t id) const {
      const nonterminal_type &nonterm = get_nonterminal(id);
      if (nonterm.m_exp_len == 255) {

        // Binary search in m_long_exp_len.
        std::uint64_t beg = 0;
        std::uint64_t end = m_long_exp_len.size();
        while (beg + 1 < end) {
          const std::uint64_t mid = (beg + end) / 2;
          const std::uint64_t mid_id = m_long_exp_len[mid].first;
          if (mid_id <= id)
            beg = mid;
          else end = mid;
        }
        return (std::uint64_t)m_long_exp_len[beg].second;
      } else return nonterm.m_exp_len;
    }

    //=========================================================================
    // Return the Karp-Rabin hash a given nonterminal.
    //=========================================================================
    std::uint64_t get_kr_hash(const std::uint64_t id) const {

      // Check, if the value is in cache.
      const std::uint64_t *cache_ret =
        m_kr_hash_cache->lookup((text_offset_type)id);
      if (cache_ret != NULL)
        return *cache_ret;

      // The value is not in cache.
      const nonterminal_type &nonterm = get_nonterminal(id);
      std::uint64_t ret = 0;
      if (nonterm.m_exp_len < 255) {

        // Recompute the hash from scratch.
        (void) nonterm.write_expansion(id, m_snippet, this);
        ret = karp_rabin_hashing::hash_string<char_type>(
            m_snippet, nonterm.m_exp_len);
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
        ret = m_long_exp_hashes[beg].second;
      }

      // Update cache.
      m_kr_hash_cache->insert(id, ret);

      // Return the result.
      return ret;
    }

    //=========================================================================
    // Return the ID of the left child a given nonterminal.
    //=========================================================================
    std::uint64_t get_left_id(const std::uint64_t id) const {
      const nonterminal_type &nonterm = get_nonterminal(id);
      return nonterm.m_left;
    }

    //=========================================================================
    // Return the ID of a right child a given nonterminal.
    //=========================================================================
    std::uint64_t get_right_id(const std::uint64_t id) const {
      const nonterminal_type &nonterm = get_nonterminal(id);
      return nonterm.m_right;
    }

    //=========================================================================
    // Add nonterminal expanding to single symbol.
    //=========================================================================
    std::uint64_t add_nonterminal(const nonterminal_type &nonterm) {
      const std::uint64_t id = m_nonterminals.size();
      m_nonterminals.push_back(nonterm);

      // With probability 1/16 add to hash table.
      if (utils::random_int<std::uint64_t>(
            (std::uint64_t)0,
            (std::uint64_t)7) == 0) {
        const std::uint64_t kr_hash = get_kr_hash(id);
        text_offset_type * const ret = m_hashes.find(kr_hash);
        if (ret == NULL)
          m_hashes.insert(kr_hash, id);
        else *ret = (text_offset_type)id;
      }

      // Return the id of the nonterminal.
      return id;
    }

    //=========================================================================
    // Add a new binary nonterminal.
    //=========================================================================
    std::uint64_t add_nonterminal(
        const std::uint64_t left_id,
        const std::uint64_t right_id) {

      // Compute values for the new nonterminal.
      const std::uint64_t left_exp_len = get_exp_len(left_id);
      const std::uint64_t right_exp_len = get_exp_len(right_id);
      const std::uint64_t new_exp_len = left_exp_len + right_exp_len;
      const std::uint8_t left_height = get_height(left_id);
      const std::uint8_t right_height = get_height(right_id);
      const std::uint8_t new_height = std::max(left_height, right_height) + 1;
      const std::uint64_t new_id = m_nonterminals.size();

      // Create and add new nonterminal.
      nonterminal_type new_nonterm;
      new_nonterm.m_height = new_height;
      new_nonterm.m_exp_len = std::min(255UL, new_exp_len);
      new_nonterm.m_left = left_id;
      new_nonterm.m_right = right_id;
      m_nonterminals.push_back(new_nonterm);

      // With probability 1/16 add to hash table.
      std::uint64_t new_kr_hash = 0;
      bool hash_computed = false;
      if (utils::random_int<std::uint64_t>(
            (std::uint64_t)0,
            (std::uint64_t)7) == 0) {
        const std::uint64_t left_hash = get_kr_hash(left_id);
        const std::uint64_t right_hash = get_kr_hash(right_id);
        new_kr_hash =
          karp_rabin_hashing::concat(left_hash, right_hash, right_exp_len);
        hash_computed = true;
        text_offset_type * const ret = m_hashes.find(new_kr_hash);
        if (ret == NULL)
          m_hashes.insert(new_kr_hash, new_id);
        else *ret = (text_offset_type)new_id;
      }

      // Update list of long nonterminals.
      if (new_exp_len >= 255) {
        m_long_exp_len.push_back(
            pair_type(
              (text_offset_type)new_id,
              (text_offset_type)new_exp_len));
      }

      if (new_exp_len >= 255) {
        if (!hash_computed) {
          const std::uint64_t left_hash = get_kr_hash(left_id);
          const std::uint64_t right_hash = get_kr_hash(right_id);
          new_kr_hash =
            karp_rabin_hashing::concat(left_hash, right_hash, right_exp_len);
        }

        m_long_exp_hashes.push_back(
            hash_pair_type(
              (text_offset_type)new_id,
              new_kr_hash));
      }

      // Return the id of the new nonterminal.
      return new_id;
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
        const std::uint64_t id = m_roots_vec[i].second;
        if (preflen != 0)
          text_length += get_exp_len(id);
      }
      text = new char_type[text_length];
      std::uint64_t ptr = 0;
      for (std::uint64_t i = roots_begin();
          i != roots_end(); i = roots_next(i)) {
        const std::uint64_t preflen = m_roots_vec[i].first;
        const std::uint64_t id = m_roots_vec[i].second;
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
        std::uint64_t preflen = m_roots_vec[i].first;
        std::uint64_t id = m_roots_vec[i].second;
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
        const std::uint64_t id = m_roots_vec[i].second;
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
        hash_table<text_offset_type, std::uint64_t> &hashes) {
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
        hash_table<text_offset_type, std::uint64_t> &hashes,
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
        std::vector<text_offset_type> &pointers) {
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
        v.push_back(triple_type(
              m_roots_vec[range_end].second,
              (text_offset_type)cur_exp_size,
              get_kr_hash(m_roots_vec[range_end].second)));
        m_roots_vec[range_end].second =
          std::numeric_limits<text_offset_type>::max();
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
          const std::uint64_t newroot_id =
            add_nonterminal(left_id, right_id);
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
              const std::uint64_t X_id = add_nonterminal(
                  get_left_id(left_id),
                  get_left_id(get_left_id(newright_id)));
              const std::uint64_t Z_id = add_nonterminal(
                  get_right_id(get_left_id(newright_id)),
                  get_right_id(newright_id));
              const std::uint64_t Y_id = add_nonterminal(
                  X_id, Z_id);
              return Y_id;
            } else {

              // Single (left) rotation.
              const std::uint64_t X_id = add_nonterminal(
                  get_left_id(left_id), get_left_id(newright_id));
              const std::uint64_t Y_id = add_nonterminal(
                  X_id, get_right_id(newright_id));
              return Y_id;
            }
          } else {

            // No need to rebalance.
            const std::uint64_t newroot_id = add_nonterminal(
                get_left_id(left_id), newright_id);
            return newroot_id;
          }
        }
      } else {
        if (get_height(right_id) - get_height(left_id) <= 1) {

          // Heights are close. Just merge and return.
          const std::uint64_t newroot_id =
            add_nonterminal(left_id, right_id);
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
              const std::uint64_t X_id = add_nonterminal(
                  get_left_id(newleft_id),
                  get_left_id(get_right_id(newleft_id)));
              const std::uint64_t Z_id = add_nonterminal(
                  get_right_id(get_right_id(newleft_id)),
                  get_right_id(right_id));
              const std::uint64_t Y_id = add_nonterminal(
                  X_id, Z_id);
              return Y_id;
            } else {

              // Single (right) rotation.
              const std::uint64_t Y_id = add_nonterminal(
                  get_right_id(newleft_id), get_right_id(right_id));
              const std::uint64_t X_id = add_nonterminal(
                  get_left_id(newleft_id), Y_id);
              return X_id;
            }
          } else {

            // No need to rebalance.
            const std::uint64_t newroot_id = add_nonterminal(
                newleft_id, get_right_id(right_id));
            return newroot_id;
          }
        }
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
        ret.push_back(pair_type(
              (text_offset_type)id,
              (text_offset_type)exp_len));
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
      if (seq.empty())
        return;

      // Create the vector to hold the solution.
      space_efficient_vector<pair_type> ret;

      // Allocate the arrays used in the dynamic programming.
      const std::uint64_t length = seq.size();
      std::uint64_t *kr_hashes = utils::allocate_array<std::uint64_t>(length);
      std::uint64_t **dp = utils::allocate_array<std::uint64_t*>(length);
      std::uint64_t **dp_sol = utils::allocate_array<std::uint64_t*>(length);
      std::uint64_t **dp_explen = utils::allocate_array<std::uint64_t*>(length);
      text_offset_type **dp_nonterm = utils::allocate_array<text_offset_type*>(length);
      for (std::uint64_t i = 0; i < length; ++i) {
        dp[i] = utils::allocate_array<std::uint64_t>(length);
        dp_sol[i] = utils::allocate_array<std::uint64_t>(length);
        dp_explen[i] = utils::allocate_array<std::uint64_t>(length);
        dp_nonterm[i] = utils::allocate_array<text_offset_type>(length);
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
            const text_offset_type *nonterm_id_ptr = m_hashes.find(h);
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
        ret.push_back(pair_type(
              (text_offset_type)id,
              (text_offset_type)exp_len));
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

  private:

    //=========================================================================
    // Heap down routine.
    //========================================================================
    void heap_down(
        std::uint64_t i,
        const space_efficient_vector<triple_type> &seq,
        text_offset_type * const heap,
        const std::uint64_t heap_size) const {
      ++i;
      std::uint64_t min_pos = i;
      while (true) {
        if ((i << 1) <= heap_size &&
            get_height(seq[heap[(i << 1) - 1]].first) <
            get_height(seq[heap[min_pos - 1]].first))
          min_pos = (i << 1);
        if ((i << 1) + 1 <= heap_size &&
            get_height(seq[heap[i << 1]].first) <
            get_height(seq[heap[min_pos - 1]].first))
          min_pos = (i << 1) + 1;
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
    // Make heap rountine.
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
    std::uint64_t greedy_merge(
        space_efficient_vector<triple_type> &seq) {

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
      std::uint64_t ret = 0;
      while (true) {
        const std::uint64_t min_elem = heap[0];

        // If the element was already deleted, skip it.
        if (deleted[min_elem]) {
          (void) extract_min(seq, heap, heap_size);
          continue;
        }

        // Merge min_elem with one of its
        // beighbors (whichever is shorter).
        if ((std::uint64_t)prev[min_elem] == sentinel &&
            (std::uint64_t)next[min_elem] == sentinel) {

          // We have reached the only nonterminal.
          // Store it as output and exit the loop.
          ret = seq[min_elem].first;
          break;
        } else if ((std::uint64_t)prev[min_elem] == sentinel ||
            ((std::uint64_t)next[min_elem] != sentinel &&
             get_height(seq[next[min_elem]].first) <=
             get_height(seq[prev[min_elem]].first))) {

          // Only right neighbor exists, or both exist
          // and the right one is not taller than the left
          // one. End result: merge with the right neighbor.
          const std::uint64_t right_elem = next[min_elem];
          const std::uint64_t left_id = seq[min_elem].first;
          const std::uint64_t right_id = seq[right_elem].first;
          const std::uint64_t left_hash = seq[min_elem].third;
          const std::uint64_t right_hash = seq[right_elem].third;
          const std::uint64_t left_len = seq[min_elem].second;
          const std::uint64_t right_len = seq[right_elem].second;
          const std::uint64_t merged_len = left_len + right_len;
          const std::uint64_t h =
            karp_rabin_hashing::concat(left_hash, right_hash, right_len);
          std::uint64_t id_merged = 0;
          const text_offset_type * const hash_ret = m_hashes.find(h);
          if (hash_ret != NULL)
            id_merged = *hash_ret;
          else
            id_merged = add_concat_nonterminal(left_id, right_id);
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
          const std::uint64_t left_id = seq[left_elem].first;
          const std::uint64_t right_id = seq[min_elem].first;
          const std::uint64_t left_hash = seq[left_elem].third;
          const std::uint64_t right_hash = seq[min_elem].third;
          const std::uint64_t left_len = seq[left_elem].second;
          const std::uint64_t right_len = seq[min_elem].second;
          const std::uint64_t merged_len = left_len + right_len;
          const std::uint64_t h =
            karp_rabin_hashing::concat(left_hash, right_hash, right_len);
          std::uint64_t id_merged = 0;
          const text_offset_type * const hash_ret = m_hashes.find(h);
          if (hash_ret != NULL)
            id_merged = *hash_ret;
          else
            id_merged = add_concat_nonterminal(left_id, right_id);
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
    m_exp_len(1),
    m_left(std::numeric_limits<text_offset_type>::max()),
    m_right(std::numeric_limits<text_offset_type>::max()) {}

//=============================================================================
// Constructor for a nonterminal expanding to a single symbol.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(const char_type c)
  : m_height(0),
    m_exp_len(1),
    m_left((text_offset_type)c),
    m_right(std::numeric_limits<text_offset_type>::max()) {}

//=============================================================================
// Copy constructor.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(
    const nonterminal<char_type, text_offset_type> &x)
  : m_height(x.m_height),
    m_exp_len(x.m_exp_len),
    m_left(x.m_left),
    m_right(x.m_right) {}

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
std::uint64_t nonterminal<char_type, text_offset_type>::write_expansion(
    const std::uint64_t id,
    char_type * const text,
    const avl_grammar_multiroot<char_type, text_offset_type> * const g) const {
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  const std::uint64_t height = g->get_height(id);
  if (height == 0) {
    const char_type my_char = g->get_char(id);
    text[0] = my_char;
    return 1;
  } else {
    const std::uint64_t left_id = g->get_left_id(id);
    const std::uint64_t right_id = g->get_right_id(id);
    const nonterminal_type &left = g->get_nonterminal(left_id);
    const nonterminal_type &right = g->get_nonterminal(right_id);
    const std::uint64_t left_exp_len =
      left.write_expansion(left_id, text, g);
    const std::uint64_t right_exp_len =
      right.write_expansion(right_id, text + left_exp_len, g);
    return left_exp_len + right_exp_len;
  }
}

//=============================================================================
// Compare the expansion of the nonterminal to the given text.
//=============================================================================
template<typename char_type, typename text_offset_type>
bool nonterminal<char_type, text_offset_type>::compare_expansion_to_text(
    const std::uint64_t id,
    const char_type * const text,
    const avl_grammar_multiroot<char_type, text_offset_type> * const g) const {
  typedef nonterminal<char_type, text_offset_type> nonterminal_type;
  const std::uint64_t height = g->get_height(id);
  if (height == 0) {
    const char_type my_char = g->get_char(id);
    return (text[0] == my_char);
  } else {
    const std::uint64_t left_id = g->get_left_id(id);
    const std::uint64_t right_id = g->get_right_id(id);
    const std::uint64_t left_exp_len = g->get_exp_len(left_id);
    const nonterminal_type &left = g->get_nonterminal(left_id);
    const nonterminal_type &right = g->get_nonterminal(right_id);
    if (!left.compare_expansion_to_text(left_id, text, g))
      return false;
    if (!right.compare_expansion_to_text(right_id, text + left_exp_len, g))
      return false;
    return true;
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
void nonterminal<char_type, text_offset_type>::decomposition(
    const std::uint64_t id,
    const std::uint64_t begin,
    const std::uint64_t end,
    space_efficient_vector<pair_type> &ret,
    const avl_grammar_multiroot<char_type, text_offset_type> * const g) const {

  // Handle boundary case.
  if (begin == end)
    return;

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
    const std::uint64_t exp_len = cur_range_end - cur_range_beg;
    ret.push_back(pair_type(
          (text_offset_type)x,
          (text_offset_type)exp_len));
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
          ret.push_back(pair_type(
                (text_offset_type)y,
                (text_offset_type)y_exp_len));
          suffix_length -= y_exp_len;
        } else {
          const std::uint64_t y_right_exp_len = g->get_exp_len(y_right);
          if (suffix_length > y_right_exp_len) {
            ret.push_back(pair_type(
                  (text_offset_type)y_right,
                  (text_offset_type)y_right_exp_len));
            suffix_length -= y_right_exp_len;
            y = y_left;
          } else y = y_right;
        }
      }
    }

    // Reverse the first sequence of nonterminals
    // collected during the left downward traversal.
    ret.reverse();

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
          ret.push_back(pair_type(
                (text_offset_type)y,
                (text_offset_type)y_exp_len));
          prefix_length -= y_exp_len;
        } else {
          const std::uint64_t y_left_exp_len = g->get_exp_len(y_left);
          if (prefix_length > y_left_exp_len) {
            ret.push_back(pair_type(
                  (text_offset_type)y_left,
                  (text_offset_type)y_left_exp_len));
            prefix_length -= y_left_exp_len;
            y = y_right;
          } else y = y_left;
        }
      }
    }
  }
}

#endif  // __AVL_GRAMMAR_MULTIROOT_HPP_INCLUDED
