#ifndef __AVL_GRAMMAR_HPP_INCLUDED
#define __AVL_GRAMMAR_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "../utils/hash_table.hpp"
#include "../utils/karp_rabin_hashing.hpp"
#include "../utils/space_efficient_vector.hpp"


//=============================================================================
// Class used to represent multiroot AVL grammar. Forward declaration.
//=============================================================================
template<typename char_type, typename text_offset_type>
struct avl_grammar;

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
    typedef avl_grammar<char_type, text_offset_type> grammar_type;

    //=========================================================================
    // Class members.
    //=========================================================================
    std::uint8_t m_height;
    text_offset_type m_exp_len;
    ptr_type m_left_p;
    ptr_type m_right_p;

  public:

    //=========================================================================
    // Constructors.
    //=========================================================================
    nonterminal();
    nonterminal(const char_type);
    nonterminal(const std::uint8_t, const text_offset_type,
      const ptr_type, const ptr_type);

    //=========================================================================
    // Access methods.
    //=========================================================================
    std::uint64_t get_height() const;
    std::uint64_t get_exp_len() const;
    char_type get_char() const;
    ptr_type get_left_p() const;
    ptr_type get_right_p() const;

    //=========================================================================
    // Key methods.
    //=========================================================================
    void decomposition(const ptr_type, const std::uint64_t,
        const std::uint64_t, space_efficient_vector<ptr_type> &,
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
// Implementation of the avl_grammar class.
//=============================================================================
template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
struct avl_grammar {
  static_assert(sizeof(char_type) <= sizeof(text_offset_type),
      "Error: sizeof(char_type) > sizeof(text_offset_type)!");

  private:

    //=========================================================================
    // Declare typedefs.
    //=========================================================================
    typedef nonterminal<char_type, text_offset_type> nonterminal_type;
    typedef text_offset_type ptr_type;

    //=========================================================================
    // Class members.
    //=========================================================================
    space_efficient_vector<nonterminal_type> m_nonterminals;
    ptr_type m_root;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    avl_grammar() :
      m_root(std::numeric_limits<ptr_type>::max()) {}

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~avl_grammar() {
    }

    //=========================================================================
    // Print the string encoded by the grammar.
    //=========================================================================
    void print_expansion() const {
      const nonterminal_type &root = get_nonterminal(m_root);
      root.print_expansion(m_root, this);
    }

    //=========================================================================
    // Return the number of nonterminals.
    //=========================================================================
    std::uint64_t size() const {
      return m_nonterminals.size();
    }

    //=========================================================================
    // Get the root.
    //=========================================================================
    ptr_type get_root() const {
      return m_root;
    }

    //=========================================================================
    // Set the root.
    //=========================================================================
    void set_root(
        const ptr_type id) {
      m_root = id;
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
      const std::uint64_t exp_len = nonterm.get_exp_len();
      return exp_len;
    }

    //=========================================================================
    // Add a nonterminal.
    //=========================================================================
    ptr_type add_nonterminal(const nonterminal_type &nonterm) {
      const ptr_type new_nonterm_p = m_nonterminals.size();
      m_nonterminals.push_back(nonterm);

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
      const std::uint8_t left_height = left.get_height();
      const std::uint8_t right_height = right.get_height();
      const std::uint8_t height = std::max(left_height, right_height) + 1;

      // Create and add new nonterminal.
      const ptr_type new_nonterm_p = m_nonterminals.size();
      nonterminal_type new_nonterm(height, exp_len, left_p, right_p);
      m_nonterminals.push_back(new_nonterm);

      // Return the ptr to the new nonterminal.
      return new_nonterm_p;
    }

    //=========================================================================
    // Decode the text and write to a given array.
    //=========================================================================
    void decode(
        char_type* &text,
        std::uint64_t &text_length) {
      text_length = get_exp_len(m_root);
      text = new char_type[text_length];
      const nonterminal_type &nonterm = get_nonterminal(m_root);
      nonterm.write_expansion(m_root, text, this);
    }

    //=========================================================================
    // Check if the grammar expands to the given string.
    //=========================================================================
    bool compare_expansion_to_text(
        const char_type * const text,
        const std::uint64_t text_length) {

      // Compute length of expansion.
      std::uint64_t exp_total_len = get_exp_len(m_root);

      // If they differ, return false.
      if (text_length != exp_total_len)
        return false;

      // Otherwise, compare the generated string with `text'.
      const nonterminal_type &nonterm = get_nonterminal(m_root);
      return nonterm.compare_expansion_to_text(m_root, text, this);
    }

    //=========================================================================
    // Test the AVL property of all nonterminals.
    //=========================================================================
    bool test_avl_property() {
      const nonterminal_type &nonterm = get_nonterminal(m_root);
      return nonterm.test_avl_property(m_root, this);
    }

    //=========================================================================
    // Collect Mersenne Karp-Rabin hashes in a vector.
    //=========================================================================
    void collect_mersenne_karp_rabin_hashes(
        std::vector<std::uint64_t> &hashes) {
      const nonterminal_type &nonterm = get_nonterminal(m_root);
      (void) nonterm.collect_mersenne_karp_rabin_hashes(m_root, hashes, this);
    }

    //=========================================================================
    // Collect Mersenne Karp-Rabin hashes in a hash table.
    //=========================================================================
    void collect_mersenne_karp_rabin_hashes_2(
        hash_table<ptr_type, std::uint64_t> &hashes) {
      const nonterminal_type &nonterm = get_nonterminal(m_root);
      (void) nonterm.collect_mersenne_karp_rabin_hashes_2(
        m_root, hashes, this);
    }

    //=========================================================================
    // Count nonterminals in the pruned grammar.
    //=========================================================================
    void count_nonterminals_in_pruned_grammar(
        hash_table<ptr_type, std::uint64_t> &hashes,
        hash_table<std::uint64_t, bool> &seen_hashes,
        std::uint64_t &current_count) {
      const nonterminal_type &nonterm = get_nonterminal(m_root);
      nonterm.count_nonterminals_in_pruned_grammar(
          m_root, hashes, seen_hashes, current_count, this);
    }

    //=========================================================================
    // Collect pointers to all nonterminals reachable from the root.
    //=========================================================================
    void collect_nonterminal_pointers(
        std::vector<ptr_type> &pointers) {
      const nonterminal_type &nonterm = get_nonterminal(m_root);
      nonterm.collect_nonterminal_pointers(m_root, pointers, this);
    }

    //=========================================================================
    // Add a nonterminal expanding to a substring of a given nonterminal.
    //=========================================================================
    ptr_type add_substring_nonterminal(
        const ptr_type x_p,
        const std::uint64_t begin,
        const std::uint64_t end) {
      space_efficient_vector<ptr_type> v;
      const nonterminal_type &nonterm = get_nonterminal(x_p);
      nonterm.decomposition(x_p, begin, end, v, this);
      return greedy_merge(v);
    }

    //=========================================================================
    // Add a substring expanding to a substring of grammar.
    //=========================================================================
    ptr_type add_substring_nonterminal(
        const std::uint64_t begin,
        const std::uint64_t end) {
      return add_substring_nonterminal(m_root, begin, end);
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

  private:

    //=========================================================================
    // Heap down routine.
    //========================================================================
    void heap_down(
        std::uint64_t i,
        const space_efficient_vector<ptr_type> &seq,
        text_offset_type * const heap,
        const std::uint64_t heap_size) const {
      ++i;
      std::uint64_t min_pos = i;
      std::uint64_t height = 0;
      const std::uint64_t val = heap[min_pos - 1];
      {
        const ptr_type nonterm_p = seq[val];
        const nonterminal_type &nonterm = get_nonterminal(nonterm_p);
        height = nonterm.get_height();
      }
      while (true) {
        std::uint64_t min_height = height;
        std::uint64_t min_val = val;
        if ((i << 1) <= heap_size) {
          const std::uint64_t left_val = heap[(i << 1) - 1];
          const ptr_type left_p = seq[left_val];
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
          const ptr_type right_p = seq[right_val];
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
        const space_efficient_vector<ptr_type> &seq,
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
        const space_efficient_vector<ptr_type> &seq,
        text_offset_type * const heap,
        const std::uint64_t heap_size) const {
      for (std::uint64_t i = heap_size / 2; i > 0; --i)
        heap_down(i - 1, seq, heap, heap_size);
    }

    //=========================================================================
    // Merge greedily (shortest first) sequence of nonterminals.
    // Uses binary heap to achieve O(m log m) time.
    //=========================================================================
    ptr_type greedy_merge(space_efficient_vector<ptr_type> &seq) {

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
          const ptr_type prev_nonterm_p = seq[idx];
          const nonterminal_type &prev_nonterm =
            get_nonterminal(prev_nonterm_p);
          prev_height = prev_nonterm.get_height();
        }
        if ((std::uint64_t)next[min_elem] != sentinel) {
          const std::uint64_t idx = next[min_elem];
          const ptr_type next_nonterm_p = seq[idx];
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
          ret = seq[min_elem];
          break;
        } else if ((std::uint64_t)prev[min_elem] == sentinel ||
            ((std::uint64_t)next[min_elem] != sentinel &&
             next_height <= prev_height)) {

          // Only right neighbor exists, or both exist
          // and the right one is not taller than the left
          // one. End result: merge with the right neighbor.
          const std::uint64_t right_elem = next[min_elem];
          const ptr_type left_p = seq[min_elem];
          const ptr_type right_p = seq[right_elem];
          seq[min_elem] = add_concat_nonterminal(left_p, right_p);
          deleted[right_elem] = true;
          next[min_elem] = next[right_elem];
          prev[next[min_elem]] = min_elem;
          heap_down(0, seq, heap, heap_size);
        } else {

          // Only left neighbor exists, or both exists
          // and the left one is not taller than the
          // right one. End result: merge with left neighbor.
          const std::uint64_t left_elem = prev[min_elem];
          const ptr_type left_p = seq[left_elem];
          const ptr_type right_p = seq[min_elem];
          seq[min_elem] = add_concat_nonterminal(left_p, right_p);
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
    m_left_p(std::numeric_limits<text_offset_type>::max()),
    m_right_p(std::numeric_limits<text_offset_type>::max()) {}

//=============================================================================
// Constructor for a nonterminal expanding to a single symbol.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(const char_type c)
  : m_height(0),
    m_exp_len(1),
    m_left_p((text_offset_type)((std::uint64_t)c)),
    m_right_p(std::numeric_limits<text_offset_type>::max()) {}

//=============================================================================
// Constructor for non-single-symbol nonterminal.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(
      const std::uint8_t height,
      const text_offset_type exp_len,
      const ptr_type left_p,
      const ptr_type right_p)
  : m_height(height),
    m_exp_len(exp_len),
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
std::uint64_t nonterminal<char_type, text_offset_type>::get_exp_len() const {
  return (std::uint64_t)m_exp_len;
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
    space_efficient_vector<ptr_type> &ret,
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
    ret.push_back(x_p);
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
          ret.push_back(y_p);
          suffix_length -= y_exp_len;
        } else {
          const std::uint64_t y_right_exp_len = g->get_exp_len(y_right_p);
          if (suffix_length > y_right_exp_len) {
            ret.push_back(y_right_p);
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
          ret.push_back(y_p);
          prefix_length -= y_exp_len;
        } else {
          const std::uint64_t y_left_exp_len = g->get_exp_len(y_left_p);
          if (prefix_length > y_left_exp_len) {
            ret.push_back(y_left_p);
            prefix_length -= y_left_exp_len;
            y_p = y_right_p;
          } else y_p = y_left_p;
        }
      }
    }
  }
}

#endif  // __AVL_GRAMMAR_HPP_INCLUDED
