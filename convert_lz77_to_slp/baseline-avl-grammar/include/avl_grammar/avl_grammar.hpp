#ifndef __AVL_GRAMMAR_HPP_INCLUDED
#define __AVL_GRAMMAR_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "../utils/hash_table.hpp"
#include "../utils/karp_rabin_hashing.hpp"


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
    typedef avl_grammar<char_type, text_offset_type> grammar_type;

    //=========================================================================
    // Class members.
    //=========================================================================
    const std::uint8_t m_height;
    const text_offset_type m_exp_len;
    const std::uint64_t m_kr_hash;
    const nonterminal_type * m_left_p;
    const nonterminal_type * m_right_p;

  public:

    //=========================================================================
    // Constructors.
    //=========================================================================
    nonterminal();
    nonterminal(const char_type);
    nonterminal(const std::uint8_t, const text_offset_type, const std::uint64_t,
      const nonterminal_type * const, const nonterminal_type * const);
    nonterminal(const nonterminal_type * const, const nonterminal_type * const);

    //=========================================================================
    // Access methods.
    //=========================================================================
    std::uint64_t get_height() const;
    std::uint64_t get_exp_len() const;
    std::uint64_t get_kr_hash() const;
    char_type get_char() const;
    const nonterminal_type * get_left_p() const;
    const nonterminal_type * get_right_p() const;

    //=========================================================================
    // Key methods.
    //=========================================================================
    void decomposition(const std::uint64_t, const std::uint64_t,
        std::vector<const nonterminal_type*> &) const;
    std::uint64_t write_expansion(char_type * const) const;

    //=========================================================================
    // Mostly unused.
    //=========================================================================
    void print_expansion() const;
    bool compare_expansion_to_text(const char_type * const) const;
    bool test_avl_property() const;
    std::uint64_t collect_mersenne_karp_rabin_hashes(
        std::vector<std::uint64_t> &) const;
    void collect_nonterminal_pointers(
        std::vector<const nonterminal_type*> &) const;
    std::uint64_t collect_mersenne_karp_rabin_hashes_2(
        hash_table<const nonterminal_type*, std::uint64_t> &) const;
    void count_nonterminals_in_pruned_grammar(
        hash_table<const nonterminal_type*, std::uint64_t> &h,
        hash_table<std::uint64_t, bool> &,
        std::uint64_t &) const;
} __attribute__((packed));

//=============================================================================
// Hash functions of the appropriate type.
// Used in the hash table used to prune the grammar.
//=============================================================================
typedef const nonterminal<std::uint8_t, uint40>* get_hash_ptr_type;

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
    // Get the root.
    //=========================================================================
    const nonterminal_type * get_root() const {
      return m_root;
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
    // Add a new binary nonterminal.
    //=========================================================================
    const nonterminal_type* add_nonterminal(
        const nonterminal_type *left_p,
        const nonterminal_type *right_p) {
      typedef const nonterminal_type * ptr_type;
      const nonterminal_type &left = *left_p;
      const nonterminal_type &right = *right_p;

      // Compute values for the new nonterminal.
      const std::uint64_t left_exp_len = left.get_exp_len();
      const std::uint64_t right_exp_len = right.get_exp_len();
      const std::uint64_t exp_len = left_exp_len + right_exp_len;
      const std::uint8_t left_height = left.get_height();
      const std::uint8_t right_height = right.get_height();
      const std::uint8_t height = std::max(left_height, right_height) + 1;
      const std::uint64_t left_kr_hash = left.get_kr_hash();
      const std::uint64_t right_kr_hash = right.get_kr_hash();
      const std::uint64_t new_kr_hash = karp_rabin_hashing::concat(
          left_kr_hash, right_kr_hash, right_exp_len);

      // Create and add new nonterminal.
      const ptr_type new_nonterm_p = new nonterminal_type(
         (std::uint8_t)height,
         (text_offset_type)exp_len,
         new_kr_hash, left_p, right_p);
      m_nonterminals.push_back(new_nonterm_p);

      // Return the ptr to the new nonterminal.
      return new_nonterm_p;
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
    // Check if the grammar expands to the given string.
    //=========================================================================
    bool compare_expansion_to_text(
        const char_type * const text,
        const std::uint64_t text_length) {

      // Compute length of expansion.
      std::uint64_t exp_total_len = m_root->get_exp_len();

      // If they differ, return false.
      if (text_length != exp_total_len)
        return false;

      // Otherwise, compare the generated string with `text'.
      return m_root->compare_expansion_to_text(text);
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
      std::vector<const nonterminal_type*> v;
      x->decomposition(begin, end, v);
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
        const nonterminal_type * const left_p,
        const nonterminal_type * const right_p) {
      typedef const nonterminal_type * ptr_type;

      // Consider two cases, depending on whether
      // left of right nonterminal is taller.
      const nonterminal_type &left = *left_p;
      const nonterminal_type &right = *right_p;
      if (left.get_height() >= right.get_height()) {
        if (left.get_height() - right.get_height() <= 1) {

          // Height are close. Just merge and return.
          const ptr_type newroot_p = add_nonterminal(left_p, right_p);
          return newroot_p;
        } else {
          const ptr_type leftleft_p = left.get_left_p();
          const ptr_type leftright_p = left.get_right_p();
          const nonterminal_type &leftleft = *leftleft_p;
          const ptr_type newright_p =
            add_concat_nonterminal(leftright_p, right_p);
          const nonterminal_type &newright = *newright_p;
          if (newright.get_height() > leftleft.get_height() &&
              newright.get_height() - leftleft.get_height() > 1) {

            // Rebalancing needed.
            const ptr_type newright_left_p = newright.get_left_p();
            const ptr_type newright_right_p = newright.get_right_p();
            const nonterminal_type &newright_left = *newright_left_p;
            const nonterminal_type &newright_right = *newright_right_p;
            if (newright_left.get_height() > newright_right.get_height()) {

              // Double (right-left) rotation.
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
          const nonterminal_type &rightright = *rightright_p;
          const ptr_type newleft_p =
            add_concat_nonterminal(left_p, rightleft_p);
          const nonterminal_type &newleft = *newleft_p;
          if (newleft.get_height() > rightright.get_height() &&
              newleft.get_height() - rightright.get_height() > 1) {

            // Rebalancing needed.
            const ptr_type newleft_left_p = newleft.get_left_p();
            const ptr_type newleft_right_p = newleft.get_right_p();
            const nonterminal_type &newleft_left = *newleft_left_p;
            const nonterminal_type &newleft_right = *newleft_right_p;
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
    // Merge greedily (shortest first) sequence of nonterminals.
    //=========================================================================
    const nonterminal_type* greedy_merge(
        std::vector<const nonterminal_type*> &seq) {
      while (seq.size() > 1) {

        // Find the nonterminal with the smallest height.
        std::uint64_t smallest_height_id = 0;
        for (std::uint64_t i = 1; i < seq.size(); ++i) {
          if (seq[i]->get_height() < seq[smallest_height_id]->get_height())
            smallest_height_id = i;
        }

        // Merge the nonterminal with the smaller height with
        // one of its beighbors (whichever is shorter).
        if (smallest_height_id == 0 ||
            (smallest_height_id + 1 < seq.size() &&
             seq[smallest_height_id + 1]->get_height() <=
             seq[smallest_height_id - 1]->get_height())) {

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

//=============================================================================
// Default constructor.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal()
  : m_height(0),
    m_exp_len(1),
    m_kr_hash(0),
    m_left_p(NULL),
    m_right_p(NULL) {}

//=============================================================================
// Constructor for a nonterminal expanding to a single symbol.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(const char_type c)
  : m_height(0),
    m_exp_len(1),
    m_kr_hash(karp_rabin_hashing::hash_char(c)),
    m_left_p((nonterminal_type *)((std::uint64_t)c)),
    m_right_p(NULL) {}

//=============================================================================
// Constructor for non-single-symbol nonterminal.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(
      const std::uint8_t height,
      const text_offset_type exp_len,
      const std::uint64_t kr_hash,
      const nonterminal_type * const left_p,
      const nonterminal_type * const right_p)
  : m_height(height),
    m_exp_len(exp_len),
    m_kr_hash(kr_hash),
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
// Get nonterminal expansion length.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>::get_kr_hash() const {
  return (std::uint64_t)m_kr_hash;
}

//=============================================================================
// Get nonterminal left ptr.
//=============================================================================
template<typename char_type, typename text_offset_type>
const nonterminal<char_type, text_offset_type>* nonterminal<char_type, text_offset_type>
::get_left_p() const {
  return m_left_p;
}

//=============================================================================
// Get nonterminal left ptr.
//=============================================================================
template<typename char_type, typename text_offset_type>
const nonterminal<char_type, text_offset_type>* nonterminal<char_type, text_offset_type>
::get_right_p() const {
  return m_right_p;
}

//=============================================================================
// Print the string encoded by the grammar.
//=============================================================================
template<typename char_type, typename text_offset_type>
void nonterminal<char_type, text_offset_type>::print_expansion() const {
  typedef const nonterminal_type * ptr_type;
  const nonterminal_type &x = *this;
  const std::uint64_t height = x.get_height();
  if (height == 0) {
    const char_type my_char = x.get_char();
    fprintf(stderr, "%c", (char)my_char);
  } else {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = *x_left_p;
    const nonterminal_type &x_right = *x_right_p;
    x_left.print_expansion();
    x_right.print_expansion();
  }
}

//=============================================================================
// Write the expansion into the given array.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>
::write_expansion(char_type * const text) const {
  typedef const nonterminal_type * ptr_type;
  const ptr_type x_p = this;
  const nonterminal_type &x = *x_p;
  const std::uint64_t x_height = x.get_height();
  if (x_height == 0) {
    const char_type my_char = x.get_char();
    text[0] = my_char;
    return 1;
  } else {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = *x_left_p;
    const nonterminal_type &x_right = *x_right_p;
    const std::uint64_t x_left_exp_len =
      x_left.write_expansion(text);
    const std::uint64_t x_right_exp_len =
      x_right.write_expansion(text + x_left_exp_len);
    return x_left_exp_len + x_right_exp_len;
  }
}

//=============================================================================
// Compare the expansion of the nonterminal to the given text.
//=============================================================================
template<typename char_type, typename text_offset_type>
bool nonterminal<char_type, text_offset_type>::compare_expansion_to_text(
    const char_type * const text) const {
  typedef const nonterminal_type * ptr_type;
  const ptr_type x_p = this;
  const nonterminal_type &x = *x_p;
  const std::uint64_t x_height = x.get_height();
  if (x_height == 0) {
    const char_type my_char = x.get_char();
    return (text[0] == my_char);
  } else {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = *x_left_p;
    const nonterminal_type &x_right = *x_right_p;
    const std::uint64_t x_left_exp_len = x_left.get_exp_len();
    if (!x_left.compare_expansion_to_text(text))
      return false;
    if (!x_right.compare_expansion_to_text(text + x_left_exp_len))
      return false;
    return true;
  }
}

//=============================================================================
// Test the AVL propert of a subtree.
//=============================================================================
template<typename char_type, typename text_offset_type>
bool nonterminal<char_type, text_offset_type>::test_avl_property() const {
  typedef const nonterminal_type * ptr_type;
  const ptr_type x_p = this;
  const nonterminal_type &x = *x_p;
  const std::uint64_t x_height = x.get_height();
  if (x_height == 0)
    return true;

  const ptr_type x_left_p = x.get_left_p();
  const ptr_type x_right_p = x.get_right_p();
  const nonterminal_type &x_left = *x_left_p;
  const nonterminal_type &x_right = *x_right_p;
  if (!x_left.test_avl_property() ||
      !x_right.test_avl_property())
    return false;

  const std::uint64_t x_left_height = x_left.get_height();
  const std::uint64_t x_right_height = x_right.get_height();
  if ((x_right_height > x_left_height && x_right_height - x_left_height > 1) ||
      (x_right_height < x_left_height && x_left_height - x_right_height > 1))
    return false;

  return true;
}

//=============================================================================
// Collect Mersenne Karp-Rabin hashes of all nonterminals.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>
::collect_mersenne_karp_rabin_hashes(
    std::vector<std::uint64_t> &hashes) const {
  typedef const nonterminal_type * ptr_type;
  const ptr_type x_p = this;
  const nonterminal_type &x = *x_p;
  const std::uint64_t x_height = x.get_height();
  if (x_height == 0) {
    const char_type my_char = x.get_char();
    const std::uint64_t h = karp_rabin_hashing::hash_char(my_char);
    hashes.push_back(h);
    return h;
  } else {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = *x_left_p;
    const nonterminal_type &x_right = *x_right_p;
    const std::uint64_t x_right_len = x_right.get_exp_len();
    const std::uint64_t x_left_hash =
      x_left.collect_mersenne_karp_rabin_hashes(hashes);
    const std::uint64_t x_right_hash =
      x_right.collect_mersenne_karp_rabin_hashes(hashes);
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
    std::vector<const nonterminal_type*> &pointers) const {
  typedef const nonterminal_type * ptr_type;
  const ptr_type x_p = this;
  const nonterminal_type &x = *x_p;
  const std::uint64_t x_height = x.get_height();
  pointers.push_back(x_p);
  if (x_height > 0) {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = *x_left_p;
    const nonterminal_type &x_right = *x_right_p;
    x_left.collect_nonterminal_pointers(pointers);
    x_right.collect_nonterminal_pointers(pointers);
  }
}

//=============================================================================
// Collect Mersenne Karp-Rabin hashes of all nonterminals.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>
::collect_mersenne_karp_rabin_hashes_2(
    hash_table<const nonterminal_type*, std::uint64_t> &hashes) const {
  typedef const nonterminal_type * ptr_type;
  const ptr_type x_p = this;
  const nonterminal_type &x = *x_p;
  const std::uint64_t x_height = x.get_height();
  if (x_height == 0) {
    const char_type my_char = x.get_char();
    const std::uint64_t h = karp_rabin_hashing::hash_char(my_char);
    hashes.insert(x_p, h);
    return h;
  } else {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = *x_left_p;
    const nonterminal_type &x_right = *x_right_p;
    const std::uint64_t x_right_len = x_right.get_exp_len();
    const std::uint64_t x_left_hash =
      x_left.collect_mersenne_karp_rabin_hashes_2(hashes);
    const std::uint64_t x_right_hash =
      x_right.collect_mersenne_karp_rabin_hashes_2(hashes);
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
void nonterminal<char_type, text_offset_type>::
count_nonterminals_in_pruned_grammar(
    hash_table<const nonterminal_type*, std::uint64_t> &hashes,
    hash_table<std::uint64_t, bool> &seen_hashes,
    std::uint64_t &current_count) const {
  typedef const nonterminal_type * ptr_type;
  const std::uint64_t * const h = hashes.find(this);
  if (seen_hashes.find(*h) == NULL) {
    seen_hashes.insert(*h, true);
    ++current_count;
    const nonterminal_type &x = *this;
    const std::uint64_t x_height = x.get_height();
    if (x_height != 0) {
      const ptr_type x_left_p = x.get_left_p();
      const ptr_type x_right_p = get_right_p();
      const nonterminal_type &x_left = *x_left_p;
      const nonterminal_type &x_right = *x_right_p;
      x_left.count_nonterminals_in_pruned_grammar(
          hashes, seen_hashes, current_count);
      x_right.count_nonterminals_in_pruned_grammar(
          hashes, seen_hashes, current_count);
    }
  }
}

//=============================================================================
// Assuming S is the expansions of the nontermnal, return the
// sequence of nonterminals expanding to S[begin..end).
//=============================================================================
template<typename char_type, typename text_offset_type>
void nonterminal<char_type, text_offset_type>::decomposition(
    const std::uint64_t begin,
    const std::uint64_t end,
    std::vector<const nonterminal_type*> &ret) const {
  typedef const nonterminal_type * ptr_type;
  ptr_type x_p = this;

  // Handle boundary case.
  if (begin == end)
    return;

  // Compute height and expansion length for x.
  std::uint64_t x_height = 0;
  std::uint64_t x_exp_len = 0;
  {
    const nonterminal_type &x = *x_p;
    x_height = x.get_height();
    x_exp_len = x.get_exp_len();
  }

  // Find the deepest nonterminal in the parse tree
  // containing the range [begin..end).
  std::uint64_t cur_range_beg = 0;
  std::uint64_t cur_range_end = x_exp_len;
  while (x_height > 0) {
    const nonterminal_type &x = *x_p;
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = *x_left_p;
    const nonterminal_type &x_right = *x_right_p;
    const std::uint64_t x_left_exp_len = x_left.get_exp_len();
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
      const nonterminal_type &x = *x_p;
      const ptr_type x_left_p = x.get_left_p();
      const nonterminal_type &x_left = *x_left_p;
      const std::uint64_t x_left_exp_len = x_left.get_exp_len();
      const std::uint64_t left_range_end = cur_range_beg + x_left_exp_len;
      std::uint64_t suffix_length = left_range_end - begin;
      ptr_type y_p = x_left_p;
      while (suffix_length > 0) {
        const nonterminal_type &y = *y_p;
        const std::uint64_t y_exp_len = y.get_exp_len();
        const ptr_type y_left_p = y.get_left_p();
        const ptr_type y_right_p = y.get_right_p();
        if (y_exp_len == suffix_length) {
          ret.push_back(y_p);
          suffix_length -= y_exp_len;
        } else {
          const nonterminal_type &y_right = *y_right_p;
          const std::uint64_t y_right_exp_len = y_right.get_exp_len();
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
    std::reverse(ret.begin(), ret.end());

    // Perform the analogous operation for the right side.
    {
      const nonterminal_type &x = *x_p;
      const ptr_type x_left_p = x.get_left_p();
      const ptr_type x_right_p = x.get_right_p();
      const nonterminal_type &x_left = *x_left_p;
      const std::uint64_t x_left_exp_len = x_left.get_exp_len();
      const std::uint64_t right_range_beg = cur_range_beg + x_left_exp_len;
      std::uint64_t prefix_length = end - right_range_beg;
      ptr_type y_p = x_right_p;
      while (prefix_length > 0) {
        const nonterminal_type &y = *y_p;
        const std::uint64_t y_exp_len = y.get_exp_len();
        const ptr_type y_left_p = y.get_left_p();
        const ptr_type y_right_p = y.get_right_p();
        if (y_exp_len == prefix_length) {
          ret.push_back(y_p);
          prefix_length -= y_exp_len;
        } else {
          const nonterminal_type &y_left = *y_left_p;
          const std::uint64_t y_left_exp_len = y_left.get_exp_len();
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
