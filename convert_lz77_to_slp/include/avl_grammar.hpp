#ifndef __AVL_GRAMMAR_HPP_INCLUDED
#define __AVL_GRAMMAR_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "hash_table.hpp"
#include "karp_rabin_hashing.hpp"


//=============================================================================
// A class used to represent the nonterminal in the AVL grammar.
//=============================================================================
template<typename char_type>
struct avl_grammar_node {
  typedef avl_grammar_node<char_type> node_type;

  // Class members.
  const char_type m_char;
  const std::uint8_t m_height;
  const std::uint64_t m_exp_len;
  const node_type * const m_left;
  const node_type * const m_right;

  // Default constructor.
  avl_grammar_node() :
      m_char((char_type)0),
      m_height(0),
      m_exp_len(1),
      m_left(NULL),
      m_right(NULL) {}

  // Constructor for a node expanding to a single symbol.
  avl_grammar_node(const char_type c) :
      m_char(c),
      m_height(0),
      m_exp_len(1),
      m_left(NULL),
      m_right(NULL) {}

  // Constructor for a standard nonterminal
  // (expanding to two other nonterminals).
  avl_grammar_node(
      const node_type * const left,
      const node_type * const right) :
        m_char((char_type)0),
        m_height(std::max(left->m_height, right->m_height) + 1),
        m_exp_len(left->m_exp_len + right->m_exp_len),
        m_left(left),
        m_right(right) {}

  // Print the string encoded by the grammar.
  void print_expansion() const {
    if (m_height == 0)
      fprintf(stderr, "%c", (char)m_char);
    else {
      m_left->print_expansion();
      m_right->print_expansion();
    }
  }

  // Write the expansion into the given array.
  void write_expansion(char_type * const text) const {
    if (m_height == 0)
      text[0] = m_char;
    else {
      m_left->write_expansion(text);
      m_right->write_expansion(text + m_left->m_exp_len);
    }
  }

  // Test the AVL propert of a subtree.
  bool test_avl_property() const {
    if (m_height == 0)
      return true;
    else {
      if ((!(m_left->test_avl_property())) ||
          (!(m_right->test_avl_property())) ||
          (m_right->m_height > m_left->m_height &&
           m_right->m_height - m_left->m_height > 1) ||
          (m_left->m_height > m_right->m_height &&
           m_left->m_height - m_right->m_height > 1))
        return false;
      else return true;
    }
  }

  // Collect Karp-Rabin hashes of all nonterminals.
  std::uint64_t collect_karp_rabin_hashes(
      std::vector<std::uint64_t> &hashes,
      const std::uint64_t a,
      const std::uint64_t p) const {
    if (m_height == 0) {
      const std::uint64_t h = m_char % p;
      hashes.push_back(h);
      return h;
    } else {
      const std::uint64_t left_hash =
        m_left->collect_karp_rabin_hashes(hashes, a, p);
      const std::uint64_t right_hash =
        m_right->collect_karp_rabin_hashes(hashes, a, p);
      std::uint64_t h =
        (left_hash * mod_pow<std::uint64_t>(a, m_right->m_exp_len, p)) % p;
      h = (h + right_hash) % p;
      hashes.push_back(h);
      return h;
    }
  }

  // Collect Mersenne Karp-Rabin hashes of all nonterminals.
  std::uint64_t collect_mersenne_karp_rabin_hashes(
      std::vector<std::uint64_t> &hashes,
      const std::uint64_t hash_variable,
      const std::uint64_t mersenne_prime_exponent) const {
    if (m_height == 0) {
      const std::uint64_t h = mod_mersenne(m_char, mersenne_prime_exponent);
      hashes.push_back(h);
      return h;
    } else {
      const std::uint64_t left_hash =
        m_left->collect_mersenne_karp_rabin_hashes(hashes,
            hash_variable, mersenne_prime_exponent);
      const std::uint64_t right_hash =
        m_right->collect_mersenne_karp_rabin_hashes(hashes,
            hash_variable, mersenne_prime_exponent);
      std::uint64_t h = 0;
      {
        const std::uint64_t pow = pow_mod_mersenne(hash_variable,
            m_right->m_exp_len, mersenne_prime_exponent);
        const std::uint64_t tmp = mul_mod_meresenne(left_hash,
            pow, mersenne_prime_exponent);
        h = mod_mersenne(tmp + right_hash, mersenne_prime_exponent);
      }
      hashes.push_back(h);
      return h;
    }
  }

  // Collect reachable nonterminals.
  void collect_nonterminal_pointers(
      std::vector<const node_type*> &pointers) const {
    pointers.push_back(this);
    if (m_height > 0) {
      m_left->collect_nonterminal_pointers(pointers);
      m_right->collect_nonterminal_pointers(pointers);
    }
  }

  // Collect Mersenne Karp-Rabin hashes of all nonterminals.
  std::uint64_t collect_mersenne_karp_rabin_hashes_2(
      hash_table<const node_type*, std::uint64_t> &hashes,
      const std::uint64_t hash_variable,
      const std::uint64_t mersenne_prime_exponent) const {
    if (m_height == 0) {
      const std::uint64_t h = mod_mersenne(m_char, mersenne_prime_exponent);
      hashes.insert(this, h);
      return h;
    } else {
      const std::uint64_t left_hash =
        m_left->collect_mersenne_karp_rabin_hashes_2(hashes,
            hash_variable, mersenne_prime_exponent);
      const std::uint64_t right_hash =
        m_right->collect_mersenne_karp_rabin_hashes_2(hashes,
            hash_variable, mersenne_prime_exponent);
      std::uint64_t h = 0;
      {
        const std::uint64_t pow = pow_mod_mersenne(hash_variable,
            m_right->m_exp_len, mersenne_prime_exponent);
        const std::uint64_t tmp = mul_mod_meresenne(left_hash,
            pow, mersenne_prime_exponent);
        h = mod_mersenne(tmp + right_hash, mersenne_prime_exponent);
      }
      hashes.insert(this, h);
      return h;
    }
  }

  void count_nodes_in_pruned_grammar(
      hash_table<const node_type*, std::uint64_t> &hashes,
      hash_table<std::uint64_t, bool> &seen_hashes,
      std::uint64_t &current_count) const {
    if (m_height == 0) {
      const std::uint64_t *h = hashes.find(this);
      if (seen_hashes.find(*h) == NULL) {
        seen_hashes.insert(*h, true);
        ++current_count;
      }
    } else {
      const std::uint64_t *h = hashes.find(this);
      if (seen_hashes.find(*h) == NULL) {
        seen_hashes.insert(*h, true);
        ++current_count;
        m_left->count_nodes_in_pruned_grammar(hashes,
            seen_hashes, current_count);
        m_right->count_nodes_in_pruned_grammar(hashes,
            seen_hashes, current_count);
      }
    }
  }
};

//=============================================================================
// Hash functions of the appropriate type.
// Used in the hash table used to prune the grammar.
//=============================================================================
typedef const avl_grammar_node<std::uint8_t>* get_hash_ptr_type;

template<>
std::uint64_t get_hash(const get_hash_ptr_type &x) {
  return (std::uint64_t)x * (std::uint64_t)29996224275833;
}

template<>
std::uint64_t get_hash(const std::uint64_t &x) {
  return (std::uint64_t)x * (std::uint64_t)4972548694736365;
}

//=============================================================================
// A class storing AVL grammar.
//=============================================================================
template<typename char_type>
struct avl_grammar {
  typedef avl_grammar_node<char_type> node_type;

  // Class members.
  std::vector<const node_type*> m_nonterminals;
  const node_type *m_root;

  // Constructor.
  avl_grammar() :
    m_root(NULL) {}

  // Print the string encoded by the grammar.
  void print_expansion() const {
    m_root->print_expansion();
  }

  // Return the number of nonterminals.
  std::uint64_t size() const {
    return m_nonterminals.size();
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

  // Collect Karp-Rabin hashes in a vector.
  void collect_karp_rabin_hashes(
      std::vector<std::uint64_t> &hashes,
      const std::uint64_t a = (std::uint64_t)999285268,
      const std::uint64_t p = (std::uint64_t)1000000009) const {
    (void) m_root->collect_karp_rabin_hashes(hashes, a, p);
  }

  // Collect Mersenne Karp-Rabin hashes in a vector.
  // Allows specifying variable and prime exponent.
  void collect_mersenne_karp_rabin_hashes(
      std::vector<std::uint64_t> &hashes,
      const std::uint64_t hash_variable,
      const std::uint64_t mersenne_prime_exponent) const {
    (void) m_root->collect_mersenne_karp_rabin_hashes(hashes,
        hash_variable, mersenne_prime_exponent);
  }

  // Collect Mersenne Karp-Rabin hashes in a hash table.
  // Allows specifying variable and prime exponent.
  void collect_mersenne_karp_rabin_hashes_2(
      hash_table<const node_type*, std::uint64_t> &hashes,
      const std::uint64_t hash_variable,
      const std::uint64_t mersenne_prime_exponent) const {
    (void) m_root->collect_mersenne_karp_rabin_hashes_2(hashes,
        hash_variable, mersenne_prime_exponent);
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
    m_root->count_nodes_in_pruned_grammar(hashes,
        seen_hashes, current_count);
  }

  // Collect pointers to all nonterminals reachable from the root.
  void collect_nonterminal_pointers(
      std::vector<const node_type*> &pointers) const {
    m_root->collect_nonterminal_pointers(pointers);
  }
};

//=============================================================================
// Given two nonterminals `left' and `right' of the same grammar
// expanding respectively to strings X and Y, add to grammar a
// nonterminals that expands to XY, and return the pointer to it.
// The root of the original grammar remains unchanged.
//=============================================================================
template<typename char_type>
const avl_grammar_node<char_type> *add_concat_nonterminal(
    std::vector<const avl_grammar_node<char_type> *> &nonterminals,
    const avl_grammar_node<char_type> * const left,
    const avl_grammar_node<char_type> * const right) {
  
  // Declare type.
  typedef avl_grammar_node<char_type> node_type;

  // Consider two cases, depending on whether
  // left of right nonterminal is taller.
  if (left->m_height >= right->m_height) {
    if (left->m_height - right->m_height <= 1) {

      // Height are close. Just merge and return.
      const node_type * const newroot = new node_type(left, right);
      nonterminals.push_back(newroot);
      return newroot;
    } else {
      const node_type * const newright =
        add_concat_nonterminal<char_type>(nonterminals, left->m_right, right);
      if (newright->m_height > left->m_left->m_height &&
          newright->m_height - left->m_left->m_height > 1) {
        
        // Rebalancing needed.
        if (newright->m_left->m_height > newright->m_right->m_height) {

          // Double (right-left) rotation.
          const node_type * const X =
            new node_type(left->m_left, newright->m_left->m_left);
          const node_type * const Z =
            new node_type(newright->m_left->m_right, newright->m_right);
          const node_type * const Y = new node_type(X, Z);
          nonterminals.push_back(X);
          nonterminals.push_back(Y);
          nonterminals.push_back(Z);
          return Y;
        } else {

          // Single (left) rotation.
          const node_type * const X =
            new node_type(left->m_left, newright->m_left);
          const node_type * const Y = new node_type(X, newright->m_right);
          nonterminals.push_back(X);
          nonterminals.push_back(Y);
          return Y;
        }
      } else {

        // No need to rebalance.
        const node_type * const newroot =
          new node_type(left->m_left, newright);
        nonterminals.push_back(newroot);
        return newroot;
      }
    }
  } else {
    if (right->m_height - left->m_height <= 1) {

      // Heights are close. Just merge and return.
      const node_type * const newroot = new node_type(left, right);
      nonterminals.push_back(newroot);
      return newroot;
    } else {
      const node_type * const newleft =
        add_concat_nonterminal<char_type>(nonterminals, left, right->m_left);
      if (newleft->m_height > right->m_right->m_height &&
          newleft->m_height - right->m_right->m_height > 1) {

        // Rebalancing needed.
        if (newleft->m_right->m_height > newleft->m_left->m_height) {

          // Double (left-right) rotation.
          const node_type * const X =
            new node_type(newleft->m_left, newleft->m_right->m_left);
          const node_type * const Z =
            new node_type(newleft->m_right->m_right, right->m_right);
          const node_type * const Y = new node_type(X, Z);
          nonterminals.push_back(X);
          nonterminals.push_back(Y);
          nonterminals.push_back(Z);
          return Y;
        } else {

          // Single (right) rotation.
          const node_type * const Y =
            new node_type(newleft->m_right, right->m_right);
          const node_type * const X = new node_type(newleft->m_left, Y);
          nonterminals.push_back(X);
          nonterminals.push_back(Y);
          return X;
        }
      } else {

        // No need to rebalance.
        const node_type * const newroot =
          new node_type(newleft, right->m_right);
        nonterminals.push_back(newroot);
        return newroot;
      }
    }
  }
}

//=============================================================================
// Given the AVL grammar expanding to string T, create and return
// the nonterminal expanding to T[begin..end). The root of the
// input grammar remains unchanged.
//=============================================================================
template<typename char_type>
const avl_grammar_node<char_type> *add_substring_nonterminal(
    std::vector<const avl_grammar_node<char_type> *> &nonterminals,
    const avl_grammar_node<char_type> * const root,
    const std::uint64_t begin,
    const std::uint64_t end) {

  // Check input correctness.
  if (begin > end ||
      end > root->m_exp_len) {
    fprintf(stderr, "\nError: extract: end > root->m_exp_len!\n");
    std::exit(EXIT_FAILURE);
  }
  
  // Declare types.
  typedef avl_grammar_node<char_type> node_type;

  // Handle boundary case.
  if (begin == end)
    return (node_type *)NULL;

  // Find the deepest node in the parse tree containing the range
  // [begin..end).
  const node_type *x = root;
  std::uint64_t cur_range_beg = 0;
  std::uint64_t cur_range_end = x->m_exp_len;
  while (x->m_height > 0 &&
      (end <= cur_range_beg + x->m_left->m_exp_len ||
       begin >= cur_range_beg + x->m_left->m_exp_len)) {
    if (end <= cur_range_beg + x->m_left->m_exp_len) {
      cur_range_end = cur_range_beg + x->m_left->m_exp_len;
      x = x->m_left;
    } else {
      cur_range_beg += x->m_left->m_exp_len;
      x = x->m_right;
    }
  }

  // Check if the range of x is exactly [begin..end).
  if (cur_range_beg == begin && cur_range_end == end) {

    // If yes, return x as the answer.
    return x;
  } else {

    // Otherwise, we perform two traversals in the tree. If by G we
    // denote the AVL grammar x->m_left, we first create the AVL
    // grammar expanding to suffix of length (cur_range_beg +
    // x->m_left->m_exp_len - begin) of the string exp(G). The
    // root of the resulting grammar is stored in left_grammar.
    const node_type * left_grammar = NULL;
    {
      const std::uint64_t left_range_end =
        cur_range_beg + x->m_left->m_exp_len;
      std::uint64_t suffix_length = left_range_end - begin;
      const node_type *y = x->m_left;

      // The tricky thing here is that we cannot merge the nodes right
      // away. For the merge to be O(log n) time, we need to merge
      // grammar short-to-tall and hence we first collect their roots
      // in a vector.
      std::vector<const node_type*> grammars_to_merge;
      while (suffix_length > 0) {
        if (y->m_exp_len == suffix_length) {
          grammars_to_merge.push_back(y);
          suffix_length -= y->m_exp_len;
        } else if (suffix_length > y->m_right->m_exp_len) {
          grammars_to_merge.push_back(y->m_right);
          suffix_length -= y->m_right->m_exp_len;
          y = y->m_left;
        } else y = y->m_right;
      }

      // Merge AVL grammars in grammars_to_merge into a single AVL grammar.
      while (grammars_to_merge.size() > 1) {
        const node_type * const left = grammars_to_merge.back();
        grammars_to_merge.pop_back();
        const node_type * const right = grammars_to_merge.back();
        grammars_to_merge.pop_back();
        grammars_to_merge.push_back(
            add_concat_nonterminal<char_type>(nonterminals, left, right));
      }
      left_grammar = grammars_to_merge.back();
    }

    // Perform the analogous operation for the right side.
    const node_type *right_grammar = NULL;
    {
      const std::uint64_t right_range_beg =
        cur_range_beg + x->m_left->m_exp_len;
      std::uint64_t prefix_length = end - right_range_beg;
      const node_type *y = x->m_right;
      
      // Collect the roots of grammars to merge.
      std::vector<const node_type*> grammars_to_merge;
      while (prefix_length > 0) {
        if (y->m_exp_len == prefix_length) {
          grammars_to_merge.push_back(y);
          prefix_length -= y->m_exp_len;
        } else if (prefix_length > y->m_left->m_exp_len) {
          grammars_to_merge.push_back(y->m_left);
          prefix_length -= y->m_left->m_exp_len;
          y = y->m_right;
        } else y = y->m_left;
      }

      // Merge AVL grammars in grammars_to_merge into a single AVL grammar.
      while (grammars_to_merge.size() > 1) {
        const node_type * const right = grammars_to_merge.back();
        grammars_to_merge.pop_back();
        const node_type * const left = grammars_to_merge.back();
        grammars_to_merge.pop_back();
        grammars_to_merge.push_back(
            add_concat_nonterminal<char_type>(nonterminals, left, right));
      }
      right_grammar = grammars_to_merge.back();
    }

    // Merge left_grammar with right_grammar and return.
    // Both are guaranteed to be non-NULL.
    const node_type * const final_grammar =
      add_concat_nonterminal<char_type>(nonterminals,
          left_grammar, right_grammar);

    // Return the result.
    return final_grammar;
  }
}

//=============================================================================
// Given the LZ77 parsing (consisting of z phrases) of text T, compute
// the AVL grammar of size O(z log n) expanding to T.  The parsing is
// given as a sequence of pairs (pos, len), where either len > 0 and
// pos encodes the position of the previous occurrence in the string,
// or len = 0 and then pos contain the text symbol.
//
// TODO: the grammar at this point is guaranteed to be of size O(z log
//       n), but there might be some unused nonterminals. They should
//       be removed. Anyway, at this point, I just want to test the
//       correctness of the conversion.
//=============================================================================
template<
  typename char_type,
  typename text_offset_type>
avl_grammar<char_type> *convert_lz77_to_avl_grammar(
    const std::vector<
      std::pair<text_offset_type, text_offset_type> > &parsing) {

  // Declare types.
  typedef avl_grammar_node<char_type> node_type;
  typedef avl_grammar<char_type> grammar_type;

  // Compute the AVL grammar expanding to T.
  grammar_type *grammar = new grammar_type();
  std::uint64_t prefix_length = 0;
  for (std::uint64_t phrase_id = 0; phrase_id < parsing.size();
      ++phrase_id) {
    std::pair<text_offset_type, text_offset_type> p = parsing[phrase_id];
    std::uint64_t pos = p.first;
    std::uint64_t len = p.second;
    
    // Compute the AVL grammar expanding to phrase p.
    const node_type *phrase_root = NULL;
    if (len == 0) {

      // If this is a literal phrase, create a trivial grammar.
      phrase_root = new node_type((char_type)pos);
      grammar->m_nonterminals.push_back(phrase_root);
    } else {

      // We proceed differently, depending on whether
      // the phrase is self-overlapping. This part is
      // unadressed in the original Rytter's paper. The
      // solution is described in the proof of Theorem 6.1
      // in https://arxiv.org/abs/1910.10631v3.
      if (pos + len > prefix_length) {

        // If the phase is self-overlapping, we create the
        // nonterminal expanding to text[pos..prefix_length).
        const node_type * const suffix_nonterminal =
          add_substring_nonterminal<char_type>(
              grammar->m_nonterminals, grammar->m_root, pos, prefix_length);

        // Square the above nonterminal until
        // it reaches length >= len.
        const node_type *suffix_pow_nonterminal = suffix_nonterminal;
        std::uint64_t curlen = prefix_length - pos;
        while (curlen < len) {
          const node_type * const square =
            new node_type(suffix_pow_nonterminal, suffix_pow_nonterminal);
          grammar->m_nonterminals.push_back(square);
          curlen <<= 1;
          suffix_pow_nonterminal = square;
        }

        // Create a nonterminal expanding to the prefix
        // of exp(suffix_pow_nonterminal) of length len.
        phrase_root = add_substring_nonterminal<char_type>(
            grammar->m_nonterminals, suffix_pow_nonterminal, 0, len);
      } else {

        // Add the nonterminal expanding to phrase p.
        std::uint64_t begin = pos;
        std::uint64_t end = begin + len;
        phrase_root = add_substring_nonterminal<char_type>(
            grammar->m_nonterminals, grammar->m_root, begin, end);
      }
    }

    // Update prefix_grammar to encode the longer prefix.
    if (grammar->m_root == NULL)
      grammar->m_root = phrase_root;
    else
      grammar->m_root =
        add_concat_nonterminal<char_type>(grammar->m_nonterminals,
            grammar->m_root, phrase_root);

    // Update prefix_length.
    prefix_length += std::max(len, (std::uint64_t)1);
  }

  // Return the result.
  return grammar;
}

#endif  // __AVL_GRAMMAR_HPP_INCLUDED
