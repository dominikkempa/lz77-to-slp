#ifndef __AVL_GRAMMAR_ADD_CONCAT_NONTERMINAL_HPP
#define __AVL_GRAMMAR_ADD_CONCAT_NONTERMINAL_HPP

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "avl_grammar_node.hpp"


//=============================================================================
// Given two nonterminals `left' and `right' of the same grammar
// expanding respectively to strings X and Y, add to the set of
// nonterminals, the nonterminals that expands to XY, and return the
// pointer to it. The root of the original grammar remains unchanged.
// The function resembles the procedure to merge two AVL trees.
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
        add_concat_nonterminal<char_type>(
            nonterminals, left->m_right, right);
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

#endif  // __AVL_GRAMMAR_ADD_CONCAT_NONTERMINAL_HPP
