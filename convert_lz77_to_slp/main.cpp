#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <string>
#include <ctime>
#include <unistd.h>
#include <getopt.h>

#include "include/utils.hpp"
#include "include/avl_grammar.hpp"


int main(int, char **) {

  typedef std::uint8_t char_type;
  typedef avl_grammar_node<char_type> node_type;

  // Manually create an AVL grammar encoding the string
  // Fib_7 = abaababaabaab (the 7-th fibonacci string).
  // Nonterminals are names as in the example in Rytter's paper.
  std::vector<node_type*> nonterminals;
  node_type *X1 = new node_type((char_type)'b');
  nonterminals.push_back(X1);
  node_type *X2 = new node_type((char_type)'a');
  nonterminals.push_back(X2);
  node_type *X3 = new node_type(X2, X1);
  nonterminals.push_back(X3);
  node_type *X4 = new node_type(X3, X2);
  nonterminals.push_back(X4);
  node_type *X5 = new node_type(X4, X3);
  nonterminals.push_back(X5);
  node_type *X6 = new node_type(X5, X4);
  nonterminals.push_back(X6);
  node_type *X7 = new node_type(X6, X5);
  nonterminals.push_back(X7);

  // Print expansions.
  for (std::uint64_t i = 0; i < nonterminals.size(); ++i) {
    fprintf(stderr, "exp(X%lu) = ", i + 1);
    nonterminals[i]->print_expansion();
    fprintf(stderr, ", height = %lu, explen = %lu\n",
        (std::uint64_t)nonterminals[i]->m_height,
        nonterminals[i]->m_exp_len);
  }

  // Test merging.
  {
    node_type *X_7_6 = merge_avl_grammars<char_type>(nonterminals, X7, X6);
    X_7_6->print_expansion();
    fprintf(stderr, "\n");
  }

  // Clean up.
  for (std::uint64_t i = 0; i < nonterminals.size(); ++i)
    delete nonterminals[i];
}

