Lazy-AVLG
=========


Description
-----------

A C++ implementation of a fast and space-efficient algorithm
converting the LZ77 parsing into a straight-line grammar. The
algorithm is based on AVL grammars introduced by Rytter [TCS,
2003]. Our main improvements are: delayed ("lazy") merging of
nonterminals, and utilizing Karp-Rabin fingerprints. The details of
our algorithm are described in
[https://arxiv.org/abs/2105.11052](https://arxiv.org/abs/2105.11052).

The latest version of this code is available from
[https://github.com/dominikkempa/lz77-to-slp](https://github.com/dominikkempa/lz77-to-slp).



Compilation and usage
---------------------

Lazy-AVLG is compiled by simply typing `make` in the directory
containing this README. This will build the executable called
`lz_to_grammar`. The simplest usage of Lazy-AVLG is as follows.
Suppose the input LZ77 parsing (which can be computed using the
program in `tools/text_to_lz`) is located in `/data/input.txt.lz77`.
Then, to compute the lazy AVL grammar of `input.txt`, type:

    $ ./lz_to_grammar /data/input.txt.lz77

This will write the output grammar to `/data/input.txt.lz77.slg`.  For
details of the encoding of the grammar, see the function
`write_to_file` in `include/lazy_avl_grammar/lazy_avl_grammar.hpp`.
