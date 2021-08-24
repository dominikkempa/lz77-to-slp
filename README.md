Lazy-AVLG
=========


Description
-----------

A C++ implementation of a fast and space-efficient algorithm
converting the LZ77 parsing into a straight-line grammar. The
algorithm is based on AVL grammars introduced by Rytter [TCS,
2003]. Our main improvements are: delayed ("lazy") merging of
nonterminals, and utilizing Karp-Rabin fingerprints.

For more detailed description of the algorithm and reports of
experimental evaluation, refer to the following paper.

    @inproceedings{lazyavl,
      author =      {Dominik Kempa and Ben Langmead},
      title =       {Fast and Space-Efficient Construction of {AVL}
                     Grammars from the {LZ77} Parsing},
      booktitle =   {29th Annual European Symposium on Algorithms
                     (ESA 2021)},
      pages =       {56:1--56:14},
      year =        {2021},
      doi =         {10.4230/LIPIcs.ESA.2021.56},
    }

The above paper is available on arXiv at
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



Terms of use
------------

If you use this code, please cite the paper mentioned above.
Lazy-AVLG is released under the MIT/X11 license. See
the file LICENCE for more details.

