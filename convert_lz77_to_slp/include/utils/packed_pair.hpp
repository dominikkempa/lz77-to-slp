/**
 * @file    packed_pair.hpp
 * @section LICENCE
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __PACKED_PAIR_HPP_INCLUDED
#define __PACKED_PAIR_HPP_INCLUDED


template<typename S, typename T>
struct packed_pair {
  typedef packed_pair<S, T> pair_type;

  packed_pair() {}
  packed_pair(S &f, T &s) {
    first = f;
    second = s;
  }

  packed_pair(S f, T s) {
    first = f;
    second = s;
  }

  inline bool operator == (const pair_type &p) const {
    return (first == p.first) && (second == p.second);
  }

  S first;
  T second;
} __attribute__((packed));

#endif  // __PACKED_PAIR_HPP_INCLUDED
