/**
 * @file    compute_sa.hpp
 * @section LICENCE
 *
 * Copyright (C) 2017
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
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

#ifndef __COMPUTE_SA_HPP_INCLUDED
#define __COMPUTE_SA_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>

#include "sais.hxx"


template<typename char_type,
  typename text_offset_type>
void compute_sa(
    const char_type *,
    std::uint64_t,
    text_offset_type *) {
  fprintf(stderr, "\nError: compute_sa not supported for these types!\n");
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr ,"sizeof(text_offset_type) = %lu\n",
      sizeof(text_offset_type));
  std::exit(EXIT_FAILURE);
}

template<>
void compute_sa(
    const std::uint8_t *text,
    std::uint64_t text_length,
    std::uint32_t *sa) {
  saisxx((const unsigned char *)text, (int *)sa, (int)text_length);
}

template<>
void compute_sa(
    const std::uint8_t *text,
    std::uint64_t text_length,
    std::uint64_t *sa) {
  saisxx((const unsigned char *)text, (long *)sa, (long)text_length);
}

#endif  // __COMPUTE_SA_HPP_INCLUDED
