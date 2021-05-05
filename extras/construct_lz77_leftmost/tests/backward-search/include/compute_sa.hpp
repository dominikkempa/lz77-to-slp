/**
 * @file    compute_sa.hpp
 * @section LICENCE
 *
 * Copyright (C) 2021
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

#include "uint40.hpp"
#include "uint48.hpp"
#include "sais.hxx"
#include "naive_compute_sa.hpp"


//=============================================================================
// Compute SA of a given text[0..text_length)
// and write to sa[0..text_length).
//=============================================================================
template<
  typename char_type,
  typename text_offset_type>
void compute_sa(
    const char_type * const text,
    const std::uint64_t text_length,
    text_offset_type * const sa) {

  naive_compute_sa(text, text_length, sa);
}

//=============================================================================
// Instantiation of compute_sa for
// text_offset_type == std::uint32_t
// and char_type == std::uint8_t.
//=============================================================================
template<>
void compute_sa(
    const std::uint8_t * const text,
    const std::uint64_t text_length,
    std::uint32_t * const sa) {

  // Run sais.
  saisxx<const std::uint8_t *, std::int32_t*, std::int32_t>(
      text, (std::int32_t *)sa, (std::int32_t)text_length);
}

//=============================================================================
// Instantiation of compute_sa for
// text_offset_type == std::uint64_t
// and char_type == std::uint8_t.
//=============================================================================
template<>
void compute_sa(
    const std::uint8_t * const text,
    const std::uint64_t text_length,
    std::uint64_t * const sa) {

  // Run sais.
  saisxx<const std::uint8_t *, std::int64_t*, std::int64_t>(text,
      (std::int64_t *)sa, (std::int64_t)text_length);
}

//=============================================================================
// Instantiation (not space efficient) of compute_sa for
// text_offset_type == uint40 and char_type == std::uint8_t.
//=============================================================================
template<>
void compute_sa(
    const std::uint8_t * const text,
    const std::uint64_t text_length,
    uint40 * const sa) {

  std::uint64_t * const sa64 = new std::uint64_t[text_length];
  compute_sa(text, text_length, sa64);
  for (std::uint64_t i = 0; i < text_length; ++i)
    sa[i] = sa64[i];
  delete[] sa64;
}

//=============================================================================
// Instantiation (not space efficient) of compute_sa for
// text_offset_type == uint48 and char_type == std::uint8_t.
//=============================================================================
template<>
void compute_sa(
    const std::uint8_t * const text,
    const std::uint64_t text_length,
    uint48 * const sa) {

  std::uint64_t * const sa64 = new std::uint64_t[text_length];
  compute_sa(text, text_length, sa64);
  for (std::uint64_t i = 0; i < text_length; ++i)
    sa[i] = sa64[i];
  delete[] sa64;
}

#endif  // __COMPUTE_SA_HPP_INCLUDED
