/**
 * @file    repair.hpp
 * @section LICENCE
 *
 * This file is part of Lazy-AVLG v0.1.0
 * See: https://github.com/dominikkempa/lz77-to-slp
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

#ifndef __REPAIR_HPP_INCLUDED
#define __REPAIR_HPP_INCLUDED

#include <cstdlib>
#include <cstdint>
#include <algorithm>

#include "uint40.hpp"
#include "uint48.hpp"
#include "utils.hpp"
#include "hash_table.hpp"


template<typename char_type>
struct pair_record {
  std::pair<char_type, char_type> m_chars;
  std::uint64_t m_freq;
  std::uint64_t m_leftmost_occ;
  std::uint64_t m_heap_ptr;

  pair_record() {
    m_chars = std::make_pair(
        (char_type)0,
        (char_type)0);
    m_freq = 0;
    m_leftmost_occ = 0;
    m_heap_ptr = 0;
  }

  // Pairs are sorted by frequency and in case of ties by reversed
  // lex order. This is so that the value of maximal pair has the
  // largest frequency, and among pairs with maximal frequency, is
  // lexicographically first.
  inline bool operator < (const pair_record &rec) const {
    return (m_freq < rec.m_freq) ||
      (m_freq == rec.m_freq && rec.m_chars < m_chars);
  }
};

template<typename text_offset_type>
void heapify(
    std::vector<pair_record<text_offset_type> > &keys,
    std::vector<text_offset_type> &maxheap,
    std::uint64_t i) {
  std::uint64_t max_elem = i;
  while (true) {
    if (((i + 1) << 1) <= maxheap.size() &&
        keys[maxheap[max_elem]] < keys[maxheap[((i + 1) << 1) - 1]])
      max_elem = ((i + 1) << 1) - 1;
    if (((i + 1) << 1) + 1 <= maxheap.size() &&
        keys[maxheap[max_elem]] < keys[maxheap[((i + 1) << 1)]])
      max_elem = ((i + 1) << 1);
    if (max_elem != i) {
      std::swap(
          keys[maxheap[i]].m_heap_ptr,
          keys[maxheap[max_elem]].m_heap_ptr);
      std::swap(maxheap[i], maxheap[max_elem]);
      i = max_elem;
    } else break;
  }
}

template<typename text_offset_type>
void heapup(
    std::vector<pair_record<text_offset_type> > &keys,
    std::vector<text_offset_type> &maxheap,
    std::uint64_t i) {
  while (i != 0 && keys[maxheap[((i + 1) >> 1) - 1]] < keys[maxheap[i]]) {
    std::swap(
        keys[maxheap[i]].m_heap_ptr,
        keys[maxheap[((i + 1) >> 1) - 1]].m_heap_ptr);
    std::swap(maxheap[i], maxheap[((i + 1) >> 1) - 1]);
    i = ((i + 1) >> 1) - 1;
  }
}

template<>
std::uint64_t get_hash(const std::pair<std::uint32_t, std::uint32_t> &x) {
  return (std::uint64_t)x.first * (std::uint64_t)29996224275833 +
    (std::uint64_t)x.second * (std::uint64_t)14638944639703;
}

template<>
std::uint64_t get_hash(const std::pair<uint40, uint40> &x) {
  return (std::uint64_t)x.first * (std::uint64_t)29996224275833 +
    (std::uint64_t)x.second * (std::uint64_t)14638944639703;
}

template<>
std::uint64_t get_hash(const std::pair<uint48, uint48> &x) {
  return (std::uint64_t)x.first * (std::uint64_t)29996224275833 +
    (std::uint64_t)x.second * (std::uint64_t)14638944639703;
}

template<>
std::uint64_t get_hash(const std::pair<std::uint64_t, std::uint64_t> &x) {
  return (std::uint64_t)x.first * (std::uint64_t)29996224275833 +
    (std::uint64_t)x.second * (std::uint64_t)14638944639703;
}

//=============================================================================
// Returns the number of rounds + length of remaining string - 1.
//=============================================================================
template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint32_t>
std::uint64_t repair(
    const char_type * const text,
    const std::uint64_t text_length,
    std::vector<std::pair<text_offset_type, text_offset_type> > *nonterminals = NULL) {

  // Allocate current text.
  text_offset_type *cur_text = utils::allocate_array<text_offset_type>(text_length);
  for (std::uint64_t i = 0; i < text_length; ++i)
    cur_text[i] = (text_offset_type)text[i];

  // Rename symbols in the string to be
  // consecutive and start from 0.
  std::uint64_t alphabet_size = 0;
  {
    std::vector<char_type> symbols;
    for (std::uint64_t i = 0; i < text_length; ++i)
      symbols.push_back(text[i]);
    std::sort(symbols.begin(), symbols.end());
    symbols.erase(
        std::unique(symbols.begin(), symbols.end()),
        symbols.end());
    alphabet_size = symbols.size();
    for (std::uint64_t i = 0; i < text_length; ++i)
      cur_text[i] = std::lower_bound(symbols.begin(),
          symbols.end(), text[i]) - symbols.begin();
  }

  // Store original alphabet size.
  std::uint64_t original_alphabet_size = alphabet_size;

  // Allocate navigational arrays.
  text_offset_type *next = utils::allocate_array<text_offset_type>(text_length);
  text_offset_type *prev = utils::allocate_array<text_offset_type>(text_length);
  typedef pair_record<text_offset_type> record_type;
  typedef std::pair<text_offset_type, text_offset_type> key_type;
  typedef hash_table<key_type, text_offset_type> hash_table_type;
  std::vector<record_type> pair_records;
  hash_table_type pair_mapping;
  std::vector<text_offset_type> maxheap;

  // Initialize navigational structures.
  {
    std::uint64_t cur_run_length = 0;
    for (std::uint64_t i = 0; i + 1 < text_length; ++i) {
      key_type key = std::make_pair(
          (text_offset_type)cur_text[i],
          (text_offset_type)cur_text[i + 1]);
      text_offset_type *value = pair_mapping.find(key);
      ++cur_run_length;

      // Invariant: cur_run_length holds the length of
      // maximal segment ending at position i.
      if (value == NULL) {
        record_type rec;
        rec.m_chars = key;
        rec.m_freq = 1;
        rec.m_leftmost_occ = i;
        rec.m_heap_ptr = maxheap.size();
        pair_records.push_back(rec);
        maxheap.push_back(pair_records.size() - 1);
        heapup<text_offset_type>(pair_records, maxheap, rec.m_heap_ptr);
        pair_mapping.insert(key, pair_records.size() - 1);
        prev[i] = i;
        next[i] = i;
      } else {
        if ((cur_text[i] != cur_text[i + 1]) ||
            (cur_run_length & 1)) {
          pair_records[*value].m_freq += 1;
          heapup<text_offset_type>(pair_records, maxheap,
              pair_records[*value].m_heap_ptr);
        }
        std::uint64_t leftmost = pair_records[*value].m_leftmost_occ;
        std::uint64_t rightmost = prev[leftmost];
        next[rightmost] = i;
        prev[i] = rightmost;
        next[i] = leftmost;
        prev[leftmost] = i;
      }
      if (cur_text[i] != cur_text[i + 1])
        cur_run_length = 0;
    }
    prev[text_length - 1] =
      std::numeric_limits<text_offset_type>::max();
    next[text_length - 1] =
      std::numeric_limits<text_offset_type>::max();
  }

  // Main loop of the algorithm.
  std::uint64_t cur_text_length = text_length;
  std::uint64_t n_replacements = 0;
  while (cur_text_length > 1) {

    // Find the most frequent pair.
    std::uint64_t max_freq_pair_first = pair_records[maxheap[0]].m_chars.first;
    std::uint64_t max_freq_pair_second = pair_records[maxheap[0]].m_chars.second;
    std::uint64_t max_freq_pair_ptr =  pair_records[maxheap[0]].m_leftmost_occ;
    std::uint64_t max_freq = pair_records[maxheap[0]].m_freq;

    if (max_freq > 1) {
      ++n_replacements;
      if (nonterminals != NULL)
        nonterminals->push_back(std::make_pair(
              max_freq_pair_first,
              max_freq_pair_second));
      std::uint64_t left_char_pos = max_freq_pair_ptr;
      std::uint64_t cur_run_length = 0;

      // Iterate through all pairs
      // (max_freq_pair_first, max_freq_pair_second).
      do {
        std::uint64_t right_char_pos =
          (cur_text[left_char_pos + 1] ==
           std::numeric_limits<text_offset_type>::max()) ?
          (next[left_char_pos + 1] + 1) : (left_char_pos + 1);
        std::uint64_t next_pair_occ = next[left_char_pos];
        if (next_pair_occ == right_char_pos)
          next_pair_occ = next[right_char_pos];

        bool left_context_exists = 
          ((left_char_pos != 0) &&
           (cur_text[left_char_pos - 1] !=
            std::numeric_limits<text_offset_type>::max() ||
            prev[left_char_pos - 1] != (text_offset_type)0));

        // Decrease frequency of xa.
        if (left_context_exists) {
          std::uint64_t left_context_pos = 
            ((cur_text[left_char_pos - 1] ==
              std::numeric_limits<text_offset_type>::max()) ?
             (prev[left_char_pos - 1] - 1) : (left_char_pos - 1));
          if ((std::uint64_t)cur_text[left_context_pos] != max_freq_pair_first) {
            key_type key_xa = std::make_pair(
                cur_text[left_context_pos],
                cur_text[left_char_pos]);
            text_offset_type *value = pair_mapping.find(key_xa);
            pair_records[*value].m_freq -= 1;
            heapify<text_offset_type>(pair_records, maxheap,
                pair_records[*value].m_heap_ptr);
            if (pair_records[*value].m_freq > 0) {
              if (pair_records[*value].m_leftmost_occ == left_context_pos)
                pair_records[*value].m_leftmost_occ = next[left_context_pos];
              next[prev[left_context_pos]] = next[left_context_pos];
              prev[next[left_context_pos]] = prev[left_context_pos];
            }
          } else {

            // Invariant: max_freq_pair_first != max_freq_pair_second.
            // Case xa = aa. In order to know whether we need to
            // decrease the frequency of xa, we need to know
            // the (parity of the) length of the run containing
            // this pair.
            // It is easy to see that this scan (and the analogous
            // for the other side) will not encounter more pairs
            // than the total number of elements in this do-while loop,
            // so we can affort to do this naively.
            std::uint64_t run_length = 2;
            {
              std::uint64_t ptr = left_context_pos;
              bool exists_left_neighbor = ((ptr > 0) && 
                  (cur_text[ptr - 1] !=
                   std::numeric_limits<text_offset_type>::max() ||
                   (std::uint64_t)prev[ptr - 1] != 0));
              while (exists_left_neighbor) {
                std::uint64_t left_neighbor_pos =
                  (cur_text[ptr - 1] !=
                   std::numeric_limits<text_offset_type>::max()) ?
                  (ptr - 1) : (prev[ptr - 1] - 1);
                if ((std::uint64_t)cur_text[left_neighbor_pos] ==
                    max_freq_pair_first) {
                  ++run_length;
                  ptr = left_neighbor_pos;
                  exists_left_neighbor = ((ptr > 0) && 
                      (cur_text[ptr - 1] !=
                       std::numeric_limits<text_offset_type>::max() ||
                       (std::uint64_t)prev[ptr - 1] != 0));
                } else break;
              }
            }

            if (!(run_length & 1)) {
              key_type key_xa = std::make_pair(
                  cur_text[left_context_pos],
                  cur_text[left_char_pos]);
              text_offset_type *value = pair_mapping.find(key_xa);
              pair_records[*value].m_freq -= 1;
              heapify<text_offset_type>(pair_records, maxheap,
                  pair_records[*value].m_heap_ptr);
              if (pair_records[*value].m_freq > 0) {
                if (pair_records[*value].m_leftmost_occ == left_context_pos)
                  pair_records[*value].m_leftmost_occ = next[left_context_pos];
                next[prev[left_context_pos]] = next[left_context_pos];
                prev[next[left_context_pos]] = prev[left_context_pos];
              }
            } else {
              next[prev[left_context_pos]] = next[left_context_pos];
              prev[next[left_context_pos]] = prev[left_context_pos];
            }
          }
        }

        bool right_context_exists = 
          ((right_char_pos + 1 != text_length) &&
           (cur_text[right_char_pos + 1] !=
            std::numeric_limits<text_offset_type>::max() ||
            next[right_char_pos + 1] + 1 != (text_offset_type)text_length));

        // Decrease frequency of by.
        if (right_context_exists) {
          std::uint64_t right_context_pos = 
            ((cur_text[right_char_pos + 1] ==
              std::numeric_limits<text_offset_type>::max()) ?
             (next[right_char_pos + 1] + 1) : (right_char_pos + 1));

          // Note: if a == b == y, we do not need to do anything. Erasing
          // a pair by will be accompliehsed later since we anyway remove
          // all occurrences of pair ab = by.
          if (max_freq_pair_second != cur_text[right_context_pos]) {
            key_type key_by = std::make_pair(
                cur_text[right_char_pos],
                cur_text[right_context_pos]);
            text_offset_type *value = pair_mapping.find(key_by);
            pair_records[*value].m_freq -= 1;
            heapify<text_offset_type>(pair_records, maxheap,
                pair_records[*value].m_heap_ptr);
            if (pair_records[*value].m_freq > 0) {
              if (pair_records[*value].m_leftmost_occ == right_char_pos)
                pair_records[*value].m_leftmost_occ = next[right_char_pos];
              next[prev[right_char_pos]] = next[right_char_pos];
              prev[next[right_char_pos]] = prev[right_char_pos];
            }
          } else {
            if (max_freq_pair_first != max_freq_pair_second) {

              // Case by = bb. In order to know whether we need to
              // decrease the frequency of by, we need to know
              // the (parity of the) length of the run containing
              // this pair.
              std::uint64_t run_length = 2;
              {
                std::uint64_t ptr = right_context_pos;
                bool exists_right_neighbor = ((ptr + 1 != text_length) && 
                    (cur_text[ptr + 1] !=
                     std::numeric_limits<text_offset_type>::max() ||
                     next[ptr + 1] + 1 != text_length));
                while (exists_right_neighbor) {
                  std::uint64_t right_neighbor_pos =
                    (cur_text[ptr + 1] !=
                     std::numeric_limits<text_offset_type>::max()) ?
                    (ptr + 1) : (next[ptr + 1] + 1);
                  if ((std::uint64_t)cur_text[right_neighbor_pos] ==
                      max_freq_pair_second) {
                    ++run_length;
                    ptr = right_neighbor_pos;
                    exists_right_neighbor = ((ptr + 1 != text_length) && 
                        (cur_text[ptr + 1] !=
                         std::numeric_limits<text_offset_type>::max() ||
                         next[ptr + 1] + 1 != text_length));
                  } else break;
                }
              }

              if (!(run_length & 1)) {
                key_type key_by = std::make_pair(
                    cur_text[right_char_pos],
                    cur_text[right_context_pos]);
                text_offset_type *value = pair_mapping.find(key_by);
                pair_records[*value].m_freq -= 1;
                heapify<text_offset_type>(pair_records, maxheap,
                    pair_records[*value].m_heap_ptr);
                if (pair_records[*value].m_freq > 0) {
                  if (pair_records[*value].m_leftmost_occ == right_char_pos)
                    pair_records[*value].m_leftmost_occ = next[right_char_pos];
                  next[prev[right_char_pos]] = next[right_char_pos];
                  prev[next[right_char_pos]] = prev[right_char_pos];
                }
              } else {
                key_type key_by = std::make_pair(
                    cur_text[right_char_pos],
                    cur_text[right_context_pos]);
                text_offset_type *value = pair_mapping.find(key_by);
                if (pair_records[*value].m_leftmost_occ == right_char_pos)
                  pair_records[*value].m_leftmost_occ = next[right_char_pos];
                next[prev[right_char_pos]] = next[right_char_pos];
                prev[next[right_char_pos]] = prev[right_char_pos];
              }
            }
          }
        }

        // Update cur_run_length.
        if (left_context_exists) {
          std::uint64_t left_context_pos = 
            ((cur_text[left_char_pos - 1] ==
              std::numeric_limits<text_offset_type>::max()) ?
             (prev[left_char_pos - 1] - 1) : (left_char_pos - 1));
          if ((std::uint64_t)cur_text[left_context_pos] == alphabet_size)
            ++cur_run_length;
          else cur_run_length = 1;
        } else cur_run_length = 1;

        // Increase the frequency of xc.
        if (left_context_exists) {
          std::uint64_t left_context_pos = 
            ((cur_text[left_char_pos - 1] ==
              std::numeric_limits<text_offset_type>::max()) ?
             (prev[left_char_pos - 1] - 1) : (left_char_pos - 1));
          key_type key_xc =
            std::make_pair(cur_text[left_context_pos], alphabet_size);
          text_offset_type *value = pair_mapping.find(key_xc);
          if (value == NULL || pair_records[*value].m_freq == 0) {
            if (value == NULL) {
              record_type rec;
              rec.m_chars = key_xc;
              rec.m_freq = 1;
              rec.m_leftmost_occ = left_context_pos;
              rec.m_heap_ptr = maxheap.size();
              pair_records.push_back(rec);
              maxheap.push_back(pair_records.size() - 1);
              heapup<text_offset_type>(pair_records, maxheap, rec.m_heap_ptr);
              pair_mapping.insert(key_xc, pair_records.size() - 1);
            } else {
              pair_records[*value].m_freq = 1;
              pair_records[*value].m_leftmost_occ = left_context_pos;
              heapup<text_offset_type>(pair_records, maxheap,
                  pair_records[*value].m_heap_ptr);
            }
            prev[left_context_pos] = left_context_pos;
            next[left_context_pos] = left_context_pos;
          } else {
            if ((std::uint64_t)cur_text[left_context_pos] != alphabet_size ||
                (!(cur_run_length & 1))) {
              pair_records[*value].m_freq += 1;
              heapup<text_offset_type>(pair_records, maxheap,
                  pair_records[*value].m_heap_ptr);
            }
            prev[left_context_pos] =
              prev[pair_records[*value].m_leftmost_occ];
            next[left_context_pos] =
              pair_records[*value].m_leftmost_occ;
            next[prev[left_context_pos]] = left_context_pos;
            prev[next[left_context_pos]] = left_context_pos;
          }
        }

        // Increase the frequency of cy.
        if (right_context_exists) {
          std::uint64_t right_context_pos = 
            ((cur_text[right_char_pos + 1] ==
              std::numeric_limits<text_offset_type>::max()) ?
             (next[right_char_pos + 1] + 1) : (right_char_pos + 1));
          key_type key_cy =
            std::make_pair(alphabet_size, cur_text[right_context_pos]);
          text_offset_type *value = pair_mapping.find(key_cy);
          if (value == NULL || pair_records[*value].m_freq == 0) {
            if (value == NULL) {
              record_type rec;
              rec.m_chars = key_cy;
              rec.m_freq = 1;
              rec.m_leftmost_occ = left_char_pos;
              rec.m_heap_ptr = maxheap.size();
              pair_records.push_back(rec);
              maxheap.push_back(pair_records.size() - 1);
              heapup<text_offset_type>(pair_records, maxheap, rec.m_heap_ptr);
              pair_mapping.insert(key_cy, pair_records.size() - 1);
            } else {
              pair_records[*value].m_freq = 1;
              pair_records[*value].m_leftmost_occ = left_char_pos;
              heapup<text_offset_type>(pair_records, maxheap,
                  pair_records[*value].m_heap_ptr);
            }
            prev[left_char_pos] = left_char_pos;
            next[left_char_pos] = left_char_pos;
          } else {
            pair_records[*value].m_freq += 1;
            heapup<text_offset_type>(pair_records, maxheap,
                pair_records[*value].m_heap_ptr);
            prev[left_char_pos] =
              prev[pair_records[*value].m_leftmost_occ];
            next[left_char_pos] =
              pair_records[*value].m_leftmost_occ;
            next[prev[left_char_pos]] = left_char_pos;
            prev[next[left_char_pos]] = left_char_pos;
          }
        }

        // Set c at position cur_text[left_char_pos]
        // and erase cur_text[right_char_pos].
        cur_text[left_char_pos] = alphabet_size;
        cur_text[right_char_pos] =
          std::numeric_limits<text_offset_type>::max();

        // Update the threading of gaps.
        if (right_char_pos + 1 == text_length ||
            cur_text[right_char_pos + 1] !=
            std::numeric_limits<text_offset_type>::max()) {
          next[left_char_pos + 1] = right_char_pos;
          prev[next[left_char_pos + 1]] = left_char_pos + 1;
        } else {
          next[left_char_pos + 1] = next[right_char_pos + 1];
          prev[next[left_char_pos + 1]] = left_char_pos + 1;
        }

        --cur_text_length;
        left_char_pos = next_pair_occ;
      } while (left_char_pos != max_freq_pair_ptr);

      ++alphabet_size;
      {

        // Set to zero the frequency of pair
        // (max_freq_pair_first, max_freq_pair_second).
        key_type key_ab = std::make_pair(
            max_freq_pair_first,
            max_freq_pair_second);
        text_offset_type *value = pair_mapping.find(key_ab);
        pair_records[*value].m_freq = 0;
        heapify<text_offset_type>(pair_records, maxheap,
            pair_records[*value].m_heap_ptr);
      }
    } else break;
  }

  // Clean up.
  utils::deallocate(prev);
  utils::deallocate(next);
  utils::deallocate(cur_text);

  // Print stats.
  fprintf(stderr, "Original alphabet size = %lu\n", original_alphabet_size);
  fprintf(stderr, "Number of binary nonterminals = %lu\n", n_replacements);
  fprintf(stderr, "Root expansion length = %lu\n", cur_text_length);

  // Return result.
  std::uint64_t grammar_size =
    original_alphabet_size +
    2 * n_replacements +
    (cur_text_length - 1);
  return grammar_size;
}

#endif  // __REPAIR_HPP_INCLUDED

