/**
 * @file    hash_table.hpp
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

#ifndef __HASH_TABLE_HPP_INCLUDED
#define __HASH_TABLE_HPP_INCLUDED

#include <limits>
#include <algorithm>

#include "utils.hpp"

template<typename key_type>
std::uint64_t get_hash(const key_type &) {
  fprintf(stderr, "\n\nError: calling a generic get_hash.\n");
  std::exit(EXIT_FAILURE);
  return (std::uint64_t)0;
}

template<typename KeyType,
  typename ValueType,
  typename SizeType = std::uint64_t>
class hash_table {
  public:
    typedef KeyType key_type;
    typedef ValueType value_type;
    typedef SizeType size_type;  // able to encode [0..2*max_items)

  private:
    template<typename T, typename S, typename U>
    struct hash_item {
      T m_key;
      S m_value;
      U m_next;
    };

    typedef hash_item<key_type, value_type, size_type> item_type;

    // After the first item has been inserted into hash table:
    // the following invariant holds at all times:
    // m_item_count <= m_bucket_count < 2 * m_item_count.
    std::uint64_t m_bucket_count;
    std::uint64_t m_item_count;

    item_type *m_items;
    size_type *m_buckets;

  public:
    hash_table() {
      m_bucket_count = 1;
      m_item_count = 0;

      m_items = utils::allocate_array<item_type>(m_bucket_count);
      m_buckets = utils::allocate_array<size_type>(m_bucket_count);
      std::fill(m_buckets, m_buckets + m_bucket_count,
          std::numeric_limits<size_type>::max());
    }

  private:
    void enlarge() {

      // Allocate new arrays.
      std::uint64_t new_bucket_count = m_bucket_count * 2;
      item_type *new_items =
        utils::allocate_array<item_type>(new_bucket_count);
      size_type *new_buckets =
        utils::allocate_array<size_type>(new_bucket_count);
      std::fill(new_buckets, new_buckets + new_bucket_count,
          std::numeric_limits<size_type>::max());

      // Rehash all items.
      std::uint64_t item_count = 0;
      for (std::uint64_t i = 0; i < m_bucket_count; ++i) {
        std::uint64_t j = m_buckets[i];
        std::uint64_t size_type_max =
          (std::uint64_t)std::numeric_limits<size_type>::max();
        while (j != size_type_max) {
          std::uint64_t hash =
            get_hash(m_items[j].m_key) & (new_bucket_count - 1);
          std::uint64_t new_item = item_count++;
          new_items[new_item].m_next = new_buckets[hash];
          new_items[new_item].m_key = m_items[j].m_key;
          new_items[new_item].m_value = m_items[j].m_value;
          new_buckets[hash] = new_item;
          j = m_items[j].m_next;
        }
      }

      // Update arrays.
      m_bucket_count = new_bucket_count;
      utils::deallocate(m_buckets);
      utils::deallocate(m_items);
      m_items = new_items;
      m_buckets = new_buckets;
    }

  public:
    void insert(const key_type &key, const value_type &value) {
      if (m_item_count == m_bucket_count)
        enlarge();

      std::uint64_t hash = get_hash(key) & (m_bucket_count - 1);
      std::uint64_t new_item = m_item_count++;
      m_items[new_item].m_next = m_buckets[hash];
      m_items[new_item].m_key = key;
      m_items[new_item].m_value = value;
      m_buckets[hash] = new_item;
    }

    value_type* find(const key_type &key) {
      std::uint64_t hash = get_hash(key) & (m_bucket_count - 1);
      std::uint64_t j = m_buckets[hash];
      while (j != std::numeric_limits<size_type>::max()) {
        if (m_items[j].m_key == key) return &(m_items[j].m_value);
        j = m_items[j].m_next;
      }
      return NULL;
    }


    const value_type* find(const key_type &key) const {
      std::uint64_t hash = get_hash(key) & (m_bucket_count - 1);
      std::uint64_t j = m_buckets[hash];
      while (j != std::numeric_limits<size_type>::max()) {
        if (m_items[j].m_key == key) return &(m_items[j].m_value);
        j = m_items[j].m_next;
      }
      return NULL;
    }

    ~hash_table() {
      utils::deallocate(m_buckets);
      utils::deallocate(m_items);
    }
};

#endif  // __HASH_TABLE_HPP_INCLUDED
