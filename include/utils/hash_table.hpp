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

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <algorithm>

#include "utils.hpp"
#include "space_efficient_vector.hpp"


//=============================================================================
// Generic hash function. Needs to be specialized for the key type.
//=============================================================================
template<typename key_type>
std::uint64_t get_hash(const key_type &) {
  fprintf(stderr, "\n\nError: calling a generic get_hash.\n");
  std::exit(EXIT_FAILURE);
  return (std::uint64_t)0;
}

//=============================================================================
// Simple chaining hash table.
//=============================================================================
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
    } __attribute__((packed));

    typedef hash_item<key_type, value_type, size_type> item_type;

    //=========================================================================
    // After the first item has been inserted into hash table:
    // the following invariant holds at all times:
    // m_item.size() <= m_bucket_count < 2 * m_item.size().
    //=========================================================================
    std::uint64_t m_bucket_count;

    space_efficient_vector<item_type> m_items;
    size_type *m_buckets;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    hash_table() {
      m_bucket_count = 1;
      m_buckets = utils::allocate_array<size_type>(m_bucket_count);
      std::fill(m_buckets, m_buckets + m_bucket_count,
          std::numeric_limits<size_type>::max());
    }

  private:

    //=========================================================================
    // Double the capacity of the hash table.
    //=========================================================================
    void enlarge() {

      // Allocate new arrays.
      m_bucket_count <<= 1;
      utils::deallocate(m_buckets);
      m_buckets = utils::allocate_array<size_type>(m_bucket_count);
      std::fill(m_buckets, m_buckets + m_bucket_count,
          std::numeric_limits<size_type>::max());

      // Rehash all items.
      for (std::uint64_t i = 0; i < m_items.size(); ++i) {
        const key_type &key = m_items[i].m_key;
        const std::uint64_t hash = get_hash(key) & (m_bucket_count - 1);
        m_items[i].m_next = m_buckets[hash];
        m_buckets[hash] = i;
      }
    }

  public:

    //=========================================================================
    // Insert a new key (and the associated value).
    //=========================================================================
    void insert(const key_type &key, const value_type &value) {
      if (m_items.size() == m_bucket_count)
        enlarge();

      const std::uint64_t hash = get_hash(key) & (m_bucket_count - 1);
      const std::uint64_t new_item_pos = m_items.size();
      item_type new_item;
      new_item.m_next = m_buckets[hash];
      new_item.m_key = key;
      new_item.m_value = value;
      m_items.push_back(new_item);
      m_buckets[hash] = new_item_pos;
    }

    //=========================================================================
    // Find the value associated with a given key.
    //=========================================================================
    value_type* find(const key_type &key) {
      const std::uint64_t hash = get_hash(key) & (m_bucket_count - 1);
      std::uint64_t j = m_buckets[hash];
      while (j != std::numeric_limits<size_type>::max()) {
        if (m_items[j].m_key == key) return &(m_items[j].m_value);
        j = m_items[j].m_next;
      }
      return NULL;
    }

    //=========================================================================
    // Find a value associated with a given key.
    //=========================================================================
    const value_type* find(const key_type &key) const {
      const std::uint64_t hash = get_hash(key) & (m_bucket_count - 1);
      std::uint64_t j = m_buckets[hash];
      while (j != std::numeric_limits<size_type>::max()) {
        if (m_items[j].m_key == key) return &(m_items[j].m_value);
        j = m_items[j].m_next;
      }
      return NULL;
    }

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~hash_table() {
      utils::deallocate(m_buckets);
    }

    //=========================================================================
    // Return the number of items in hash table.
    //=========================================================================
    std::uint64_t size() const {
      return m_items.size();
    }

    //=========================================================================
    // Erase all elements, but keep current capacity.
    //=========================================================================
    void reset() {
      m_items.set_empty();
      std::fill(m_buckets, m_buckets + m_bucket_count,
          std::numeric_limits<size_type>::max());
    }

    //=========================================================================
    // Return RAM used by hash table.
    //=========================================================================
    std::uint64_t ram_use() const {
      const std::uint64_t m_items_ram_use =
        m_items.ram_use();
      const std::uint64_t m_buckets_ram_use =
        sizeof(size_type) * m_bucket_count;
      const std::uint64_t total =
        m_items_ram_use +
        m_buckets_ram_use;
      return total;
    }
};

#endif  // __HASH_TABLE_HPP_INCLUDED
