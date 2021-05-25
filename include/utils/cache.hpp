/**
 * @file    cache.hpp
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

#ifndef __CACHE_HPP_INCLUDED
#define __CACHE_HPP_INCLUDED

#include <cstdint>
#include <algorithm>

#include "utils.hpp"
#include "hash_table.hpp"
#include "packed_pair.hpp"


template<
  typename key_type,
  typename value_type>
struct cache {

  //===========================================================================
  // Declare types.
  //===========================================================================
  typedef cache<key_type, value_type> cache_type;
  typedef packed_pair<key_type, value_type> pair_type;

  private:

    //=========================================================================
    // Class members.
    //=========================================================================
    pair_type *m_tab;
    hash_table<key_type, value_type> m_hash_table;
    const std::uint64_t m_max_size;
    std::uint64_t m_size;
    std::uint64_t m_query_counter;
    std::uint64_t m_cache_misses;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    cache(const std::uint64_t max_size = (1 << 17))
        : m_max_size(max_size) {
      m_size = 0;
      m_query_counter = 0;
      m_cache_misses = 0;
      m_tab = utils::allocate_array<pair_type>(m_max_size);
    }

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~cache() {
      utils::deallocate(m_tab);
    }

    //=========================================================================
    // Lookup an element in cache.
    //=========================================================================
    const value_type *lookup(const key_type &key) {
      const value_type * const ret = m_hash_table.find(key);
      ++m_query_counter;
      if (ret == NULL)
        ++m_cache_misses;
      return ret;
    }

    //=========================================================================
    // Inset a new element into cache.
    //=========================================================================
    void insert(const key_type &key, const value_type &value) {
      if (m_size == m_max_size)
        reduce();
      m_tab[m_size++] = pair_type(key, value);
      m_hash_table.insert(key, value);
    }

    //=========================================================================
    // Return RAM use.
    //=========================================================================
    std::uint64_t ram_use() const {
      const std::uint64_t tab_ram_use =
        sizeof(pair_type) * m_max_size;
      const std::uint64_t hash_table_ram_use =
        m_hash_table.ram_use();
      const std::uint64_t total =
        tab_ram_use +
        hash_table_ram_use;
      return total;
    }

    //=========================================================================
    // Return the number of queries.
    //=========================================================================
    inline std::uint64_t get_query_counter() const {
      return m_query_counter;
    }

    //=========================================================================
    // Return the number of cache misses.
    //=========================================================================
    inline std::uint64_t get_cache_misses() const {
      return m_cache_misses;
    }

  private:

    //=========================================================================
    // Discard half of the oldest elements in cache.
    //=========================================================================
    void reduce() {
      std::uint64_t new_size = (m_size + 1) / 2;
      std::uint64_t to_discard = m_size - new_size;
      for (std::uint64_t i = 0; i < new_size; ++i)
        m_tab[i] = m_tab[to_discard + i];
      m_size = new_size;
      m_hash_table.reset();
      for (std::uint64_t i = 0; i < m_size; ++i)
        m_hash_table.insert(m_tab[i].first, m_tab[i].second);
    }
};

#endif  // __CACHE_HPP_INCLUDED
