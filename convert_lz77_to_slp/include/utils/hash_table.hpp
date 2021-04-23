#ifndef __HASH_TABLE_HPP_INCLUDED
#define __HASH_TABLE_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <algorithm>

#include "utils.hpp"


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
    };

    typedef hash_item<key_type, value_type, size_type> item_type;

    //=========================================================================
    // After the first item has been inserted into hash table:
    // the following invariant holds at all times:
    // m_item_count <= m_bucket_count < 2 * m_item_count.
    //=========================================================================
    std::uint64_t m_bucket_count;
    std::uint64_t m_item_count;

    item_type *m_items;
    size_type *m_buckets;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    hash_table() {
      m_bucket_count = 1;
      m_item_count = 0;

      m_items = utils::allocate_array<item_type>(m_bucket_count);
      m_buckets = utils::allocate_array<size_type>(m_bucket_count);
      std::fill(m_buckets, m_buckets + m_bucket_count,
          std::numeric_limits<size_type>::max());
    }

  private:

    //=========================================================================
    // Double the capacity of the hash table. There is room for
    // improvement in the way we rehash all elements.
    //=========================================================================
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

    //=========================================================================
    // Insert a new key (and the associated value).
    //=========================================================================
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

    //=========================================================================
    // Find the value associated with a given key.
    //=========================================================================
    value_type* find(const key_type &key) {
      std::uint64_t hash = get_hash(key) & (m_bucket_count - 1);
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
      std::uint64_t hash = get_hash(key) & (m_bucket_count - 1);
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
      utils::deallocate(m_items);
    }

    //=========================================================================
    // Return the number of items in hash table.
    //=========================================================================
    std::uint64_t size() const {
      return m_item_count;
    }

    //=========================================================================
    // Return RAM used by hash table.
    //=========================================================================
    std::uint64_t ram_use() const {
      const std::uint64_t m_items_ram_use =
        sizeof(item_type) * m_bucket_count;
      const std::uint64_t m_buckets_ram_use =
        sizeof(size_type) * m_bucket_count;
      const std::uint64_t total =
        m_items_ram_use +
        m_buckets_ram_use;
      return total;
    }
};

#endif  // __HASH_TABLE_HPP_INCLUDED
