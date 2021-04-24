/**
 * @file    space_efficient_vector.hpp
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

#ifndef __SPACE_EFFICIENT_VECTOR_HPP_INCLUDED
#define __SPACE_EFFICIENT_VECTOR_HPP_INCLUDED

#include <cstdint>
#include <algorithm>

#include "utils.hpp"


//=============================================================================
// A class that implements essentially the same functionality as std::vector,
// but with much smaller peak RAM usage during reallocation (which in most
// implementations is  not inplace). Slowdown of access it very moderate.
//=============================================================================
template<typename ValueType>
class space_efficient_vector {
  public:

    //=========================================================================
    // Declare types.
    //=========================================================================
    typedef ValueType value_type;

  private:

    //=========================================================================
    // Set the number of blocks. This controls the peak RAM use.
    //=========================================================================
    static const std::uint64_t max_blocks = 32;

    //=========================================================================
    // Class members.
    //=========================================================================
    value_type **m_blocks;
    std::uint64_t m_block_size_log;
    std::uint64_t m_block_size_mask;
    std::uint64_t m_block_size;
    std::uint64_t m_allocated_blocks;
    std::uint64_t m_cur_block_id;
    std::uint64_t m_cur_block_filled;
    std::uint64_t m_size;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    space_efficient_vector() {
      m_size = 0;
      m_block_size_log = 0;
      m_block_size_mask = 0;
      m_block_size = 1;
      m_allocated_blocks = 1;
      m_cur_block_filled = 0;
      m_cur_block_id = 0;
      m_blocks = utils::allocate_array<value_type*>(max_blocks);
      m_blocks[0] = utils::allocate_array<value_type>(m_block_size);
    }

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~space_efficient_vector() {
      for (std::uint64_t block_id = 0;
          block_id < m_allocated_blocks; ++block_id)
        utils::deallocate(m_blocks[block_id]);
      utils::deallocate(m_blocks);
    }

    //=========================================================================
    // Return vector size.
    //=========================================================================
    inline std::uint64_t size() const {
      return m_size;
    }

    //=========================================================================
    // Check if the vector is empty.
    //=========================================================================
    inline bool empty() const {
      return m_size == 0;
    }

    //=========================================================================
    // Set the vector as empty (without deallocating space).
    //=========================================================================
    void set_empty() {
      m_size = 0;
      m_cur_block_filled = 0;
      m_cur_block_id = 0;
    }

    //=========================================================================
    // Remove the last element (without reducing capacity).
    //=========================================================================
    inline void pop_back() {
      --m_size;
      --m_cur_block_filled;
      if (m_cur_block_filled == 0 && m_cur_block_id > 0) {
        --m_cur_block_id;
        m_cur_block_filled = m_block_size;
      }
    }

    //=========================================================================
    // Push a new item at the end.
    //=========================================================================
    inline void push_back(const value_type &value) {
      if (m_cur_block_filled == m_block_size &&
          m_cur_block_id + 1 == max_blocks) {
        const std::uint64_t new_block_size = m_block_size * 2;
        for (std::uint64_t block_id = 0;
            block_id < m_allocated_blocks; block_id += 2) {
          value_type *newblock =
            utils::allocate_array<value_type>(new_block_size);
          std::copy(m_blocks[block_id],
              m_blocks[block_id] + m_block_size, newblock);
          std::copy(m_blocks[block_id + 1],
              m_blocks[block_id + 1] + m_block_size,
              newblock + m_block_size);
          utils::deallocate(m_blocks[block_id]);
          utils::deallocate(m_blocks[block_id + 1]);
          m_blocks[block_id / 2] = newblock;
        }
        m_allocated_blocks = max_blocks / 2;
        m_block_size = new_block_size;
        m_block_size_mask = new_block_size - 1;
        ++m_block_size_log;
        m_cur_block_id = m_allocated_blocks - 1;
        m_cur_block_filled = new_block_size;
      }

      if (m_cur_block_filled == m_block_size) {
        ++m_cur_block_id;
        m_cur_block_filled = 0;
        if (m_cur_block_id == m_allocated_blocks) {
          ++m_allocated_blocks;
          m_blocks[m_cur_block_id] =
            utils::allocate_array<value_type>(m_block_size);
        }
      }

      m_blocks[m_cur_block_id][m_cur_block_filled++] = value;
      ++m_size;
    }

    //=========================================================================
    // Return the reference to the last element.
    //=========================================================================
    inline value_type& back() {
      return m_blocks[m_cur_block_id][m_cur_block_filled - 1];
    }

    //=========================================================================
    // Return the reference to the last element.
    //=========================================================================
    inline const value_type& back() const {
      return m_blocks[m_cur_block_id][m_cur_block_filled - 1];
    }

    //=========================================================================
    // Access operator.
    //=========================================================================
    inline value_type& operator[] (const std::uint64_t i) {
      const std::uint64_t block_id = (i >> m_block_size_log);
      const std::uint64_t block_offset = (i & m_block_size_mask);
      return m_blocks[block_id][block_offset];
    }

    //=========================================================================
    // Access operator.
    //=========================================================================
    inline const value_type& operator[] (const std::uint64_t i) const {
      const std::uint64_t block_id = (i >> m_block_size_log);
      const std::uint64_t block_offset = (i & m_block_size_mask);
      return m_blocks[block_id][block_offset];
    }

    //=========================================================================
    // Remove all elements and deallocate.
    //=========================================================================
    void clear() {
      for (std::uint64_t block_id = 0;
          block_id < m_allocated_blocks; ++block_id)
        utils::deallocate(m_blocks[block_id]);

      m_size = 0;
      m_block_size_log = 0;
      m_block_size_mask = 0;
      m_block_size = 1;
      m_allocated_blocks = 1;
      m_cur_block_filled = 0;
      m_cur_block_id = 0;
      m_blocks[0] = utils::allocate_array<value_type>(m_block_size);
    }

    //=========================================================================
    // Write contents to given file.
    //=========================================================================
    void write_to_file(const std::string filename) const {
      std::FILE * const f = utils::file_open_nobuf(filename, "w");
      for (std::uint64_t block_id = 0; block_id < m_cur_block_id; ++block_id)
        utils::write_to_file(m_blocks[block_id], m_block_size, f);
      if (m_cur_block_filled > 0)
        utils::write_to_file(m_blocks[m_cur_block_id], m_cur_block_filled, f);
      std::fclose(f);
    }

    //=========================================================================
    // Return used RAM.
    //=========================================================================
    std::uint64_t ram_use() const {
      const std::uint64_t blocks_ram =
        sizeof(value_type) * m_allocated_blocks * m_block_size;
      const std::uint64_t total = blocks_ram;
      return total;
    }

    //=========================================================================
    // Reverse the vector.
    //=========================================================================
    void reverse() {
      const std::uint64_t all = size();
      const std::uint64_t half = all / 2;
      for (std::uint64_t i = 0; i < half; ++i)
        std::swap((*this)[i], (*this)[all - 1 - i]);
    }
};

#endif  // __SPACE_EFFICIENT_VECTOR_HPP_INCLUDED
