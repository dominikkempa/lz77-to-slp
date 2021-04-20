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


template<typename ValueType>
class space_efficient_vector {
  public:
    typedef ValueType value_type;

  private:
    static const std::uint64_t max_blocks = 16;

    value_type **m_blocks;
    std::uint64_t m_block_size_log;
    std::uint64_t m_block_size_mask;
    std::uint64_t m_block_size;
    std::uint64_t m_allocated_blocks;
    std::uint64_t m_cur_block_id;
    std::uint64_t m_cur_block_filled;
    std::uint64_t m_size;

  public:
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

    ~space_efficient_vector() {
      for (std::uint64_t block_id = 0;
          block_id < m_allocated_blocks; ++block_id)
        utils::deallocate(m_blocks[block_id]);
      utils::deallocate(m_blocks);
    }

    inline std::uint64_t size() const {
      return m_size;
    }

    inline bool empty() const {
      return m_size == 0;
    }

    inline void pop_back() {
      --m_size;
      --m_cur_block_filled;
      if (m_cur_block_filled == 0 && m_cur_block_id > 0) {
        --m_cur_block_id;
        m_cur_block_filled = m_block_size;
      }
    }

    void push_back(const value_type &value) {
      if (m_cur_block_filled == m_block_size &&
          m_cur_block_id + 1 == max_blocks) {
        std::uint64_t new_block_size = m_block_size * 2;
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

    inline value_type& back() {
      return m_blocks[m_cur_block_id][m_cur_block_filled - 1];
    }

    inline const value_type& back() const {
      return m_blocks[m_cur_block_id][m_cur_block_filled - 1];
    }

    inline value_type& operator[] (std::uint64_t i) {
      std::uint64_t block_id = (i >> m_block_size_log);
      std::uint64_t block_offset = (i & m_block_size_mask);
      return m_blocks[block_id][block_offset];
    }

    inline const value_type& operator[] (std::uint64_t i) const {
      std::uint64_t block_id = (i >> m_block_size_log);
      std::uint64_t block_offset = (i & m_block_size_mask);
      return m_blocks[block_id][block_offset];
    }

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

    void write_to_file(std::string filename) const {
      std::FILE *f = utils::file_open_nobuf(filename, "w");
      for (std::uint64_t block_id = 0; block_id < m_cur_block_id; ++block_id)
        utils::write_to_file(m_blocks[block_id], m_block_size, f);
      if (m_cur_block_filled > 0)
        utils::write_to_file(m_blocks[m_cur_block_id], m_cur_block_filled, f);
      std::fclose(f);
    }
};

#endif  // __SPACE_EFFICIENT_VECTOR_HPP_INCLUDED
