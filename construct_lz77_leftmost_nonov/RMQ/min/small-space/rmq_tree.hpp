/**
 * @file    rmq_tree.hpp
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

#ifndef __RMQ_TREE_HPP_INCLUDED
#define __RMQ_TREE_HPP_INCLUDED

#include <cstdint>
#include <algorithm>


// With the default block size (256), the data structure uses
// (assuming ValueType = std::uint64_t) at most n bits. For
// smaller ValueType it uses even less space.
template<typename ValueType,
  std::uint64_t BlockSize = 256>
struct rmq_tree {
  public:
    typedef ValueType value_type;

  private:
    std::uint64_t m_size;
    std::uint64_t m_blocks;
    std::uint64_t m_leaves2;

    const value_type *m_tab;
    value_type *m_data;
    std::uint64_t *m_pos;

    const std::uint64_t block_size = BlockSize;

  public:

    // Constructor.
    rmq_tree(const value_type *tab, std::uint64_t size) {
      m_size = size;
      m_tab = tab;
      m_blocks = (size + block_size - 1) / block_size;

      if (m_blocks > 1) {

        // Compute the number of leaves.
        m_leaves2 = 1;
        while (2 * m_leaves2 < m_blocks)
          m_leaves2 <<= 1;

        // Allocate data.
        m_data = new value_type[m_leaves2 * 2];
        m_pos = new std::uint64_t[m_leaves2 * 2];

        // Encode bottom level of the tree.
        for (std::uint64_t i = 0; i < m_leaves2; ++i) {
          std::uint64_t left_val = 0;
          std::uint64_t right_val = 0;
          std::uint64_t left_pos = m_size;
          std::uint64_t right_pos = m_size;

          // Check left child.
          if ((i << 1) < m_blocks) {
            std::uint64_t block_id = (i << 1);
            std::uint64_t block_beg = block_id * block_size;
            std::uint64_t block_end =
              std::min(m_size, block_beg + block_size);
            left_pos = block_beg;
            left_val = m_tab[left_pos];
            for (std::uint64_t j = block_beg + 1; j < block_end; ++j) {
              if ((std::uint64_t)m_tab[j] < left_val) {
                left_val = m_tab[j];
                left_pos = j;
              }
            }
          }

          // Check right child.
          if ((i << 1) + 1 < m_blocks) {
            std::uint64_t block_id = (i << 1) + 1;
            std::uint64_t block_beg = block_id * block_size;
            std::uint64_t block_end =
              std::min(m_size, block_beg + block_size);
            right_pos = block_beg;
            right_val = m_tab[right_pos];
            for (std::uint64_t j = block_beg + 1; j < block_end; ++j) {
              if ((std::uint64_t)m_tab[j] < right_val) {
                right_val = m_tab[j];
                right_pos = j;
              }
            }
          }

          // Store the answer.
          if (left_val < right_val) {
            m_data[m_leaves2 + i] = left_val;
            m_pos[m_leaves2 + i] = left_pos;
          } else {
            m_data[m_leaves2 + i] = right_val;
            m_pos[m_leaves2 + i] = right_pos;
          }
        }

        // Encode remaining levels of the tree.
        for (std::uint64_t i = m_leaves2 - 1; i > 0; --i) {
          std::uint64_t left = m_data[i << 1];
          std::uint64_t right = m_data[(i << 1) + 1];
          if (left < right) {
            m_data[i] = left;
            m_pos[i] = m_pos[i << 1];
          } else {
            m_data[i] = right;
            m_pos[i] = m_pos[(i << 1) + 1];
          }
        }
      } else m_data = NULL;
    }

    // Destructor.
    ~rmq_tree() {
      if (m_data != NULL) {
        delete[] m_data;
        delete[] m_pos;
      }
    }

    // Return the boolean value telling whether there is any item
    // in the range [beg..end) that is < than given threshold.
    inline bool less(std::uint64_t beg, std::uint64_t end,
        std::uint64_t threshold) const {

      // Handle special cases.
      beg = std::min(beg, m_size);
      end = std::min(end, m_size);
      if (beg >= end)
        return false;
      if (beg + 1 == end)
        return (std::uint64_t)m_tab[beg] < threshold;

      // Compute block IDs.
      std::uint64_t left_block_id = beg / block_size;
      std::uint64_t right_block_id = (end - 1) / block_size;

      // Handle another special case.
      if (left_block_id == right_block_id) {
        for (std::uint64_t j = beg; j < end; ++j)
          if ((std::uint64_t)m_tab[j] < threshold)
            return true;

        return false;
      }

      // Check leftmost and rightmost blocks.
      {

        // Check leftmost block.
        std::uint64_t left_block_beg = left_block_id * block_size;
        std::uint64_t left_block_end = std::min(m_size,
            left_block_beg + block_size);
        left_block_beg = std::max(left_block_beg, beg);
        for (std::uint64_t j = left_block_beg; j < left_block_end; ++j)
          if ((std::uint64_t)m_tab[j] < threshold)
            return true;

        // Check rightmost block.
        std::uint64_t right_block_beg = right_block_id * block_size;
        std::uint64_t right_block_end = std::min(m_size,
            right_block_beg + block_size);
        right_block_end = std::min(right_block_end, end);
        for (std::uint64_t j = right_block_beg; j < right_block_end; ++j)
          if ((std::uint64_t)m_tab[j] < threshold)
            return true;
      }

      // Prepare pointers in the tree.
      std::uint64_t off = 2 * m_leaves2;
      left_block_id += off;
      right_block_id += off;

      // Traverse un in the tree until pointers meet.
      while (left_block_id + 1 != right_block_id) {

        // Process right sibling of the left pointer.
        if (!(left_block_id & 1)) {
          if (left_block_id >= off) {
            std::uint64_t block_id = (left_block_id + 1) - off;
            std::uint64_t block_beg = block_id * block_size;
            std::uint64_t block_end = std::min(m_size,
                block_beg + block_size);
            for (std::uint64_t j = block_beg; j < block_end; ++j)
              if ((std::uint64_t)m_tab[j] < threshold)
                return true;
          } else {
            std::uint64_t val = m_data[left_block_id + 1];
            if (val < threshold)
              return true;
          }
        }

        // Process left sibling of the right pointer.
        if (right_block_id & 1) {
          if (right_block_id >= off) {
            std::uint64_t block_id = (right_block_id - 1) - off;
            std::uint64_t block_beg = block_id * block_size;
            std::uint64_t block_end = std::min(m_size,
                block_beg + block_size);
            for (std::uint64_t j = block_beg; j < block_end; ++j)
              if ((std::uint64_t)m_tab[j] < threshold)
                return true;
          } else {
            std::uint64_t val = m_data[right_block_id - 1];
            if (val < threshold)
              return true;
          }
        }

        // Update pointers.
        left_block_id >>= 1;
        right_block_id >>= 1;
      }

      // No value was found.
      return false;
    }

    // Return position of min in the range [beg..end).
    inline std::uint64_t rmq(std::uint64_t beg, std::uint64_t end) const {

      // Handle special cases.
      beg = std::min(beg, m_size);
      end = std::min(end, m_size);
      if (beg >= end)
        return m_size;
      if (beg + 1 == end)
        return beg;

      // Compute block IDs.
      std::uint64_t left_block_id = beg / block_size;
      std::uint64_t right_block_id = (end - 1) / block_size;
      std::uint64_t ret_pos = m_size;
      std::uint64_t ret_val = 0;

      // Handle another special case.
      if (left_block_id == right_block_id) {
        ret_val = m_tab[beg];
        ret_pos = beg;
        for (std::uint64_t j = beg + 1; j < end; ++j) {
          if ((std::uint64_t)m_tab[j] < ret_val) {
            ret_val = m_tab[j];
            ret_pos = j;
          }
        }
        return ret_pos;
      }


      // Check leftmost and rightmost blocks.
      {

        // Check leftmost block.
        std::uint64_t left_block_beg = left_block_id * block_size;
        std::uint64_t left_block_end = std::min(m_size,
            left_block_beg + block_size);
        left_block_beg = std::max(left_block_beg, beg);
        ret_val = m_tab[left_block_beg];
        ret_pos = left_block_beg;
        for (std::uint64_t j = left_block_beg + 1; j < left_block_end; ++j) {
          if ((std::uint64_t)m_tab[j] < ret_val) {
            ret_val = m_tab[j];
            ret_pos = j;
          }
        }

        // Check rightmost block.
        std::uint64_t right_block_beg = right_block_id * block_size;
        std::uint64_t right_block_end = std::min(m_size,
            right_block_beg + block_size);
        right_block_end = std::min(right_block_end, end);
        for (std::uint64_t j = right_block_beg; j < right_block_end; ++j) {
          if ((std::uint64_t)m_tab[j] < ret_val) {
            ret_val = m_tab[j];
            ret_pos = j;
          }
        }
      }

      // Prepare pointers in the tree.
      std::uint64_t off = 2 * m_leaves2;
      left_block_id += off;
      right_block_id += off;

      // Traverse up in the tree until pointers meet.
      while (left_block_id + 1 != right_block_id) {

        // Process right sibling of the left pointer.
        if (!(left_block_id & 1)) {
          if (left_block_id >= off) {
            std::uint64_t block_id = (left_block_id + 1) - off;
            std::uint64_t block_beg = block_id * block_size;
            std::uint64_t block_end = std::min(m_size,
                block_beg + block_size);
            for (std::uint64_t j = block_beg; j < block_end; ++j) {
              if ((std::uint64_t)m_tab[j] < ret_val) {
                ret_val = m_tab[j];
                ret_pos = j;
              }
            }
          } else {
            if ((std::uint64_t)m_data[left_block_id + 1] < ret_val) {
              ret_val = m_data[left_block_id + 1];
              ret_pos = m_pos[left_block_id + 1];
            }
          }
        }

        // Process left sibling of the right pointer.
        if (right_block_id & 1) {
          if (right_block_id >= off) {
            std::uint64_t block_id = (right_block_id - 1) - off;
            std::uint64_t block_beg = block_id * block_size;
            std::uint64_t block_end = std::min(m_size,
                block_beg + block_size);
            for (std::uint64_t j = block_beg; j < block_end; ++j) {
              if ((std::uint64_t)m_tab[j] < ret_val) {
                ret_val = m_tab[j];
                ret_pos = j;
              }
            }
          } else {
            if ((std::uint64_t)m_data[right_block_id - 1] < ret_val) {
              ret_val = m_data[right_block_id - 1];
              ret_pos = m_pos[right_block_id - 1];
            }
          }
        }

        // Update pointers.
        left_block_id >>= 1;
        right_block_id >>= 1;
      }

      // Return the position of minimum.
      return ret_pos;
    }
};

#endif  // __RMQ_TREE_HPP_INCLUDED
