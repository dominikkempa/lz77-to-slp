/**
 * @file    byte_rank.hpp
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

#ifndef __BYTE_RANK_HPP_INCLUDED
#define __BYTE_RANK_HPP_INCLUDED

#include <cstdint>
#include <algorithm>


// With the default block size (1024), the data
// structure uses a little over n bytes.
template<std::uint64_t block_size_log = 10,
  std::uint64_t sblock_size_log = 32>
struct byte_rank {
  static_assert(block_size_log <= sblock_size_log,
      "byte_rank requires: block_size_log <= sblock_size_log");

  private:
    const std::uint8_t *m_text;
    std::uint64_t m_length;
    std::uint64_t n_blocks;
    std::uint64_t n_sblocks;
    std::uint64_t *m_count;
    std::uint32_t *m_block_rank;
    std::uint64_t *m_sblock_rank;

    static const std::uint64_t block_size =
      ((std::uint64_t)1 << block_size_log);
    static const std::uint64_t sblock_size =
      ((std::uint64_t)1 << sblock_size_log);

  public:

    // Constructor.
    byte_rank(const std::uint8_t *text, std::uint64_t length) {

      // Compute basic parameters.
      m_text = text;
      m_length = length;
      n_blocks = (m_length + block_size - 1) / block_size;
      n_sblocks = (m_length + sblock_size - 1) / sblock_size;
      m_block_rank = new std::uint32_t[256 * n_blocks];
      m_sblock_rank = new std::uint64_t[256 * n_sblocks];
      m_count = new std::uint64_t[256];

      // Process all blocks.
      std::fill(m_count, m_count + 256, (std::uint64_t)0);
      std::uint64_t blocks_per_sblock = sblock_size / block_size;
      std::uint64_t sblock_id = 0;
      for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
        std::uint64_t block_beg = block_id * block_size;
        std::uint64_t block_end = std::min(m_length, block_beg + block_size);

        // Update sblock counts.
        if ((block_id % blocks_per_sblock) == 0) {
          for (std::uint64_t i = 0; i < 256; ++i)
            m_sblock_rank[sblock_id * 256 + i] = m_count[i];
          ++sblock_id;
        }

        // Update block counts.
        for (std::uint64_t i = 0; i < 256; ++i)
          m_block_rank[block_id * 256 + i] =
            m_count[i] - m_sblock_rank[(sblock_id - 1) * 256 + i];

        // Update symbol counts.
        for (std::uint64_t i = block_beg; i < block_end; ++i)
          ++m_count[m_text[i]];
      }
    }

    // Return number of occurrences of c in text[0..i).
    inline std::uint64_t query(std::uint64_t i, std::uint8_t c) const {

      // Handle special case.
      if (i >= m_length)
        return m_count[c];
      
      // Compute rank at block and superblock boundary.
      std::uint64_t block_id = (i >> block_size_log);
      std::uint64_t block_beg = (block_id << block_size_log);
      std::uint64_t sblock_id = (i >> sblock_size_log);
      std::uint64_t block_offset = (i & (block_size - 1));
      std::uint64_t sblock_rank = m_sblock_rank[(sblock_id << 8) + c];
      std::uint64_t block_rank = m_block_rank[(block_id << 8) + c];
      std::uint64_t ret = 0;

      // Compute rank within block.
      for (std::uint64_t j = 0; j < block_offset; ++j)
        if (m_text[block_beg + j] == c) ++ret;

      // Return result.
      return sblock_rank + block_rank + ret;
    }

    // Destructor.
    ~byte_rank() {
      delete[] m_block_rank;
      delete[] m_sblock_rank;
      delete[] m_count;
    }
};

#endif  // __BYTE_RANK_HPP_INCLUDED
