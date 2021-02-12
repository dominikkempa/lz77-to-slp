#ifndef __LARGE_RMQ_TREE_HPP_INCLUDED
#define __LARGE_RMQ_TREE_HPP_INCLUDED

#include <cstdint>
#include <algorithm>


template<typename ValueType>
struct large_rmq_tree {
  public:
    typedef ValueType value_type;

  private:
    std::uint64_t m_size;
    std::uint64_t m_leaves2;

    const value_type *m_tab;
    value_type *m_data;
    std::uint64_t *m_pos;

  public:
    large_rmq_tree(const value_type *tab, std::uint64_t size) {
      m_size = size;
      m_tab = tab;

      if (m_size > 1) {
        m_leaves2 = 1;
        while (2 * m_leaves2 < m_size)
          m_leaves2 <<= 1;
        m_data = new value_type[m_leaves2 * 2];
        m_pos = new std::uint64_t[m_leaves2 * 2];

        for (std::uint64_t i = 0; i < m_leaves2; ++i) {
          std::uint64_t left = 0;
          std::uint64_t right = 0;
          if ((i << 1) < m_size) left = m_tab[i << 1];
          if ((i << 1) + 1 < m_size) right = m_tab[(i << 1) + 1];
          if (left > right) {
            m_data[m_leaves2 + i] = left;
            m_pos[m_leaves2 + i] = (i << 1);
          } else {
            m_data[m_leaves2 + i] = right;
            m_pos[m_leaves2 + i] = (i << 1) + 1;
          }
        }

        for (std::uint64_t i = m_leaves2 - 1; i > 0; --i) {
          std::uint64_t left = m_data[i << 1];
          std::uint64_t right = m_data[(i << 1) + 1];
          if (left > right) {
            m_data[i] = left;
            m_pos[i] = m_pos[i << 1];
          } else {
            m_data[i] = right;
            m_pos[i] = m_pos[(i << 1) + 1];
          }
        }
      } else m_data = NULL;
    }

    ~large_rmq_tree() {
      if (m_data != NULL) {
        delete[] m_data;
        delete[] m_pos;
      }
    }

    // Return the boolean value telling whether there is any item
    // in the range [beg..end) that is >= than given threshold.
    inline bool geq(
        std::uint64_t beg,
        std::uint64_t end,
        std::uint64_t threshold) const {

      beg = std::min(beg, m_size);
      end = std::min(end, m_size);

      if (beg >= end) return false;
      else if (m_size == 1) {
        if (beg != 0 || end != 1) return false;
        else return (std::uint64_t)m_tab[0] >= threshold;
      } else if (beg + 1 == end)
        return (std::uint64_t)m_tab[beg] >= threshold;

      std::uint64_t off = 2 * m_leaves2;
      std::uint64_t left = off + beg;
      std::uint64_t right = off + end - 1;
      std::uint64_t ret = std::max(
          (std::uint64_t)m_tab[left - off],
          (std::uint64_t)m_tab[right - off]);

      if (ret >= threshold)
        return true;

      static const std::uint64_t small_range = 512;
      if (right - left <= small_range) {
        for (std::uint64_t j = beg + 1; j + 1 < end; ++j)
          if ((std::uint64_t)m_tab[j] >= threshold) return true;
        return false;
      }

      while (left + 1 != right) {
        if (!(left & 1)) {
          if (left >= off)
            ret = std::max(ret, (std::uint64_t)m_tab[left - off + 1]);
          else ret = std::max(ret, (std::uint64_t)m_data[left + 1]);
        }

        if (right & 1) {
          if (right >= off)
            ret = std::max(ret, (std::uint64_t)m_tab[right - off - 1]);
          else ret = std::max(ret, (std::uint64_t)m_data[right - 1]);
        }

        if (ret >= threshold)
          return true;

        left >>= 1;
        right >>= 1;
      }

      return false;
    }

    // Return position of max in the range [beg..end).
    inline std::uint64_t rmq(std::uint64_t beg, std::uint64_t end) const {
      beg = std::min(beg, m_size);
      end = std::min(end, m_size);

      if (beg >= end)
        return m_size;
      else if (m_size == 1) {
        if (beg != 0 || end != 1) return m_size;
        else return 0;
      } else if (beg + 1 == end)
        return beg;

      std::uint64_t off = 2 * m_leaves2;
      std::uint64_t left = off + beg;
      std::uint64_t right = off + end - 1;
      std::uint64_t ret_pos = m_size;
      std::uint64_t ret_val = 0;
      if ((std::uint64_t)m_tab[left - off] >
          (std::uint64_t)m_tab[right - off]) {
        ret_val = m_tab[left - off];
        ret_pos = beg;
      } else {
        ret_val = m_tab[right - off];
        ret_pos = end - 1;
      }

      while (left + 1 != right) {
        if (!(left & 1)) {
          if (left >= off) {
            if ((std::uint64_t)m_tab[left - off + 1] > ret_val) {
              ret_val = m_tab[left - off + 1];
              ret_pos = left - off + 1;
            }
          } else {
            if ((std::uint64_t)m_data[left + 1] > ret_val) {
              ret_val = m_data[left + 1];
              ret_pos = m_pos[left + 1];
            }
          }
        }

        if (right & 1) {
          if (right >= off) {
            if ((std::uint64_t)m_tab[right - off - 1] > ret_val) {
              ret_val = m_tab[right - off - 1];
              ret_pos = right - off - 1;
            }
          } else {
            if ((std::uint64_t)m_data[right - 1] > ret_val) {
              ret_val = m_data[right - 1];
              ret_pos = m_pos[right - 1];
            }
          }
        }

        left >>= 1;
        right >>= 1;
      }

      return ret_pos;
    }
};

#endif  // __LARGE_RMQ_TREE_HPP_INCLUDED
