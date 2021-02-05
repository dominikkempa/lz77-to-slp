#ifndef __NAIVE_COMPUTE_SA_HPP_INCLUDED
#define __NAIVE_COMPUTE_SA_HPP_INCLUDED

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>


template<typename char_type>
class substring {
  public:
    typedef substring<char_type> substring_type;

    substring() {}
    substring(
        const char_type * const text,
        const std::uint64_t beg,
        const std::uint64_t length,
        const std::uint64_t text_length) {

      m_beg = beg;
      m_text_length = text_length;
      for (std::uint64_t j = 0; j < length; ++j)
        m_data.push_back(text[beg + j]);
    }

    inline bool operator < (const substring_type &s) const {
      std::uint64_t lcp = 0;
      while (m_beg + lcp < m_text_length &&
          s.m_beg + lcp < m_text_length &&
          m_data[lcp] == s.m_data[lcp])
        ++lcp;

      return (m_beg + lcp == m_text_length ||
          (s.m_beg + lcp < m_text_length &&
           (std::uint64_t)m_data[lcp] < (std::uint64_t)s.m_data[lcp]));
    }

    std::uint64_t m_beg;
    std::uint64_t m_text_length;
    std::vector<char_type> m_data;
};

template<
  typename char_type,
  typename text_offset_type>
void naive_compute_sa(
    const char_type * const text,
    const std::uint64_t text_length,
    text_offset_type * const sa) {

  typedef substring<char_type> substring_type;
  std::vector<substring_type> substrings;

  for (std::uint64_t i = 0; i < text_length; ++i)
    substrings.push_back(
        substring_type(text, i, text_length - i, text_length));

  std::sort(substrings.begin(), substrings.end());
  for (std::uint64_t i = 0; i < text_length; ++i)
    sa[i] = substrings[i].m_beg;
}

#endif  // __NAIVE_COMPUTE_SA_HPP_INCLUDED
