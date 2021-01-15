#ifndef __COMPUTE_LZ77_HPP_INCLUDED
#define __COMPUTE_LZ77_HPP_INCLUDED

#include <cstdint>
#include <vector>
#include <algorithm>


namespace compute_lz77 {

template<
  typename char_type,
  typename text_offset_type>
std::uint64_t parse_phrase(
    std::uint64_t, std::uint64_t,
    std::uint64_t, std::uint64_t,
    const char_type * const,
    std::vector<std::pair<text_offset_type, text_offset_type> > &);

template<
  typename char_type,
  typename text_offset_type>
void kkp2n(
    const char_type * const text,
    const std::uint64_t text_length,
    const text_offset_type * const sa,
    std::vector<std::pair<text_offset_type, text_offset_type> > &parsing) {

  // Handle special case.
  if (text_length == 0)
    return;

  // Compute the PSV array.
  text_offset_type *psv = new text_offset_type[text_length];
  {
    std::uint64_t prev_plus = 0;
    for (std::uint64_t i = 0; i < text_length; ++i) {
      const std::uint64_t cur = sa[i];
      while (prev_plus > 0 && prev_plus - 1 > cur) {
        const std::uint64_t j = prev_plus - 1;
        prev_plus = ((std::uint64_t)psv[j] == j) ?
          (std::uint64_t)0 : (std::uint64_t)psv[j] + 1;
      }
      psv[cur] = (prev_plus == 0) ? cur : prev_plus - 1;
      prev_plus = cur + 1;
    }
  }

  // Compute LZ77 phrases using the PSV values.
  // NSV values are computed on the fly from PSV.
  {
    parsing.push_back(
        std::make_pair(
          (text_offset_type)text[0],
          (text_offset_type)0));
    std::uint64_t next = 1;
    std::uint64_t list_head = 0;
    text_offset_type *invphi = psv;
    for (std::uint64_t i = 1; i < text_length; ++i) {
      const std::uint64_t psv_pos = psv[i];
      std::uint64_t nsv_pos = i;
      if (psv_pos == i) {
        nsv_pos = invphi[list_head];
        invphi[i] = nsv_pos;
        invphi[list_head] = i;
      } else {
        if (psv_pos != list_head) {
          nsv_pos = invphi[psv_pos];
          invphi[i] = invphi[psv_pos];
        } else {
          invphi[i] = invphi[list_head];
          list_head = i;
        }
        invphi[psv_pos] = i;
      }
      if (i == next)
        next = parse_phrase(i, text_length, psv_pos, nsv_pos, text, parsing);
    }
  }

  // Clean up.
  delete[] psv;
}

template<
  typename char_type,
  typename text_offset_type>
std::uint64_t parse_phrase(
    const std::uint64_t i,
    const std::uint64_t text_length,
    const std::uint64_t psv,
    const std::uint64_t nsv,
    const char_type * const text,
    std::vector<std::pair<text_offset_type, text_offset_type> > &parsing) {

  std::uint64_t pos = 0;
  std::uint64_t len = 0;

  if (nsv == i) {
    while (text[i + len] == text[psv + len])
      ++len;
    pos = psv;
  } else if (psv == i) {
    while (i + len < text_length && text[i + len] == text[nsv + len])
      ++len;
    pos = nsv;
  } else {
    while (text[psv + len] == text[nsv + len])
      ++len;
    if (text[i + len] == text[psv + len]) {
      ++len;
      while (text[i + len] == text[psv + len])
        ++len;
      pos = psv;
    } else {
      while (i + len < text_length && text[i + len] == text[nsv + len])
        ++len;
      pos = nsv;
    }
  }

  if (len == 0)
    parsing.push_back(
        std::make_pair(
          (text_offset_type)text[i],
          (text_offset_type)((std::uint64_t)0)));
  else parsing.push_back(
      std::make_pair(
        (text_offset_type)pos,
        (text_offset_type)len));

  return i + std::max((std::uint64_t)1, len);
}

}  // namespace compute_lz77

#endif  // __COMPUTE_LZ77_HPP_INCLUDED
