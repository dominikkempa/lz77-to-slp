#ifndef __UINT48_HPP_INCLUDED
#define __UINT48_HPP_INCLUDED

#include <cstdint>
#include <limits>


class uint48 {
  private:
    std::uint32_t low;
    std::uint16_t high;

  public:
    uint48() {}
    uint48(std::uint32_t l, std::uint16_t h) : low(l), high(h) {}
    uint48(const uint48& a) : low(a.low), high(a.high) {}
    uint48(const std::int32_t& a) : low(a), high(0) {}
    uint48(const std::uint32_t& a) : low(a), high(0) {}
    uint48(const std::uint64_t& a) :
      low(a & 0xFFFFFFFF), high((a >> 32) & 0xFFFF) {}
    uint48(const std::int64_t& a) :
      low(a & 0xFFFFFFFFL), high((a >> 32) & 0xFFFF) {}

    inline operator uint64_t() const {
      return (((std::uint64_t)high) << 32) | (std::uint64_t)low; }
    inline bool operator == (const uint48& b) const {
      return (low == b.low) && (high == b.high); }
    inline bool operator != (const uint48& b) const {
      return (low != b.low) || (high != b.high); }
} __attribute__((packed));

namespace std {

template<>
struct is_unsigned<uint48> {
  public:
    static const bool value = true;
};

template<>
class numeric_limits<uint48> {
  public:
    static uint48 min() {
      return uint48(std::numeric_limits<std::uint32_t>::min(),
          std::numeric_limits<std::uint16_t>::min());
    }

    static uint48 max() {
      return uint48(std::numeric_limits<std::uint32_t>::max(),
          std::numeric_limits<std::uint16_t>::max());
    }
};

}  // namespace std

#endif  // __UINT48_HPP_INCLUDED
