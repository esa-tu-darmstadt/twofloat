#pragma once
#include <ipudef.h>

#include <limits>
#include <type_traits>

namespace std {
// Specializations for SIMD types
template <>
struct is_floating_point<float2> : std::true_type {};
template <>
struct is_floating_point<float4> : std::true_type {};

template <>
class numeric_limits<float2> : public std::numeric_limits<float> {
 public:
  static constexpr float2 max() noexcept {
    return {std::numeric_limits<float>::max(),
            std::numeric_limits<float>::max()};
  }
};
template <>
class numeric_limits<float4> : public std::numeric_limits<float> {};
}  // namespace std
