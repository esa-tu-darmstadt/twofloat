#pragma once

#include <libtwofloat/twofloat.hpp>
#include <limits>
#include <type_traits>

namespace std {
/// \brief Specialization of std::numeric_limits for twofloat::twofloat.
template <typename T>
class numeric_limits<::twofloat::two<T>> : public numeric_limits<T> {
  // Make sure that T is a floating point type
  static_assert(
      is_floating_point<T>::value,
      "twofloat::two<T> can only be instantiated with a floating point "
      "type.");

 public:
  /// \brief number of radix digits that can be represented without change
  static constexpr int digits = numeric_limits<T>::digits * 2;

  /// \brief number of decimal digits that can be represented without change
  /// \details Calculated as floor(digits * log10(2))
  static constexpr int digits10 = (int)((digits - 1) * 0.30102999566);

  /// \brief returns the difference between 1.0 and the next representable value
  /// of the given floating-point type
  static constexpr T epsilon() noexcept {
    return numeric_limits<T>::epsilon() * numeric_limits<T>::epsilon();
  }
};
}  // namespace std