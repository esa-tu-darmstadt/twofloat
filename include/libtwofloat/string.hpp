#pragma once

#include <libtwofloat/twofloat.hpp>
#include <string>

namespace std {
template <typename T>
std::string to_string(twofloat::two<T> value) {
  return "{ " + to_string(value.h) + ", " + to_string(value.l) + " }";
}
}  // namespace std