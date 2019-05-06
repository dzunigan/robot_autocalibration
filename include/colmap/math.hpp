// Minimal required math functions for COLMAP

#ifndef COLMAP_SRC_UTIL_MATH_H_
#define COLMAP_SRC_UTIL_MATH_H_

#include <algorithm>

namespace colmap {

bool IsNaN(const float x) { return x != x; }
bool IsNaN(const double x) { return x != x; }

bool IsInf(const float x) { return !IsNaN(x) && IsNaN(x - x); }
bool IsInf(const double x) { return !IsNaN(x) && IsNaN(x - x); }

template <typename T>
T Clip(const T& value, const T& low, const T& high) {
  return std::max(low, std::min(value, high));
}

constexpr inline float DegToRad(const float deg) {
  return deg * 0.0174532925199432954743716805978692718781530857086181640625f;
}

constexpr inline double DegToRad(const double deg) {
  return deg * 0.0174532925199432954743716805978692718781530857086181640625;
}

constexpr inline float RadToDeg(const float rad) {
  return rad * 57.29577951308232286464772187173366546630859375f;
}

constexpr inline double RadToDeg(const double rad) {
  return rad * 57.29577951308232286464772187173366546630859375;
}

}  // namespace colmap

#endif  // COLMAP_SRC_UTIL_MATH_H_
