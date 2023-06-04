#pragma once

// Common math utilities

#include <algorithm>
#include <array>
#include <cmath>
#include <ostream>

namespace common {

// Wrap an angle in radians to lie in the range [lower, 2*PI + lower) if lower >= 0 or
// (lower, 2*PI + lower] otherwise
double Wrap(double rads, double lower = -M_PI);

template <typename T>
bool IsNear(const T& a, const T& b, const T& tolerance = std::numeric_limits<T>::epsilon()) {
    // NOTE: <= so that we can get exact equality for integer types.
    return std::abs(a - b) <= tolerance;
}

template <typename T>
bool IsZero(const T& a, const T& tolerance = std::numeric_limits<T>::epsilon()) {
    return IsNear(a, T{}, tolerance);
}

// Check if all elements in array are equal to some value
// Use caution when comparing floating point types
template <typename T, size_t N>
bool AllEqual(const std::array<T, N>& arr, T val) {
    return std::all_of(arr.begin(), arr.end(), [val](T x) { return x == val; });
}

// Check if all elements in array are close to some value
template <typename T, size_t N>
bool AllNear(const std::array<T, N>& arr, T val, T tolerance) {
    return std::all_of(arr.begin(), arr.end(),
                       [val, tolerance](T x) { return IsNear(x, val, tolerance); });
}

// Get the clockwise distance to get from angle a [rad] to b [rad]
template <typename T>
T ClockwiseAngleDistance(T a, T b) {
    // Eq 11.8 in Small Unmanned Aircraft by Beard and McLain
    // Note that they use positive-north whereas we use positive x-axis
    a = Wrap(a, 0.0);  // [0,2pi]
    b = Wrap(b, 0.0);  // [0,2pi]
    return Wrap(2 * M_PI - b + a, 0.0);
}

// Counter clockwise distance to get from angle a [rad] to b [rad]
template <typename T>
T CounterClockwiseAngleDistance(T a, T b) {
    // Eq 11.7 in Small Unmanned Aircraft by Beard and McLain
    // Note that they use positive-north whereas we use positive x-axis
    a = Wrap(a, 0.0);  // [0,2pi]
    b = Wrap(b, 0.0);  // [0,2pi]
    return Wrap(2 * M_PI + b - a, 0.0);
}

// Get the minimum-magnitude δ such that a + δ = mod(b, 2π)
// a - The first angle [rad]
// b - The second angle [rad]
template <typename T>
T SignedAngleDist(T a, T b) {
    return std::atan2(std::sin(b - a), std::cos(b - a));
}

// Get the unsigned angle distance.
template <typename T>
T AngleDist(T a, T b) {
    return std::abs(SignedAngleDist(a, b));
}

// Both angles are in radians
template <typename T>
bool IsAngleNear(const T& a, const T& b, const T& tolerance = std::numeric_limits<T>::epsilon()) {
    return AngleDist(a, b) <= tolerance;
}

// Whether the angle a [rad] is within a tolerance distance of angle b [rad],
// or b's reciprocal.
template <typename T>
bool IsNearToAngleOrReciprocal(const T& a, const T& b,
                               const T& tolerance = std::numeric_limits<T>::epsilon()) {
    T dyaw = AngleDist(a, b);
    return std::min(dyaw, static_cast<T>(M_PI) - dyaw) <= tolerance;
}

// True if angles a and b are within 90 deg of each other
template <typename T>
bool AreInSameHemisphere(const T& a, const T& b) {
    return AngleDist(a, b) <= M_PI_2;
}

// Clamp a value to a target range.
// NOTE: std::clamp not defined until C++20
template <typename T>
T Clamp(T value, T lowerbound, T upperbound) {
    if (value > upperbound) {
        value = upperbound;
    } else if (value < lowerbound) {
        value = lowerbound;
    }
    return value;
}

// A second version of clamp for which lowerbound = -upperbound
template <typename T>
T Clamp(T value, T upperbound) {
    return Clamp(value, -upperbound, upperbound);
}

// Return the minimum of two values in an absolute sense
// AbsMin(1, 2) = 1
// AbsMin(-1, -2) = -1
template <typename T>
T AbsMin(const T& a, const T& b) {
    if (std::abs(a) <= std::abs(b)) {
        return a;
    }
    return b;
}

// Move from current value toward target at the given rate. Return the new value and an int
// indicating how the value was changed.
// value  - Current value [units]
// target - Target value [units]
// rate   - Rate of change [units/s]
// dt     - Elapsed time [s]
template <typename T>
std::pair<T, int> RateLimit(T value, T target, T rate, T dt) {
    int changed = 0;
    if (target > (value + rate * dt)) {
        target = value + rate * dt;
        changed = 1;
    } else if (target < (value - rate * dt)) {
        target = value - rate * dt;
        changed = -1;
    }
    return std::make_pair(target, changed);
}

// Move from current angle toward target angle at the given rate, compensating for wraparound
// value  - Current value [rad]
// target - Target value [rad]
// rate   - Rate of change [rad/s]
// dt     - Elapsed time [s]
template <typename T>
T RateLimitAngle(T value, T target, T rate, T dt) {
    T signed_error = SignedAngleDist(value, target);
    T angle_step = Clamp(signed_error, rate * dt);
    return value + angle_step;
}

// Return true if a and b have the same sign
bool SameSign(double a, double b);

// Calculate the determinant of the 2x2 matrix [a b]
//                                             [c d]
template <typename T>
T Det(T a, T b, T c, T d) {
    return a * d - b * c;
}

// Calculate the determinant of the 3x3 matrix [a b c]
//                                             [d e f]
//                                             [g h i]
template <typename T>
T Det(T a, T b, T c, T d, T e, T f, T g, T h, T i) {
    return a * Det(e, f, h, i) - b * Det(d, f, g, i) + c * Det(d, e, g, h);
}

// Calculate the determinant of the 4x4 matrix [a b c d]
//                                             [e f g h]
//                                             [i j k l]
//                                             [m n o p]
template <typename T>
T Det(T a, T b, T c, T d, T e, T f, T g, T h, T i, T j, T k, T l, T m, T n, T o, T p) {
    // clang-format off
    return a * Det(f, g, h,
                   j, k, l,
                   n, o, p)
         - b * Det(e, g, h,
                   i, k, l,
                   m, o, p)
         + c * Det(e, f, h,
                   i, j, l,
                   m, n, p)
         - d * Det(e, f, g,
                   i, j, k,
                   m, n, o);
    // clang-format on
}

}  // namespace common