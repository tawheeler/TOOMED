#include "math_utils.hpp"

namespace common {

// ------------------------------------------------------------------------------------------------
int Wrap(int degrees, int lower) {
    int rev = 360;
    int upper = rev + lower;

    // Wrap to within [lower, lower + 360]
    while (degrees > upper) {
        degrees -= rev;
    }
    while (degrees < lower) {
        degrees += rev;
    }

    // If lower < 0, wrap within (lower, lower + 360], otherwise [lower, lower + 360)
    if (degrees == lower and lower < 0) {
        degrees = upper;
    } else if (degrees == upper and lower >= 0) {
        degrees = lower;
    }
    return degrees;
}

// ------------------------------------------------------------------------------------------------
double Wrap(double rads, double lower) {
    double two_pi = 2 * M_PI;
    double upper = lower + two_pi;

    // Wrap to within [lower, lower + 2*pi]
    while (rads > upper) {
        rads -= two_pi;
    }
    while (rads < lower) {
        rads += two_pi;
    }

    // If lower < 0, wrap within (lower, lower + 2*pi], otherwise [lower, lower + 2*pi)
    if (rads == lower and lower < 0) {
        rads = upper;
    } else if (rads == upper and lower >= 0) {
        rads = lower;
    }
    return rads;
}

// ------------------------------------------------------------------------------------------------
bool SameSign(double a, double b) { return (a >= 0) == (b >= 0); }

}  // namespace common