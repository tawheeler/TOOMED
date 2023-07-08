#pragma once

#include <cmath>
#include <vector>

#include "math_utils.hpp"

namespace common {

struct Vec2f {
    float x;
    float y;

    // constructors
    Vec2f(float x_, float y_) : x(x_), y(y_) {}
    Vec2f() : Vec2f(0.0f, 0.0f) {}

    // operators
    bool operator==(const Vec2f& rhs) const { return x == rhs.x && y == rhs.y; }
    bool operator!=(const Vec2f& rhs) const { return !(*this == rhs); }
    Vec2f operator*(float v) const { return Vec2f({v * x, v * y}); }
    Vec2f operator/(float v) const { return Vec2f({x / v, y / v}); }
    Vec2f operator+(const Vec2f& rhs) const { return Vec2f({x + rhs.x, y + rhs.y}); }
    Vec2f operator-(const Vec2f& rhs) const { return Vec2f({x - rhs.x, y - rhs.y}); }
    Vec2f operator+(float v) const { return Vec2f({x + v, y + v}); }
    Vec2f operator-(float v) const { return Vec2f({x - v, y - v}); }
    Vec2f& operator+=(const Vec2f& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }
    Vec2f& operator-=(const Vec2f& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }
    Vec2f& operator*=(float f) {
        x *= f;
        y *= f;
        return *this;
    }
    Vec2f& operator/=(float f) {
        x /= f;
        y /= f;
        return *this;
    }
};

// Operator overload for left-multiply.
Vec2f operator*(float scalar, const Vec2f& vec);

// An integer version of the above type.
struct Vec2i {
    long x;
    long y;

    // constructors
    Vec2i(long x_, long y_) : x(x_), y(y_) {}
    Vec2i() : Vec2i(0, 0) {}

    // operators
    bool operator==(const Vec2i& rhs) const { return x == rhs.x && y == rhs.y; }
    bool operator!=(const Vec2i& rhs) const { return !(*this == rhs); }
    Vec2i operator*(long v) const { return Vec2i({v * x, v * y}); }
    Vec2i operator+(const Vec2i& rhs) const { return Vec2i({x + rhs.x, y + rhs.y}); }
    Vec2i operator-(const Vec2i& rhs) const { return Vec2i({x - rhs.x, y - rhs.y}); }
    Vec2i operator+(long v) const { return Vec2i({x + v, y + v}); }
    Vec2i operator-(long v) const { return Vec2i({x - v, y - v}); }
    Vec2i& operator+=(const Vec2i& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }
    Vec2i& operator-=(const Vec2i& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    // The lexicographic order for Vec2i is by minimum y value, then minimum x value.
    bool operator<(const Vec2i& rhs) const { return y < rhs.y || (y == rhs.y && x < rhs.x); }
    bool operator<=(const Vec2i& rhs) const { return y < rhs.y || (y == rhs.y && x <= rhs.x); }
};

// Construct a fixed-point represented Vec of the given resolution.
// For example, rather than representing a position in meters with floating point, we can represent
// it in mm in integers with a resolution factor of 1000 mm/m.
Vec2i ToFixedPoint(const Vec2f& v, float resolution);

// Obtain the floating-point represented Vec from the fixed-point version
Vec2f ToFloatingPoint(const Vec2i& v, float resolution);

// Convert a sequence of points in floating point to fixed point
std::vector<Vec2i> ToFixedPoint(const std::vector<Vec2f>& vertices, float resolution);

// Convert a sequence of points in fixed point to floating point
std::vector<Vec2f> ToFloatingPoint(const std::vector<Vec2i>& vertices, float resolution);

template <typename T>
const T& GetCyclicVertex(const std::vector<T>& polygon, size_t i) {
    return polygon.at(i % polygon.size());
}

// Vector 2-norm
float Norm(const Vec2f& v);
float Norm(const Vec2i& v);
Vec2f Normalize(const Vec2f& v);

// Squared vector 2-norm
float SqNorm(const Vec2f& v);
long long SqNorm(const Vec2i& v);

// Compute the L2 distance between a and b
float Dist(const Vec2f& a, const Vec2f& b);
float Dist(const Vec2i& a, const Vec2i& b);

// Dot product
float Dot(const Vec2f& a, const Vec2f& b);

// Vector cross product.
float Cross(const Vec2f& a, const Vec2f& b);
long long Cross(const Vec2i& a, const Vec2i& b);

// Vector projection
Vec2f VectorProjection(const Vec2f& a, const Vec2f& b);

// Clamp a vector by magnitude
Vec2f ClampNorm(const Vec2f& v, float threshold);

// Get a + (b-a)*interpolant
Vec2f Interp(const Vec2f& a, const Vec2f& b, float interpolant);

// Get a + (b-a)*interpolant, but note that because we are using fixed-point, we
// round the result to the nearest fixed-point value.
Vec2i Interp(const Vec2i& a, const Vec2i& b, float interpolant);

// Rotate the given angle right by 90 deg.
Vec2f Rotr(const Vec2f& a);

// Determine whether the given points are in right-hand order (CCW).
// Returns a value > 0 if right-hand (CCW).
// Returns a value < 0 if left-hand (CW).
// Returns zero if the points are colinear.
float GetRightHandedness(const Vec2f& a, const Vec2f& b, const Vec2f& c);

// Get the distance of point P to the line that goes through A and B.
// Note that A and B cannot be colocated.
float GetDistanceToLine(const Vec2f& p, const Vec2f& a, const Vec2f& b);

// Get the distance of point P to the line segment AB.
// Note that A and B cannot be colocated.
float GetDistanceToLineSegment(const Vec2f& p, const Vec2f& a, const Vec2f& b);

enum class Intersection : uint8_t {
    UNDEFINED = 0,
    NONE = 1,
    AT_POINT = 2,
    OVERLAP = 3,
};

// Check for an intersection P between two lines, AB and CD, returning the interpolants
// s and t such that P = A + (B-A)*s = C + (D-C)*t.
// Note that if AB and CD are line segments, we only get an intersection if s and t in [0,1].
struct LineIntersectionResult {
    Intersection intersection;
    float s;  // Intersection P = A + (B-A)*s
    float t;  // Intersection P = C + (D-C)*t
};
LineIntersectionResult CalcLineIntersection(const Vec2f& a, const Vec2f& b, const Vec2f& c,
                                            const Vec2f& d);
LineIntersectionResult CalcLineIntersection(const Vec2i& a, const Vec2i& b, const Vec2i& c,
                                            const Vec2i& d);

// Determine whether the given interpolant for a line segment produces a point inside, at an edge,
// or outside the line segment. Our point P is A + (B-A)*s for segment endpoints A and B.
// We assume A is not coincident with B.
// Returns 1 if P is inside segment AB, 0 if P is A or B, and -1 if P is outside of segment AB.
int GetLineSegmentContainment(float s);

// Determine whether the given line intersection result, given by its interpolations, is also a line
// segment intersection.
// Returns 1 if a point intersection occurs inside each segment,
//         0 if a point intersection occurs inside or on the border of each segment,
//        -1 if there is no point intersection.
int GetLineSegmentContainment(float s, float t);

// Determine whether the given point is inside (or on the border of) the triangle defined by the
// three vertices.
// Returns 1 if inside the triangle, 0 if on its perimeter, and -1 if outside.
int GetTriangleContainment(const Vec2f& p, const Vec2f& a, const Vec2f& b, const Vec2f& c);

// Determine whether the given point is inside the circle inscribed on the three vertices.
// Returns 1 if inside the circle,
//         0 if on the circle,
//    and -1 if outside the circle.
int GetCircleContainment(const Vec2f& p, const Vec2f& a, const Vec2f& b, const Vec2f& c);

// Get the lexicographically lowest vertex index for the given polygon.
size_t GetLowestVertexIndex(const std::vector<Vec2i>& polygon);

// Compute the Minkowski sum of two convex polygons P and Q, whose vertices are assumed to be in CCW
// order.
std::vector<Vec2i> MinkowskiSum(const std::vector<Vec2i>& P, const std::vector<Vec2i>& Q);

// Compute the union of two polygons P & Q, whose vertices are assumed to be in CCW order.
// Returns an empty vector if no joining occurred.
std::vector<Vec2i> Union(const std::vector<Vec2i>& P, const std::vector<Vec2i>& Q);

// Get the centroid (arithmetic mean of its points) of a polygon.
Vec2f GetCentroid(const std::vector<Vec2f>& polygon);

}  // namespace common