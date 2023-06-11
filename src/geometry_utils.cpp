#include "geometry_utils.hpp"

// #include <clipper2/clipper.h>

namespace common {

// ------------------------------------------------------------------------------------------------
Vec2i ToFixedPoint(const Vec2f& v, float resolution) {
    long x = static_cast<long>(v.x * resolution);
    long y = static_cast<long>(v.y * resolution);
    return Vec2i({x, y});
}

// ------------------------------------------------------------------------------------------------
Vec2f ToFloatingPoint(const Vec2i& v, float resolution) {
    float x = v.x / resolution;
    float y = v.y / resolution;
    return Vec2f({x, y});
}

// ------------------------------------------------------------------------------------------------
std::vector<Vec2i> ToFixedPoint(const std::vector<Vec2f>& vertices, float resolution) {
    // Reserve vertices
    std::vector<Vec2i> fp_vertices;
    fp_vertices.reserve(vertices.size());

    // Copy over the vertices
    for (const Vec2f& v : vertices) {
        fp_vertices.push_back(ToFixedPoint(v, resolution));
    }

    return fp_vertices;
}

// ------------------------------------------------------------------------------------------------
std::vector<Vec2f> ToFloatingPoint(const std::vector<Vec2i>& vertices, float resolution) {
    // Reserve vertices
    std::vector<Vec2f> fp_vertices;
    fp_vertices.reserve(vertices.size());

    // Copy over the vertices
    for (const Vec2i& v : vertices) {
        fp_vertices.push_back(ToFloatingPoint(v, resolution));
    }

    return fp_vertices;
}

// ------------------------------------------------------------------------------------------------
float Norm(const Vec2f& v) { return std::hypot(v.x, v.y); }
float Norm(const Vec2i& v) { return std::hypot(v.x, v.y); }
Vec2f Normalize(const Vec2f& v) { return v / Norm(v); }

// ------------------------------------------------------------------------------------------------
float SqNorm(const Vec2f& v) { return v.x * v.x + v.y * v.y; }
long long SqNorm(const Vec2i& v) { return v.x * v.x + v.y * v.y; }

// ------------------------------------------------------------------------------------------------
float Dist(const Vec2f& a, const Vec2f& b) { return Norm(a - b); }
float Dist(const Vec2i& a, const Vec2i& b) { return Norm(a - b); }

// ------------------------------------------------------------------------------------------------
float Dot(const Vec2f& a, const Vec2f& b) { return a.x * b.x + a.y * b.y; }

// ------------------------------------------------------------------------------------------------
float Cross(const Vec2f& a, const Vec2f& b) { return a.x * b.y - a.y * b.x; }
long long Cross(const Vec2i& a, const Vec2i& b) {
    return (long long)a.x * b.y - (long long)a.y * b.x;
};

// ------------------------------------------------------------------------------------------------
Vec2f VectorProjection(const Vec2f& a, const Vec2f& b) { return b * (Dot(a, b) / Dot(b, b)); }

// ------------------------------------------------------------------------------------------------
Vec2f ClampNorm(const Vec2f& v, float threshold) {
    float norm = Norm(v);
    if (norm > threshold) {
        return v * (threshold / norm);
    }
    return v;
}

// ------------------------------------------------------------------------------------------------
Vec2f Interp(const Vec2f& a, const Vec2f& b, float interpolant) {
    return a + (b - a) * interpolant;
}

// ------------------------------------------------------------------------------------------------
Vec2i Interp(const Vec2i& a, const Vec2i& b, float interpolant) {
    float x = a.x + (b.x - a.x) * interpolant;
    float y = a.y + (b.y - a.y) * interpolant;
    return Vec2i({(int)round(x), (int)round(y)});
}

// ------------------------------------------------------------------------------------------------
Vec2f Rotr(const Vec2f& a) { return Vec2f(-a.y, a.x); }

// ------------------------------------------------------------------------------------------------
float GetRightHandedness(const Vec2f& a, const Vec2f& b, const Vec2f& c) {
    return Det(a.x, a.y, 1.0f, b.x, b.y, 1.0f, c.x, c.y, 1.0f);
}

// ------------------------------------------------------------------------------------------------
float GetDistanceToLine(const Vec2f& p, const Vec2f& a, const Vec2f& b) {
    Vec2f delta = b - a;
    return std::abs(delta.x * (a.y - p.y) - delta.y * (a.x - p.x)) / Norm(delta);
}

// ------------------------------------------------------------------------------------------------
LineIntersectionResult CalcLineIntersection(const Vec2f& a, const Vec2f& b, const Vec2f& c,
                                            const Vec2f& d) {
    // Direction vectors
    Vec2f u = b - a;
    Vec2f v = d - c;
    Vec2f w = a - c;

    // Check for degeneracy (A&B or C&D are colocated)
    if (IsZero(SqNorm(u), 1e-6f) || IsZero(SqNorm(v), 1e-6f)) {
        return LineIntersectionResult({.intersection = Intersection::UNDEFINED, .s = 0, .t = 0});
    }

    // Check to see whether the lines are parallel (or nearly so)
    float denominator = Cross(u, v);
    float v_cross_w = Cross(v, w);
    float u_cross_w = Cross(u, w);
    bool is_parallel = IsZero(denominator, 1e-6f);
    if (is_parallel) {
        // If they are not colinear, then they do not intersect
        if (!IsZero(v_cross_w, 1e-6f) || !IsZero(u_cross_w, 1e-6f)) {
            return LineIntersectionResult({.intersection = Intersection::NONE, .s = 0, .t = 0});
        }

        // They overlap.
        return LineIntersectionResult({.intersection = Intersection::OVERLAP, .s = 0, .t = 0});
    }

    // Not parallel.
    float s = v_cross_w / denominator;
    float t = u_cross_w / denominator;
    return LineIntersectionResult({.intersection = Intersection::AT_POINT, .s = s, .t = t});
}

// ------------------------------------------------------------------------------------------------
LineIntersectionResult CalcLineIntersection(const Vec2i& a, const Vec2i& b, const Vec2i& c,
                                            const Vec2i& d) {
    // Direction vectors
    Vec2i u = b - a;
    Vec2i v = d - c;
    Vec2i w = a - c;

    // Check for degeneracy (A&B or C&D are colocated)
    if (SqNorm(u) == 0 || SqNorm(v) == 0) {
        return LineIntersectionResult({.intersection = Intersection::UNDEFINED, .s = 0, .t = 0});
    }

    // Check to see whether the lines are parallel (or nearly so)
    long long denominator = Cross(u, v);
    long long v_cross_w = Cross(v, w);
    long long u_cross_w = Cross(u, w);
    bool is_parallel = (denominator == 0);
    if (is_parallel) {
        // If they are not colinear, then they do not intersect
        if (v_cross_w != 0 || u_cross_w != 0) {
            return LineIntersectionResult({.intersection = Intersection::NONE, .s = 0, .t = 0});
        }

        // They overlap.
        return LineIntersectionResult({.intersection = Intersection::OVERLAP, .s = 0, .t = 0});
    }

    // Not parallel.
    float s = static_cast<float>(v_cross_w) / denominator;
    float t = static_cast<float>(u_cross_w) / denominator;
    return LineIntersectionResult({.intersection = Intersection::AT_POINT, .s = s, .t = t});
}

// ------------------------------------------------------------------------------------------------
int GetLineSegmentContainment(float s) {
    if (s > 0 && s < 1) {
        return 1;
    } else if (s >= 0 && s <= 1) {
        return 0;
    } else {
        return -1;
    }
}

// ------------------------------------------------------------------------------------------------
int GetLineSegmentContainment(float s, float t) {
    int containment_ab = GetLineSegmentContainment(s);
    int containment_cd = GetLineSegmentContainment(t);
    return std::min(containment_ab, containment_cd);
}

// ------------------------------------------------------------------------------------------------
int GetTriangleContainment(const Vec2f& p, const Vec2f& a, const Vec2f& b, const Vec2f& c) {
    float ab = GetRightHandedness(a, b, p);
    float bc = GetRightHandedness(b, c, p);
    float ca = GetRightHandedness(c, a, p);
    if (ab > 0 && bc > 0 && ca > 0) {
        return 1;
    } else if (ab >= 0 && bc >= 0 && ca >= 0) {
        return 0;
    } else {
        return -1;
    }
}

// ------------------------------------------------------------------------------------------------
int GetCircleContainment(const Vec2f& p, const Vec2f& a, const Vec2f& b, const Vec2f& c) {
    // clang-format off
    float d = Det(a.x, a.y, a.x*a.x + a.y*a.y, 1.0f,
                  b.x, b.y, b.x*b.x + b.y*b.y, 1.0f,
                  c.x, c.y, c.x*c.x + c.y*c.y, 1.0f,
                  p.x, p.y, p.x*p.x + p.y*p.y, 1.0f);
    // clang-format on
    if (d > 0) {
        return 1;
    } else if (d < 0) {
        return -1;
    } else {
        return 0;
    }
}

// ------------------------------------------------------------------------------------------------
size_t GetLowestVertexIndex(const std::vector<Vec2i>& polygon) {
    // TODO: This could be found with a version of bisection search if the polygon is convex.
    //       That probably is overkill for us, since our polygons are usually small.
    size_t best_index = 0;
    for (size_t i = 1; i < polygon.size(); i++) {
        if (polygon[i] < polygon[best_index]) {
            best_index = i;
        }
    }
    return best_index;
}

// ------------------------------------------------------------------------------------------------
std::vector<Vec2i> MinkowskiSum(const std::vector<Vec2i>& P, const std::vector<Vec2i>& Q) {
    // We assume that both P and Q have vertices in CCW order.

    size_t i_lowest = GetLowestVertexIndex(P);
    size_t j_lowest = GetLowestVertexIndex(Q);

    std::vector<Vec2i> retval;

    size_t i = 0, j = 0;
    while (i < P.size() || j < Q.size()) {
        // Get p and q via cyclic coordinates.
        const Vec2i& p = GetCyclicVertex(P, i_lowest + i);
        const Vec2i& q = GetCyclicVertex(Q, j_lowest + j);
        const Vec2i& p_next = GetCyclicVertex(P, i_lowest + i + 1);
        const Vec2i& q_next = GetCyclicVertex(Q, j_lowest + j + 1);

        // Append p + q
        retval.push_back(p + q);

        // Use the cross product to gage polar angle
        long long cross = Cross(p_next - p, q_next - q);

        if (cross >= 0 && i < P.size()) {
            i++;
        }
        if (cross <= 0 && j < Q.size()) {
            j++;
        }
    }

    return retval;
}

// //
// ------------------------------------------------------------------------------------------------
// Clipper2Lib::Path64 MakePath64(const std::vector<Vec2i>& path) {
//     Clipper2Lib::Path64 path64;
//     path64.reserve(path.size());
//     for (const auto& p : path) {
//         path64.emplace_back(p.x, p.y);
//     }
//     return path64;
// }

// //
// ------------------------------------------------------------------------------------------------
// std::vector<Vec2i> Union(const std::vector<Vec2i>& P, const std::vector<Vec2i>& Q) {
//     std::vector<Vec2i> retval;
//     if (P.size() < 3 || Q.size() < 3) {
//         return retval;
//     }

//     Clipper2Lib::Paths64 subjects;
//     subjects.push_back(MakePath64(P));
//     subjects.push_back(MakePath64(Q));

//     // Our two polygons should be non-self-intersecting, so the fill rule should not matter.
//     Clipper2Lib::Paths64 union_result =
//         Clipper2Lib::Union(subjects, Clipper2Lib::FillRule::Positive);

//     if (union_result.size() != 1) {
//         // We did not have overlap, so no union occured.
//         return retval;
//     }

//     // Otherwise, we have one unioned polygon.
//     Clipper2Lib::Path64 trimmed = Clipper2Lib::TrimCollinear(union_result.at(0));
//     retval.reserve(trimmed.size());
//     for (const Clipper2Lib::Point64& p : trimmed) {
//         retval.emplace_back(static_cast<long>(p.x), static_cast<long>(p.y));
//     }
//     return retval;
// }

// ------------------------------------------------------------------------------------------------
Vec2f GetCentroid(const std::vector<Vec2f>& polygon) {
    Vec2f centroid = Vec2f({0.0, 0.0});
    for (const auto& v : polygon) {
        centroid += v;
    }
    return centroid / polygon.size();
}

}  // namespace common