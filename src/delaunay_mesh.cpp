#include "delaunay_mesh.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace core {

// ------------------------------------------------------------------------------------------------
bool IsValid(VertexIndex i) { return i.i != kInvalidIndex; }

// ------------------------------------------------------------------------------------------------
bool IsValid(QuarterEdgeIndex i) { return i.i != kInvalidIndex; }

// ------------------------------------------------------------------------------------------------
bool IsDualEdge(const QuarterEdge& qe) { return !IsValid(qe.i_vertex); }

// ------------------------------------------------------------------------------------------------
bool IsPrimalEdge(const QuarterEdge& qe) { return IsValid(qe.i_vertex); }

// ------------------------------------------------------------------------------------------------
DelaunayMesh::DelaunayMesh(float bounding_radius, float min_dist_to_vertex,
                           float min_dist_to_edge) :
    bounding_radius_(bounding_radius),
    min_dist_to_vertex_(min_dist_to_vertex),
    min_dist_to_edge_(min_dist_to_edge) {
    // The triangle radius is 2r + eps(), which guarantees that it is large enough.
    float r = 2 * bounding_radius + min_dist_to_edge + min_dist_to_vertex;

    VertexIndex a = AddVertex(r * std::cos(90 * M_PI / 180.0), r * std::sin(90 * M_PI / 180.0));
    VertexIndex b = AddVertex(r * std::cos(210 * M_PI / 180.0), r * std::sin(210 * M_PI / 180.0));
    VertexIndex c = AddVertex(r * std::cos(-30 * M_PI / 180.0), r * std::sin(-30 * M_PI / 180.0));

    QuarterEdgeIndex ab = AddEdge(a, b);
    QuarterEdgeIndex bc = AddEdge(b, c);
    QuarterEdgeIndex ca = AddEdge(c, a);

    Splice(Sym(ab), bc);
    Splice(Sym(bc), ca);
    Splice(Sym(ca), ab);
}

// ------------------------------------------------------------------------------------------------
DelaunayMesh::~DelaunayMesh() {}

// ------------------------------------------------------------------------------------------------
void DelaunayMesh::Clear() {
    // Move all vertices to the free list
    while (IsValid(i_vertex_alive_last_)) {
        FreeVertex(i_vertex_alive_last_);
    }
    n_vertices_ = 0;  // should not be necessary

    // Move all quarter edges to the free list
    while (IsValid(i_qe_alive_last_)) {
        FreeQuarterEdge(i_qe_alive_last_);
    }
    n_quarter_edges_ = 0;  // should not be necessary
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMesh::LoadFromData(const u8* data) {
    // Reset the mesh.
    // Note that we do this rather than call Clear() because
    // we need the indices in the save file to be the same as the indices in the new mesh.
    // If we simply call Clear(), then our free/alive list may have a weird ordering.

    vertices_.clear();
    i_vertex_alive_first_ = {kInvalidIndex};
    i_vertex_alive_last_ = {kInvalidIndex};
    i_vertex_free_first_ = {kInvalidIndex};
    i_vertex_free_last_ = {kInvalidIndex};
    n_vertices_ = 0;

    quarter_edges_.clear();
    i_qe_alive_first_ = {kInvalidIndex};
    i_qe_alive_last_ = {kInvalidIndex};
    i_qe_free_first_ = {kInvalidIndex};
    i_qe_free_last_ = {kInvalidIndex};
    n_quarter_edges_ = 0;

    u32 offset = 0;
    u32 n_vertices = *(u32*)(data + offset);
    offset += sizeof(u32);

    u32 n_quarter_edges = *(u32*)(data + offset);
    offset += sizeof(u32);

    // Load the vertices
    for (u32 i = 0; i < n_vertices; i++) {
        VertexData& vertex_data = Get(LivenVertex());
        vertex_data.v = *(common::Vec2f*)(data + offset);
        offset += sizeof(common::Vec2f);
    }

    // Load the quarter edges
    for (u32 i = 0; i < n_quarter_edges; i++) {
        QuarterEdge& qe = Get(LivenQuarterEdge());

        u32 vertex_index = *(u32*)(data + offset);
        offset += sizeof(u32);
        if (vertex_index != std::numeric_limits<u32>::max()) {
            qe.i_vertex = {vertex_index};
        } else {
            qe.i_vertex = {kInvalidIndex};
        }

        u32 nxt = *(u32*)(data + offset);
        offset += sizeof(u32);
        qe.i_nxt = {nxt};

        u32 rot = *(u32*)(data + offset);
        offset += sizeof(u32);
        qe.i_rot = {rot};
    }

    return true;
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::Next(QuarterEdgeIndex qe) const { return Get(qe).i_nxt; }

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::Rot(QuarterEdgeIndex qe) const { return Get(qe).i_rot; }

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::Sym(QuarterEdgeIndex qe) const { return Get(Get(qe).i_rot).i_rot; }

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::Tor(QuarterEdgeIndex qe) const {
    return Get(Get(Get(qe).i_rot).i_rot).i_rot;
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::Prev(QuarterEdgeIndex qe) const {
    return Get(Get(Get(qe).i_rot).i_nxt).i_rot;
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::Lnext(QuarterEdgeIndex qe) const {
    return Get(Get(Tor(qe)).i_nxt).i_rot;
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::GetQuarterEdgeRightHandClosestTo(QuarterEdgeIndex qe_src,
                                                                const common::Vec2f& dst) const {
    const common::Vec2f& a = GetVertex(qe_src);

    QuarterEdgeIndex qe_next = Next(qe_src);

    f32 rh = common::GetRightHandedness(a, GetVertex(Sym(qe_src)), dst);
    f32 rh_next = common::GetRightHandedness(a, GetVertex(Sym(qe_next)), dst);
    while (rh < 0.0 || rh_next >= 0.0) {
        rh = rh_next;
        qe_src = qe_next;
        qe_next = Next(qe_src);
        rh_next = common::GetRightHandedness(a, GetVertex(Sym(qe_next)), dst);
    }
    return qe_src;
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMesh::IsBoundaryVertex(VertexIndex i_vertex) const {
    return Get(i_vertex).i_self.i <= 2;  // We assume that we never free the first 3 vertices
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMesh::IsBoundaryVertex(const VertexData& vertex_data) const {
    return vertex_data.i_self.i <= 2;  // We assume that we never free the first 3 vertices
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMesh::IsBoundaryVertex(QuarterEdgeIndex qe_primal) const {
    return IsBoundaryVertex(GetVertexIndex(qe_primal));
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMesh::IsBoundaryEdge(QuarterEdgeIndex qe_primal) const {
    return IsBoundaryVertex(qe_primal) && IsBoundaryVertex(Sym(qe_primal));
}

// ------------------------------------------------------------------------------------------------
VertexIndex DelaunayMesh::LivenVertex() {
    VertexIndex i_vertex = {kInvalidIndex};
    if (IsValid(i_vertex_free_first_)) {
        // Grab the first item off the free list and use that.
        i_vertex = i_vertex_free_first_;

        VertexData& vertex_data = Get(i_vertex);
        i_vertex_free_first_ = vertex_data.i_next;
        vertex_data.i_next = {kInvalidIndex};
        vertex_data.i_prev = i_vertex_alive_last_;

    } else {
        // If there is no first item in the free list, append a new item
        i_vertex = {vertices_.size()};
        VertexData vertex_data = {};
        vertex_data.i_self = i_vertex;
        vertex_data.i_next = {kInvalidIndex};
        vertex_data.i_prev = i_vertex_alive_last_;
        vertices_.push_back(vertex_data);
    }

    // Hook up to previous item
    if (IsValid(i_vertex_alive_last_)) {
        Get(i_vertex_alive_last_).i_next = i_vertex;
    } else {
        i_vertex_alive_first_ = i_vertex;
    }
    i_vertex_alive_last_ = i_vertex;

    // Keep track of how many alive vertices we have
    n_vertices_ += 1;

    return i_vertex;
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::LivenQuarterEdge() {
    QuarterEdgeIndex i_qe = {kInvalidIndex};
    if (IsValid(i_qe_free_first_)) {
        // Grab the first item off the free list and use that.
        i_qe = i_qe_free_first_;

        QuarterEdge& qe = Get(i_qe);
        i_qe_free_first_ = qe.i_next;
        qe.i_vertex = {kInvalidIndex};
        qe.flags = 0;
        qe.i_next = {kInvalidIndex};
        qe.i_prev = i_qe_alive_last_;

    } else {
        // If there is no first item in the free list, append a new item
        i_qe = {quarter_edges_.size()};
        QuarterEdge qe = {};
        qe.i_vertex = {kInvalidIndex};
        qe.i_self = i_qe;
        qe.i_next = {kInvalidIndex};
        qe.i_prev = i_qe_alive_last_;
        quarter_edges_.push_back(qe);
    }

    // Hook up to previous item
    if (IsValid(i_qe_alive_last_)) {
        Get(i_qe_alive_last_).i_next = i_qe;
    } else {
        i_qe_alive_first_ = i_qe;
    }
    i_qe_alive_last_ = i_qe;

    // Keep track of how many alive quarter edges we have
    n_quarter_edges_ += 1;

    return i_qe;
}

// ------------------------------------------------------------------------------------------------
VertexIndex DelaunayMesh::AddVertex(float x, float y) {
    VertexIndex i_vertex = LivenVertex();
    Get(i_vertex).v = common::Vec2f(x, y);
    return i_vertex;
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::AddEdge(VertexIndex i_vertex_a, VertexIndex i_vertex_b) {
    // Add them all first, so we avoid issues with memory changing as we push new quarter edges.
    QuarterEdgeIndex i_ab = LivenQuarterEdge();
    QuarterEdgeIndex i_lr = LivenQuarterEdge();
    QuarterEdgeIndex i_ba = LivenQuarterEdge();
    QuarterEdgeIndex i_rl = LivenQuarterEdge();

    QuarterEdge& qe_ab = Get(i_ab);
    QuarterEdge& qe_lr = Get(i_lr);
    QuarterEdge& qe_ba = Get(i_ba);
    QuarterEdge& qe_rl = Get(i_rl);

    qe_ab.i_vertex = i_vertex_a;
    qe_lr.i_vertex = {kInvalidIndex};
    qe_ba.i_vertex = i_vertex_b;
    qe_rl.i_vertex = {kInvalidIndex};

    qe_ab.i_nxt = qe_ab.i_self;
    qe_ab.i_rot = qe_lr.i_self;
    qe_lr.i_nxt = qe_rl.i_self;
    qe_lr.i_rot = qe_ba.i_self;
    qe_ba.i_nxt = qe_ba.i_self;
    qe_ba.i_rot = qe_rl.i_self;
    qe_rl.i_nxt = qe_lr.i_self;
    qe_rl.i_rot = qe_ab.i_self;

    return qe_ab.i_self;
}

// ------------------------------------------------------------------------------------------------
void DelaunayMesh::FreeVertex(VertexIndex i) {
    VertexData& data = Get(i);

    // Hook the previous and next alive vertices
    if (IsValid(data.i_prev)) {
        Get(data.i_prev).i_next = data.i_next;
    } else {
        // This is the first alive item
        i_vertex_alive_first_ = data.i_next;
    }
    if (IsValid(data.i_next)) {
        Get(data.i_next).i_prev = data.i_prev;
    } else {
        // This is the last alive item
        i_vertex_alive_last_ = data.i_prev;
    }

    // Now append the vertex into the free list
    data.i_prev = i_vertex_free_last_;
    i_vertex_free_last_ = data.i_self;
    if (IsValid(data.i_prev)) {
        Get(data.i_prev).i_next = data.i_self;
    } else {
        // This is the first free item
        i_vertex_free_first_ = data.i_self;
    }
}

// ------------------------------------------------------------------------------------------------
void DelaunayMesh::FreeQuarterEdge(QuarterEdgeIndex i) {
    QuarterEdge& qe = Get(i);

    // Hook the previous and next alive vertices
    if (IsValid(qe.i_prev)) {
        Get(qe.i_prev).i_next = qe.i_next;
    } else {
        // This is the first alive item
        i_qe_alive_first_ = qe.i_next;
    }
    if (IsValid(qe.i_next)) {
        Get(qe.i_next).i_prev = qe.i_prev;
    } else {
        // This is the last alive item
        i_qe_alive_last_ = qe.i_prev;
    }

    // Now append the vertex into the free list
    qe.i_prev = i_qe_free_last_;
    i_qe_free_last_ = qe.i_self;
    if (IsValid(qe.i_prev)) {
        Get(qe.i_prev).i_next = qe.i_self;
    } else {
        // This is the first free item
        i_qe_free_first_ = qe.i_self;
    }
}

// ------------------------------------------------------------------------------------------------
void DelaunayMesh::SwapNexts(QuarterEdgeIndex a, QuarterEdgeIndex b) {
    QuarterEdgeIndex tmp = Get(a).i_nxt;
    Get(a).i_nxt = Get(b).i_nxt;
    Get(b).i_nxt = tmp;
}

// ------------------------------------------------------------------------------------------------
void DelaunayMesh::Splice(QuarterEdgeIndex a, QuarterEdgeIndex b) {
    SwapNexts(Get(Get(a).i_nxt).i_rot, Get(Get(b).i_nxt).i_rot);
    SwapNexts(a, b);
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::GetEnclosingTriangle(const common::Vec2f& p,
                                                    QuarterEdgeIndex qe_dual) const {
    // This should always terminate
    constexpr int kMaxIters = 1000;
    for (int iter = 0; iter < kMaxIters; iter++) {
        const auto [qe_ab, qe_bc, qe_ca] = GetTriangleQuarterEdges(qe_dual);

        const common::Vec2f& a = GetVertex(qe_ab);
        const common::Vec2f& b = GetVertex(qe_bc);
        const common::Vec2f& c = GetVertex(qe_ca);

        if (common::GetRightHandedness(a, b, p) < 0) {
            qe_dual = Rot(qe_ab);  // Move across AB
        } else if (common::GetRightHandedness(b, c, p) < 0) {
            qe_dual = Rot(qe_bc);  // Move across BC
        } else if (common::GetRightHandedness(c, a, p) < 0) {
            qe_dual = Rot(qe_ca);  // Move across CA
        } else {
            return qe_dual;
        }
    }
    return qe_dual;
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::GetEnclosingTriangle(const common::Vec2f& p) const {
    return GetEnclosingTriangle(p, Tor(quarter_edges_.front().i_self));
}

// ------------------------------------------------------------------------------------------------
std::tuple<QuarterEdgeIndex, QuarterEdgeIndex, QuarterEdgeIndex>
DelaunayMesh::GetTriangleQuarterEdges(QuarterEdgeIndex qe_dual) const {
    QuarterEdgeIndex qe_nxt1 = Get(qe_dual).i_nxt;
    QuarterEdgeIndex qe_nxt2 = Get(qe_nxt1).i_nxt;
    QuarterEdgeIndex qe_ab = Get(qe_dual).i_rot;
    QuarterEdgeIndex qe_bc = Get(qe_nxt1).i_rot;
    QuarterEdgeIndex qe_ca = Get(qe_nxt2).i_rot;
    return std::make_tuple(qe_ab, qe_bc, qe_ca);
}

// ------------------------------------------------------------------------------------------------
void DelaunayMesh::MoveVertexToward(QuarterEdgeIndex qe_primal, const common::Vec2f& pos) {
    const f32 eps = 1e-3;  // TODO - better place to define this?

    // Where the point is now
    const common::Vec2f& a = GetVertex(qe_primal);
    const common::Vec2f dir = pos - a;
    if (common::Norm(dir) < eps) {
        // We do this check to avoid issues with coincident points later on.
        // No sense doing work if we are already at pos.
        return;  // no change
    }

    // Figure out which face we are moving in by rotationg around qe_primal
    // and finding the first pair where the cone contains pos.

    // The quarter edge across from a
    QuarterEdgeIndex qe_b = Sym(qe_primal);
    common::Vec2f b = GetVertex(qe_b);

    // The quarter edge to the right of a->b
    QuarterEdgeIndex qe_c = Sym(Prev(qe_b));
    common::Vec2f c = GetVertex(qe_c);

    while (common::GetRightHandedness(a, b, pos) < 0 || common::GetRightHandedness(a, pos, c) < 0) {
        qe_primal = Next(qe_primal);  // Rotate to the next primal edge
        qe_b = Sym(qe_primal);
        b = GetVertex(qe_b);
        qe_c = Sym(Prev(qe_b));
        c = GetVertex(qe_c);
    }

    // The polygon enclosing vertex A may not be convex. We do not want to allow shifting A to
    // anywhere that would cause edges to shift their order. As such, walk around the polygon and
    // intersect each edge with both AB and AC, clipping them as appropriate.
    QuarterEdgeIndex qe_start = qe_primal;
    QuarterEdgeIndex qe = Next(qe_primal);
    while (qe != qe_start) {
        QuarterEdgeIndex qe_d = Sym(qe);
        const common::Vec2f& d = GetVertex(qe_d);
        QuarterEdgeIndex qe_e = Sym(Prev(qe_d));
        const common::Vec2f& e = GetVertex(qe_e);

        common::LineIntersectionResult result = common::CalcLineIntersection(a, b, d, e);
        if (result.intersection == common::Intersection::AT_POINT && result.s > 0.0 &&
            result.s < 1.0) {
            // Shorten.
            b = a + (b - a) * result.s;
        }

        result = common::CalcLineIntersection(a, c, d, e);
        if (result.intersection == common::Intersection::AT_POINT && result.s > 0.0 &&
            result.s < 1.0) {
            // Shorten.
            c = a + (c - a) * result.s;
        }

        // Advance
        qe = Next(qe);
    }

    // We now have a, b, c, where a is the origin and pos lies in the cone between a->b, a->c.
    // We cannot move a past the edge b->c (and have to stay a bit further back).
    // So shift b->c to the left by min_dist_to_edge.
    common::Vec2f shift = common::Rotr(common::Normalize(c - b)) * min_dist_to_edge_;
    b += shift;
    c += shift;

    common::LineIntersectionResult result = common::CalcLineIntersection(a, pos, b, c);
    if (result.intersection != common::Intersection::AT_POINT) {
        return;  // nothing we can do.
    }

    // TODO: enforce min_dist_to_vertex if it is bigger than min_dist_to_edge.

    // Pos's current interpolant is 1.0; clamp it below the s-interpolant.
    // If the successor triangle is not right-handed, we do accept movement all the way to
    // the boundary.
    if (result.s >= 1.0 && common::GetRightHandedness(pos, b, c) > 0) {
        // Accept pos
        Get(Get(qe_primal).i_vertex).v = pos;
    } else {
        // Accept the maximum movement
        Get(Get(qe_primal).i_vertex).v = a + dir * result.s;
    }
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMesh::MaybeFlipEdge(QuarterEdgeIndex qe_primal) {
    // We cannot flip the edge if it is constrained.
    if (IsConstrained(qe_primal)) {
        return false;
    }

    // Do not flip the boundary edge
    VertexData& a_data = Get(Get(qe_primal).i_vertex);
    VertexData& b_data = Get(Get(Sym(qe_primal)).i_vertex);
    int num_boundary_vertices = IsBoundaryVertex(a_data) + IsBoundaryVertex(b_data);
    if (num_boundary_vertices == 2) {  // one edge being on the boundary is okay
        return false;
    }

    const common::Vec2f& a = a_data.v;
    const common::Vec2f& b = b_data.v;
    const common::Vec2f& c = GetVertex(Sym(Prev(qe_primal)));  // right of A->B
    const common::Vec2f& d = GetVertex(Sym(Next(qe_primal)));  // left of A->B

    // We can only flip the edge if its surrounding quad is convex.
    if (common::GetRightHandedness(c, b, d) <= 0 || common::GetRightHandedness(c, d, a) <= 0) {
        return false;
    }

    // Flip!
    FlipEdgeImpl(qe_primal);

    return true;
}

// ------------------------------------------------------------------------------------------------
void DelaunayMesh::FlipEdgeImpl(QuarterEdgeIndex qe) {
    QuarterEdgeIndex qe_sym = Sym(qe);
    QuarterEdgeIndex qe_a = Prev(qe);
    QuarterEdgeIndex qe_b = Prev(qe_sym);

    Splice(qe, qe_a);
    Splice(qe_sym, qe_b);
    Splice(qe, Lnext(qe_a));
    Splice(qe_sym, Lnext(qe_b));

    Get(qe).i_vertex = Get(Sym(qe_a)).i_vertex;
    Get(qe_sym).i_vertex = Get(Sym(qe_b)).i_vertex;
}

// ------------------------------------------------------------------------------------------------
void DelaunayMesh::EnforceLocallyDelaunay(QuarterEdgeIndex qe_start) {
    const common::Vec2f& p = GetVertex(qe_start);

    // Walk around the outer edges and flip any edges that are not locally delaunay.
    // We can check an edge by seeing if the opposite vertex is within the inscribed
    // circle of the inner vertex + edge vertices.

    // We start at the given quarter edge, then walk around the polygon surrounding the source
    // vertex until we have checked all edges and return to the start.

    QuarterEdgeIndex qe_index = qe_start;
    bool done = false;
    while (!done) {
        // Get the quarter edge representing an edge along the perimeter polygon.
        QuarterEdgeIndex qe_outer_edge = Get(Sym(qe_index)).i_nxt;

        // Advance
        QuarterEdgeIndex qe_index_start_of_iter = qe_index;
        qe_index = Next(qe_index);
        done = qe_index == qe_start;

        // Do not flip the boundary edge
        VertexData& src_data = Get(Get(qe_outer_edge).i_vertex);
        VertexData& dst_data = Get(Get(Sym(qe_outer_edge)).i_vertex);

        int num_boundary_vertices = IsBoundaryVertex(src_data) + IsBoundaryVertex(dst_data);
        if (num_boundary_vertices == 2) {  // one edge being on the boundary is okay
            continue;
        }

        // Do not flip a constrained edge.
        if (IsConstrained(qe_outer_edge)) {
            continue;
        }

        // Check the edge from qe_outer_edge to sym(qe_outer_edge)
        const common::Vec2f& src = src_data.v;
        const common::Vec2f& dst = dst_data.v;

        // Get the far vertex, across the dividing edge, on the other side of P.
        VertexData& far_data = Get(Get(Sym(Get(qe_outer_edge).i_nxt)).i_vertex);

        const common::Vec2f& far = far_data.v;

        // If the edge contains a boundary vertex, don't flip it if it would produce an  inside -
        // out triangle.
        if (num_boundary_vertices == 1) {
            if (common::GetTriangleContainment(dst, far, p, src) >= 0 ||
                common::GetTriangleContainment(dst, far, src, p) >= 0 ||
                common::GetTriangleContainment(src, far, dst, p) >= 0 ||
                common::GetTriangleContainment(src, far, p, dst) >= 0) {
                continue;
            }
        }

        // If we are not locally Delaunay, flip the edge
        if (common::GetCircleContainment(p, src, dst, far) > 0 ||
            common::GetCircleContainment(far, p, dst, src) > 0) {
            // Either p is inside the circle passing through src, dst, and far, or
            //  far is inside the circle passing through p, src, and dst.
            // We have to flip the edge.
            FlipEdgeImpl(qe_outer_edge);

            // We flipped the edge, so qe_outer_edge has to be traversed again.
            // Back it up.
            qe_index = qe_index_start_of_iter;
            done = false;
        }
    }
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::GetQuarterEdge(VertexIndex a) const {
    for (const QuarterEdge& qe : quarter_edges_) {
        if (qe.i_vertex == a) {
            return qe.i_self;
        }
    }

    return QuarterEdgeIndex({kInvalidIndex});
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::GetQuarterEdge(VertexIndex a, VertexIndex b) const {
    // If we already have an edge between these two vertices, we are done.
    for (const QuarterEdge& qe : quarter_edges_) {
        if ((qe.i_vertex == a) && (Get(Sym(qe.i_self)).i_vertex == b)) {
            return qe.i_self;
        }
    }

    return QuarterEdgeIndex({kInvalidIndex});
}

// ------------------------------------------------------------------------------------------------
VertexIndex DelaunayMesh::GetNext(VertexIndex i) const {
    if (IsValid(i)) {
        return Get(i).i_next;
    }
    return i;
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::GetNext(QuarterEdgeIndex i) const {
    if (IsValid(i)) {
        return Get(i).i_next;
    }
    return i;
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMesh::IsDual(QuarterEdgeIndex i) const {
    const QuarterEdge& qe = GetQuarterEdge(i);
    return IsDualEdge(qe);
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMesh::IsPrimal(QuarterEdgeIndex i) const {
    const QuarterEdge& qe = GetQuarterEdge(i);
    return IsPrimalEdge(qe);
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMesh::IsConstrained(QuarterEdgeIndex i) const {
    const QuarterEdge& qe = GetQuarterEdge(i);
    return (qe.flags & QE_FLAG_CONSTRAINED) > 0;
}

// ------------------------------------------------------------------------------------------------
void DelaunayMesh::ConstrainEdge(QuarterEdgeIndex qe) {
    Get(qe).flags |= QE_FLAG_CONSTRAINED;
    for (int i = 0; i < 3; i++) {
        qe = Rot(qe);
        Get(qe).flags |= QE_FLAG_CONSTRAINED;
    }
}

// ------------------------------------------------------------------------------------------------
void DelaunayMesh::UnconstrainEdge(QuarterEdgeIndex qe) {
    Get(qe).flags &= ~QE_FLAG_CONSTRAINED;
    for (int i = 0; i < 3; i++) {
        qe = Rot(qe);
        Get(qe).flags &= ~QE_FLAG_CONSTRAINED;
    }
}

// ------------------------------------------------------------------------------------------------
// bool DelaunayMesh::HasEdge(int i, int j) const {
//     if (i == j) {
//         return false;  // Self edges do not exist
//     }

//     int n = NumVertices();
//     if (i < 0 || i >= n || j < 0 || j >= n) {
//         return false;  // Invalid index
//     }

//     // If we already have an edge between these two vertices, we are done.
//     VertexData* a = vertices_.at(i);
//     VertexData* b = vertices_.at(j);
//     for (QuarterEdge* qe : quarter_edges_) {
//         if (qe->vertex == a) {
//             if (Sym(qe)->vertex == b) {
//                 return true;
//             }
//         }
//     }

//     return false;
// }

// ------------------------------------------------------------------------------------------------
DelaunayMesh::InsertVertexResult DelaunayMesh::InsertVertex(const common::Vec2f& p) {
    InsertVertexResult result = {};
    result.i_vertex = {kInvalidIndex};
    result.i_qe = {kInvalidIndex};

    // Ensure that the point is within the valid bounds
    if (common::Norm(p) >= bounding_radius_) {
        result.category = InsertVertexResultCategory::OUT_OF_BOUNDS;
        return result;
    }

    // Get the enclosing triangle
    QuarterEdgeIndex qe_dual = GetEnclosingTriangle(p);
    auto [qe_ab, qe_bc, qe_ca] = GetTriangleQuarterEdges(qe_dual);

    // Grab the vertices
    common::Vec2f a = GetVertex(qe_ab);
    common::Vec2f b = GetVertex(qe_bc);
    common::Vec2f c = GetVertex(qe_ca);

    // If we are too close to an existing vertex, do nothing.
    f32 dist_ap = common::Norm(a - p);
    f32 dist_bp = common::Norm(b - p);
    f32 dist_cp = common::Norm(c - p);

    result.i_vertex = Get(qe_ab).i_vertex;
    f32 min_dist = dist_ap;
    if (dist_bp < min_dist) {
        min_dist = dist_bp;
        result.i_vertex = Get(qe_bc).i_vertex;
    }
    if (dist_cp < min_dist) {
        min_dist = dist_cp;
        result.i_vertex = Get(qe_ca).i_vertex;
    }

    if (min_dist < min_dist_to_vertex_) {
        // We are coincident with an existing vertex.
        result.category = InsertVertexResultCategory::COINCIDENT;
        // Returns the index to that vertex.
        return result;
    }

    // Add the vertex
    result.i_vertex = AddVertex(p.x, p.y);

    // Check whether we are (effectively) on an edge or inside a face.
    float dist_to_ab = common::GetDistanceToLine(p, a, b);
    float dist_to_bc = common::GetDistanceToLine(p, b, c);
    float dist_to_ca = common::GetDistanceToLine(p, c, a);
    if (std::min({dist_to_ab, dist_to_bc, dist_to_ca}) > min_dist_to_edge_) {
        // Normal case. We are not too close to an edge.
        result.category = InsertVertexResultCategory::IN_FACE;

        // Add the new edges and splice them in
        QuarterEdgeIndex ap = AddEdge(Get(qe_ab).i_vertex, result.i_vertex);
        QuarterEdgeIndex bp = AddEdge(Get(qe_bc).i_vertex, result.i_vertex);
        QuarterEdgeIndex cp = AddEdge(Get(qe_ca).i_vertex, result.i_vertex);

        Splice(ap, qe_ab);
        Splice(bp, qe_bc);
        Splice(cp, qe_ca);

        // TODO: Figure out splice such that we get this desired effect.
        //       i.e. can we replace these next calls with three calls to splice?
        QuarterEdgeIndex pa = Sym(ap);
        QuarterEdgeIndex pb = Sym(bp);
        QuarterEdgeIndex pc = Sym(cp);

        Get(pa).i_nxt = pb;
        Get(pb).i_nxt = pc;
        Get(pc).i_nxt = pa;

        Get(Get(pa).i_rot).i_nxt = Get(cp).i_rot;
        Get(Get(pb).i_rot).i_nxt = Get(ap).i_rot;
        Get(Get(pc).i_rot).i_nxt = Get(bp).i_rot;

        // Give back a quarter edge with p as its source.
        result.i_qe = pa;
    } else {
        // We are effectively on an edge.
        // Identify that edge, and cut it with the new vertex.
        result.category = InsertVertexResultCategory::ON_EDGE;

        // Let DE be the edge we are on, F be the far vertex qe, and G be the vertex qe across
        // DE from F.
        QuarterEdgeIndex qe_dp;
        QuarterEdgeIndex qe_ef;
        QuarterEdgeIndex qe_fd;
        if (dist_to_ab <= dist_to_bc && dist_to_ab <= dist_to_ca) {
            qe_dp = qe_ab;
            qe_ef = qe_bc;
            qe_fd = qe_ca;
        } else if (dist_to_bc <= dist_to_ab && dist_to_bc <= dist_to_ca) {
            qe_dp = qe_bc;
            qe_ef = qe_ca;
            qe_fd = qe_ab;
        } else {
            qe_dp = qe_ca;
            qe_ef = qe_ab;
            qe_fd = qe_bc;
        }

        QuarterEdgeIndex qe_ge = Prev(Sym(Prev(qe_dp)));

        // Our edge may not be the boundary edge
        int n_boundary_vertices = IsBoundaryVertex(Get(Get(qe_dp).i_vertex)) +
                                  IsBoundaryVertex(Get(Get(Sym(qe_dp)).i_vertex));
        if (n_boundary_vertices == 2) {
            result.category = InsertVertexResultCategory::OUT_OF_BOUNDS;
            return result;
        }

        // Grab another quarter edge we will need
        QuarterEdgeIndex qe_dg = Prev(qe_dp);

        // Reset the things that point to DE
        Get(qe_dg).i_nxt = Sym(qe_fd);
        QuarterEdgeIndex qe_ge_tor = Tor(qe_ge);
        Get(qe_ge_tor).i_nxt = Rot(Sym(qe_ef));
        Get(qe_ef).i_nxt = Sym(qe_ge);
        QuarterEdgeIndex qe_fd_tor = Tor(qe_fd);
        Get(qe_fd_tor).i_nxt = Rot(Sym(qe_dg));

        // Reset DP as a quarter edge, and change it to PD
        QuarterEdgeIndex qe_pd = Sym(qe_dp);
        Get(qe_pd).i_vertex = result.i_vertex;
        Get(qe_pd).i_nxt = qe_pd;
        Get(Get(qe_pd).i_rot).i_nxt = Rot(qe_dp);
        Get(qe_dp).i_nxt = qe_dp;
        Get(Get(qe_dp).i_rot).i_nxt = Rot(qe_pd);

        // Create three new edges EP, FP, and GP
        QuarterEdgeIndex qe_ep = AddEdge(Get(qe_ef).i_vertex, result.i_vertex);
        QuarterEdgeIndex qe_fp = AddEdge(Get(qe_fd).i_vertex, result.i_vertex);
        QuarterEdgeIndex qe_gp = AddEdge(Get(qe_ge).i_vertex, result.i_vertex);

        // Splice them all
        Splice(qe_dp, qe_dg);
        Splice(qe_gp, qe_ge);
        Splice(qe_ep, qe_ef);
        Splice(qe_fp, qe_fd);

        // TODO: Figure out splice such that we get this desired effect.
        //       i.e. can we replace these next calls with three calls to splice?
        qe_pd = Sym(qe_dp);
        QuarterEdgeIndex qe_pe = Sym(qe_ep);
        QuarterEdgeIndex qe_pf = Sym(qe_fp);
        QuarterEdgeIndex qe_pg = Sym(qe_gp);

        Get(qe_pd).i_nxt = Sym(qe_gp);
        Get(qe_pe).i_nxt = Sym(qe_fp);
        Get(qe_pf).i_nxt = Sym(qe_dp);
        Get(qe_pg).i_nxt = Sym(qe_ep);

        Get(Rot(qe_pd)).i_nxt = Rot(qe_fp);
        Get(Rot(qe_pe)).i_nxt = Rot(qe_gp);
        Get(Rot(qe_pf)).i_nxt = Rot(qe_ep);
        Get(Rot(qe_pg)).i_nxt = Rot(qe_dp);

        // Give back a quarter edge with p as its source that points along the split edge.
        result.i_qe = qe_pd;
    }

    return result;
}

// ------------------------------------------------------------------------------------------------
DelaunayMesh::EnforceEdgeInternalResult DelaunayMesh::EnforceEdgeInternal(QuarterEdgeIndex qe_a,
                                                                          VertexIndex i_vertex_b) {
    EnforceEdgeInternalResult result;
    result.qe_ab = {kInvalidIndex};
    result.progress = false;

    // Find either an existing quarter edge from A to B, or the quarter edge from A to C that is
    // immediately counter-clockwise of A->B.

    const common::Vec2f& b = GetVertex(i_vertex_b);
    qe_a = GetQuarterEdgeRightHandClosestTo(qe_a, b);

    while (Get(Sym(qe_a)).i_vertex != i_vertex_b) {
        // The edge is not pointing to b. We are point to C instead, and are immediately CCW of B.
        // We need to flip CD, where D is on the other side of AB from C.
        QuarterEdgeIndex qe_cd = Prev(Sym(qe_a));

        bool flipped = MaybeFlipEdge(qe_cd);
        if (!flipped) {
            return result;  // FAILED!
        }

        // Flipped an edge
        result.progress = true;

        // After flipping, we would have an edge pointing from A to E, where E completes the quad
        // on the other side of CD from A. It is possible that E = B, or that E lies on either side
        // of AB.
        // Effectively start over.
        qe_a = GetQuarterEdgeRightHandClosestTo(qe_a, b);
    }

    result.qe_ab = qe_a;
    return result;
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::EnforceEdge(QuarterEdgeIndex qe_a, QuarterEdgeIndex qe_b) {
    // TODO: If we have vertices intersecting the edge, we should split our constrained edge
    //       by those vertices and call this for each sub-edge.

    // TODO: If we try to constrain an edge that overlaps with an already-constrained edge,
    //       then we need to produce an intersection and fix it that way.

    if (!IsValid(qe_a) || !IsValid(qe_b)) {
        return {kInvalidIndex};
    }

    VertexIndex i_vertex_a = GetVertexIndex(qe_a);
    VertexIndex i_vertex_b = GetVertexIndex(qe_b);
    if (i_vertex_a == i_vertex_b) {
        return {kInvalidIndex};  // Cannot add self-edges
    }

    // Tackle this problem from both directions until solved or no progress is made.
    bool progress = true;
    bool success = false;
    QuarterEdgeIndex qe_ab = {kInvalidIndex};
    while (progress && !success) {
        progress = false;

        EnforceEdgeInternalResult result = EnforceEdgeInternal(qe_a, i_vertex_b);
        progress |= result.progress;
        qe_ab = result.qe_ab;
        success |= IsValid(qe_ab);

        if (success) {
            break;
        }

        result = EnforceEdgeInternal(qe_b, i_vertex_a);
        progress |= result.progress;
        qe_ab = Sym(result.qe_ab);
        success |= IsValid(qe_ab);
    }

    return qe_ab;
}

}  // namespace core