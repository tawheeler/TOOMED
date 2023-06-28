#include "mesh.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace mesh {

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

// //
// ------------------------------------------------------------------------------------------------
// void DelaunayMesh::Clear() {
//     vertices_.clear();
//     quarter_edges_.clear();
// }

// //
// ------------------------------------------------------------------------------------------------
// bool DelaunayMesh::LoadFromData(const u8* data) {
//     // Reset the mesh
//     Clear();

//     u32 offset = 0;
//     u32 n_vertices = *(u32*)(data + offset);
//     offset += sizeof(u32);

//     u32 n_quarter_edges = *(u32*)(data + offset);
//     offset += sizeof(u32);

//     // Load the vertices
//     vertices_.reserve(n_vertices);
//     for (u32 i = 0; i < n_vertices; i++) {
//         VertexData* vertex_data = new VertexData();

//         vertex_data->index = i;
//         vertex_data->vertex = *(common::Vec2f*)(data + offset);
//         offset += sizeof(common::Vec2f);

//         vertices_.push_back(vertex_data);
//     }

//     // Load the quarter edges
//     quarter_edges_.reserve(n_quarter_edges);
//     for (u32 i = 0; i < n_quarter_edges; i++) {
//         QuarterEdge* qe = new QuarterEdge();
//         qe->index = i;
//         quarter_edges_.emplace_back(qe);
//     }
//     for (u32 i = 0; i < n_quarter_edges; i++) {
//         QuarterEdge* qe = quarter_edges_[i];

//         u32 vertex_index = *(u32*)(data + offset);
//         if (vertex_index != std::numeric_limits<u32>::max()) {
//             qe->vertex = vertices_[vertex_index];
//         } else {
//             qe->vertex = nullptr;
//         }

//         offset += sizeof(u32);
//         qe->next = quarter_edges_[*(u32*)(data + offset)];
//         offset += sizeof(u32);
//         qe->rot = quarter_edges_[*(u32*)(data + offset)];
//         offset += sizeof(u32);
//     }

//     return true;
// }

// //
// ------------------------------------------------------------------------------------------------
// QuarterEdge* DelaunayMesh::Next(const QuarterEdge* qe) const { return qe->next; }

// //
// ------------------------------------------------------------------------------------------------
// QuarterEdge* DelaunayMesh::Rot(const QuarterEdge* qe) const { return qe->rot; }

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex DelaunayMesh::Sym(QuarterEdgeIndex qe) const {
    return quarter_edges_[quarter_edges_[qe.i].i_rot.i].i_rot;
}

// //
// ------------------------------------------------------------------------------------------------
// QuarterEdge* DelaunayMesh::Tor(const QuarterEdge* qe) const { return qe->rot->rot->rot; }

// //
// ------------------------------------------------------------------------------------------------
// QuarterEdge* DelaunayMesh::Prev(const QuarterEdge* qe) const { return qe->rot->next->rot; }

// //
// ------------------------------------------------------------------------------------------------
// QuarterEdge* DelaunayMesh::Lnext(const QuarterEdge* qe) const { return Tor(qe)->next->rot; }

// //
// ------------------------------------------------------------------------------------------------
// bool DelaunayMesh::IsBoundaryVertex(const common::Vec2f* v) const {
//     return v == &(vertices_[0]->vertex) || v == &(vertices_[1]->vertex) ||
//            v == &(vertices_[2]->vertex);
// }

// //
// ------------------------------------------------------------------------------------------------
// bool DelaunayMesh::IsBoundaryVertex(const VertexData& vertex_data) const {
//     return vertex_data.index <= 2;
// }
// bool DelaunayMesh::IsBoundaryVertex(const VertexData* const& vertex_data) const {
//     return vertex_data->index <= 2;
// }

// ------------------------------------------------------------------------------------------------
VertexIndex DelaunayMesh::AddVertex(float x, float y) {
    VertexIndex i_vertex = {kInvalidIndex};
    if (IsValid(i_vertex_free_first_)) {
        // Grab the first item off the free list and use that.
        i_vertex = i_vertex_free_first_;

        VertexData& vertex_data = Get(i_vertex);
        i_vertex_free_first_ = vertex_data.i_next;
        vertex_data.v = common::Vec2f(x, y);
        vertex_data.i_next = {kInvalidIndex};
        vertex_data.i_prev = i_vertex_alive_last_;

    } else {
        // If there is no first item in the free list, append a new item
        i_vertex = {vertices_.size()};
        VertexData vertex_data = {};
        vertex_data.v = common::Vec2f(x, y);
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
QuarterEdgeIndex DelaunayMesh::PrepareFreeQuarterEdge() {
    QuarterEdgeIndex i_qe = {kInvalidIndex};
    if (IsValid(i_qe_free_first_)) {
        // Grab the first item off the free list and use that.
        i_qe = i_qe_free_first_;

        QuarterEdge& qe = Get(i_qe);
        i_qe_free_first_ = qe.i_next;
        qe.i_next = {kInvalidIndex};
        qe.i_prev = i_qe_alive_last_;

    } else {
        // If there is no first item in the free list, append a new item
        i_qe = {quarter_edges_.size()};
        QuarterEdge qe = {};
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
QuarterEdgeIndex DelaunayMesh::AddEdge(VertexIndex i_vertex_a, VertexIndex i_vertex_b) {
    // Add them all first, so we avoid issues with memory changing as we push new quarter edges.
    QuarterEdgeIndex i_ab = PrepareFreeQuarterEdge();
    QuarterEdgeIndex i_lr = PrepareFreeQuarterEdge();
    QuarterEdgeIndex i_ba = PrepareFreeQuarterEdge();
    QuarterEdgeIndex i_rl = PrepareFreeQuarterEdge();

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

// //
// ------------------------------------------------------------------------------------------------
// QuarterEdge* DelaunayMesh::GetEnclosingTriangle(const common::Vec2f& p,
//                                                 QuarterEdge* qe_dual) const {
//     // This should always terminate
//     constexpr int kMaxIters = 100;
//     for (int iter = 0; iter < kMaxIters; iter++) {
//         QuarterEdge* qe_ab = qe_dual->rot;
//         QuarterEdge* qe_bc = qe_dual->next->rot;
//         QuarterEdge* qe_ca = qe_dual->next->next->rot;

//         const auto& a = (qe_ab->vertex)->vertex;
//         const auto& b = (qe_bc->vertex)->vertex;
//         const auto& c = (qe_ca->vertex)->vertex;

//         if (common::GetRightHandedness(a, b, p) < 0) {
//             qe_dual = Rot(qe_ab);  // Move across AB
//         } else if (common::GetRightHandedness(b, c, p) < 0) {
//             qe_dual = Rot(qe_bc);  // Move across BC
//         } else if (common::GetRightHandedness(c, a, p) < 0) {
//             qe_dual = Rot(qe_ca);  // Move across CA
//         } else {
//             return qe_dual;
//         }
//     }
//     return qe_dual;
// }

// //
// ------------------------------------------------------------------------------------------------
// QuarterEdge* DelaunayMesh::GetEnclosingTriangle(const common::Vec2f& p) const {
//     return GetEnclosingTriangle(p, Tor(quarter_edges_.front()));
// }

// //
// ------------------------------------------------------------------------------------------------
// QuarterEdge* DelaunayMesh::GetTriangleQuarterEdge1(const QuarterEdge* qe_dual) const {
//     return qe_dual->rot;
// }

// //
// ------------------------------------------------------------------------------------------------
// QuarterEdge* DelaunayMesh::GetTriangleQuarterEdge2(const QuarterEdge* qe_dual) const {
//     return qe_dual->next->rot;
// }

// //
// ------------------------------------------------------------------------------------------------
// QuarterEdge* DelaunayMesh::GetTriangleQuarterEdge3(const QuarterEdge* qe_dual) const {
//     return qe_dual->next->next->rot;
// }

// //
// ------------------------------------------------------------------------------------------------
// const common::Vec2f& DelaunayMesh::GetTriangleVertex1(const QuarterEdge* qe_dual) const {
//     return GetTriangleQuarterEdge1(qe_dual)->vertex->vertex;
// }

// //
// ------------------------------------------------------------------------------------------------
// const common::Vec2f& DelaunayMesh::GetTriangleVertex2(const QuarterEdge* qe_dual) const {
//     return GetTriangleQuarterEdge2(qe_dual)->vertex->vertex;
// }

// //
// ------------------------------------------------------------------------------------------------
// const common::Vec2f& DelaunayMesh::GetTriangleVertex3(const QuarterEdge* qe_dual) const {
//     return GetTriangleQuarterEdge3(qe_dual)->vertex->vertex;
// }

// //
// ------------------------------------------------------------------------------------------------
// void DelaunayMesh::FlipEdge(QuarterEdge* qe) {
//     QuarterEdge* qe_sym = Sym(qe);
//     QuarterEdge* qe_a = Prev(qe);
//     QuarterEdge* qe_b = Prev(qe_sym);

//     Splice(qe, qe_a);
//     Splice(qe_sym, qe_b);
//     Splice(qe, Lnext(qe_a));
//     Splice(qe_sym, Lnext(qe_b));

//     qe->vertex = Sym(qe_a)->vertex;
//     qe_sym->vertex = Sym(qe_b)->vertex;
// }

// //
// ------------------------------------------------------------------------------------------------
// int DelaunayMesh::AddDelaunayVertex(const common::Vec2f& p) {
//     // Ensure that the triangle is within bounds
//     if (common::Norm(p) >= bounding_radius_) {
//         return kInvalidIndex;
//     }

//     // Get the enclosing triangle
//     QuarterEdge* qe_dual = GetEnclosingTriangle(p);
//     QuarterEdge* qe_ab = qe_dual->rot;
//     QuarterEdge* qe_bc = qe_dual->next->rot;
//     QuarterEdge* qe_ca = qe_dual->next->next->rot;

//     // Grab the vertices
//     const common::Vec2f& a = qe_ab->vertex->vertex;
//     const common::Vec2f& b = qe_bc->vertex->vertex;
//     const common::Vec2f& c = qe_ca->vertex->vertex;

//     // If we are too close to an existing vertex, do nothing.
//     if (std::min({common::Norm(a - p), common::Norm(b - p), common::Norm(c - p)}) <
//         min_dist_to_vertex_) {
//         // Return an already existing vertex
//         // common::Vec2f* v = nullptr;
//         // if (common::Norm(a - p) < min_dist_to_vertex_) {
//         //     v = qe_ab->vertex;
//         // } else if (common::Norm(b - p) < min_dist_to_vertex_) {
//         //     v = qe_bc->vertex;
//         // } else {
//         //     v = qe_ca->vertex;
//         // }

//         // for (int i = 0; i < (int)NumVertices(); i++) {
//         //     if (vertices_[i] == v) {
//         //         return i;
//         //     }
//         // }

//         return kInvalidIndex;
//     }

//     // Add the vertex
//     VertexData* p_ptr = AddVertex(p.x, p.y);

//     // Check whether we are (effectively) on an edge.
//     // If we are, we delete the existing edge and add four instead.
//     float dist_to_ab = common::GetDistanceToLine(p, a, b);
//     float dist_to_bc = common::GetDistanceToLine(p, b, c);
//     float dist_to_ca = common::GetDistanceToLine(p, c, a);
//     QuarterEdge* qe_start = nullptr;
//     if (std::min({dist_to_ab, dist_to_bc, dist_to_ca}) > min_dist_to_edge_) {
//         // Normal case. We are not too close to an edge.

//         // Add the new edges and splice them in
//         QuarterEdge* ap = AddEdge(qe_ab->vertex, p_ptr);
//         QuarterEdge* bp = AddEdge(qe_bc->vertex, p_ptr);
//         QuarterEdge* cp = AddEdge(qe_ca->vertex, p_ptr);

//         Splice(ap, qe_ab);
//         Splice(bp, qe_bc);
//         Splice(cp, qe_ca);

//         // TODO: Figure out splice such that we get this desired effect.
//         //       i.e. can we replace these next calls with three calls to splice?
//         QuarterEdge* pa = Sym(ap);
//         QuarterEdge* pb = Sym(bp);
//         QuarterEdge* pc = Sym(cp);

//         pa->next = bp->rot->rot;  // sym(bp)
//         pb->next = cp->rot->rot;  // sym(bc)
//         pc->next = ap->rot->rot;  // sym(ap)

//         pa->rot->next = cp->rot;
//         pb->rot->next = ap->rot;
//         pc->rot->next = bp->rot;

//         // Set our start quarter edge
//         qe_start = pa;
//     } else {
//         // We are effectively on an edge.
//         // Identify that edge, and cut it with the new vertex.

//         // Let DE be the edge we are on, F be the far vertex qe, and G be the vertex qe across DE
//         // from F.
//         QuarterEdge* qe_dp = nullptr;
//         QuarterEdge* qe_ef = nullptr;
//         QuarterEdge* qe_fd = nullptr;
//         if (dist_to_ab < min_dist_to_edge_) {
//             qe_dp = qe_ab;
//             qe_ef = qe_bc;
//             qe_fd = qe_ca;
//         } else if (dist_to_bc < min_dist_to_edge_) {
//             qe_dp = qe_bc;
//             qe_ef = qe_ca;
//             qe_fd = qe_ab;
//         } else {
//             qe_dp = qe_ca;
//             qe_ef = qe_ab;
//             qe_fd = qe_bc;
//         }
//         QuarterEdge* qe_ge = Prev(Sym(Prev(qe_dp)));

//         // Our edge may not be the boundary edge
//         int num_boundary_vertices =
//             IsBoundaryVertex(qe_dp->vertex) + IsBoundaryVertex(Sym(qe_dp)->vertex);
//         if (num_boundary_vertices == 2) {
//             return kInvalidIndex;
//         }

//         // Grab another quarter edge we will need
//         QuarterEdge* qe_dg = Prev(qe_dp);

//         // Reset the things that point to DE
//         qe_dg->next = Rot(qe_fd)->rot;  // sym(fd)
//         QuarterEdge* qe_ge_tor = Tor(qe_ge);
//         qe_ge_tor->next = Sym(qe_ef)->rot;
//         qe_ef->next = Rot(qe_ge)->rot;  // sym(ge)
//         QuarterEdge* qe_fd_tor = Tor(qe_fd);
//         qe_fd_tor->next = Sym(qe_dg)->rot;

//         // Reset DP as a quarter edge, and change it to DP
//         QuarterEdge* qe_pd = Sym(qe_dp);
//         qe_pd->vertex = p_ptr;
//         {
//             qe_pd->next = qe_pd;
//             qe_pd->rot->next = qe_dp->rot;
//             qe_dp->next = qe_dp;
//             qe_dp->rot->next = qe_pd->rot;
//         }

//         // Create three new edges EP, FP, and F->P
//         QuarterEdge* qe_ep = AddEdge(qe_ef->vertex, p_ptr);
//         QuarterEdge* qe_fp = AddEdge(qe_fd->vertex, p_ptr);
//         QuarterEdge* qe_gp = AddEdge(qe_ge->vertex, p_ptr);

//         // Splice them all correctly
//         Splice(qe_dp, qe_dg);
//         Splice(qe_gp, qe_ge);
//         Splice(qe_ep, qe_ef);
//         Splice(qe_fp, qe_fd);

//         // TODO: Figure out splice such that we get this desired effect.
//         //       i.e. can we replace these next calls with three calls to splice?
//         qe_pd = Sym(qe_dp);
//         QuarterEdge* qe_pe = Sym(qe_ep);
//         QuarterEdge* qe_pf = Sym(qe_fp);
//         QuarterEdge* qe_pg = Sym(qe_gp);

//         qe_pd->next = Rot(qe_gp)->rot;  // sym(gp)
//         qe_pe->next = Rot(qe_fp)->rot;  // sym(fp)
//         qe_pf->next = Rot(qe_dp)->rot;  // sym(dp)
//         qe_pg->next = Rot(qe_ep)->rot;  // sym(ep)

//         Rot(qe_pd)->next = qe_fp->rot;
//         Rot(qe_pe)->next = qe_gp->rot;
//         Rot(qe_pf)->next = qe_ep->rot;
//         Rot(qe_pg)->next = qe_dp->rot;

//         // Set our start quarter edge
//         qe_start = qe_pd;
//     }

//     // Check if we are locally Delaunay.
//     // Walk around the outer edges and flip any edges that are not locally delaunay.
//     // We can check an edge by seeing if the opposite vertex is within the inscribed circle
//     // of the inner vertex + edge vertices.

//     // We start at pa and rotate around it to get all edges that we have to check.
//     // We need to walk the outer edges until we get back to pa.
//     // Each outer edge is given by next(mesh, sym(mesh, qe))

//     QuarterEdge* qe = qe_start;
//     bool done = false;
//     while (!done) {
//         QuarterEdge* qe_outer_edge = Sym(qe)->next;

//         // Advance
//         qe = qe->next;
//         done = qe == qe_start;

//         // Only consider the edge if it is not a bounding edge
//         VertexData* src_ptr = qe_outer_edge->vertex;
//         VertexData* dst_ptr = Sym(qe_outer_edge)->vertex;
//         int num_boundary_vertices = IsBoundaryVertex(src_ptr) + IsBoundaryVertex(dst_ptr);
//         if (num_boundary_vertices == 2) {  // one edge being on the boundary is okay
//             continue;
//         }

//         // Check the edge from qe_outer_edge to sym(qe_outer_edge)
//         const common::Vec2f& src = src_ptr->vertex;
//         const common::Vec2f& dst = dst_ptr->vertex;

//         // Get the far vertex across the dividing edge
//         VertexData* far_ptr = Sym(qe_outer_edge->next)->vertex;
//         const common::Vec2f& far = far_ptr->vertex;

//         // If the edge contains a boundary vertex, don't flip it if it would produce an
//         inside-out
//         // triangle.
//         if (num_boundary_vertices == 1) {
//             if (common::GetTriangleContainment(dst, far, p, src) >= 0 ||
//                 common::GetTriangleContainment(dst, far, src, p) >= 0 ||
//                 common::GetTriangleContainment(src, far, dst, p) >= 0 ||
//                 common::GetTriangleContainment(src, far, p, dst) >= 0) {
//                 continue;
//             }
//         }

//         if (common::GetCircleContainment(p, src, dst, far) > 0 ||
//             common::GetCircleContainment(far, p, dst, src) > 0) {
//             // Either p is inside the circle passing through src, dst, and far, or
//             //  far is inside the circle passing through p, src, and dst.
//             // We have to flip the edge.
//             FlipEdge(qe_outer_edge);

//             // We flipped the edge, so qe_outer_edge has to be traversed again.
//             // Back it up.
//             qe = Prev(qe);
//             if (qe != qe_start) {
//                 qe = Prev(qe);
//             }
//             done = false;
//         }
//     }

//     // Return the index of the newly added vertex
//     return vertices_.size() - 1;
// }

// //
// ------------------------------------------------------------------------------------------------
// QuarterEdge* DelaunayMesh::GetQuarterEdge(int i, int j) const {
//     if (i == j) {
//         return nullptr;  // Self edges do not exist
//     }

//     int n = NumVertices();
//     if (i < 0 || i >= n || j < 0 || j >= n) {
//         return nullptr;  // Invalid index
//     }

//     // If we already have an edge between these two vertices, we are done.
//     VertexData* a = vertices_.at(i);
//     VertexData* b = vertices_.at(j);
//     for (QuarterEdge* qe : quarter_edges_) {
//         if (qe->vertex == a) {
//             if (Sym(qe)->vertex == b) {
//                 return qe;
//             }
//         }
//     }

//     return nullptr;
// }

// //
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

// //
// ------------------------------------------------------------------------------------------------
// bool DelaunayMesh::ConstrainEdge(int i, int j) {
//     // TODO: If we have vertices intersecting the edge, we should split our constrained edge
//     //       by those vertices and call this for each sub-edge.

//     // TODO: If we try to constrain an edge that overlaps with an already-constrained edge,
//     //       then we need to produce an intersection and fix it that way.

//     if (i == j) {
//         return false;  // Cannot add self-edges
//     }

//     int n = NumVertices();
//     if (i < 0 || i >= n || j < 0 || j >= n) {
//         return false;  // Invalid index
//     }

//     // If we already have an edge between these two vertices, we are done.
//     VertexData* a_data = vertices_.at(i);
//     VertexData* b_data = vertices_.at(j);
//     QuarterEdge* qe_a = nullptr;
//     for (QuarterEdge* qe : quarter_edges_) {
//         if (qe->vertex == a_data) {
//             qe_a = qe;
//             if (Sym(qe)->vertex == b_data) {
//                 return true;  // already exists
//             }
//         }
//     }

//     const common::Vec2f& a = a_data->vertex;
//     const common::Vec2f& b = b_data->vertex;

//     std::cout << "A: (" << a.x << ", " << a.y << ")" << std::endl;
//     std::cout << "B: (" << b.x << ", " << b.y << ")" << std::endl;

//     // --------------------------------------------------------------------------
//     // The edge does not yet exist.
//     // Walk around, from A, toward B and flip any offending edges.
//     // Repeat until we have flipped our way to producing AB.

//     for (int iter = 0; iter < 100; iter++) {
//         // Rotate qe_a to the last coincident quarter edge that is CCW of the new segment.
//         while (GetRightHandedness(a, Sym(qe_a)->vertex->vertex, b) <= 0.0) {
//             qe_a = qe_a->next;
//         }
//         while (GetRightHandedness(a, Sym(qe_a->next)->vertex->vertex, b) > 0.0) {
//             qe_a = qe_a->next;
//         }

//         // qe_a will then have a CCW triangle ACD.
//         // If the far side of the triangle is B (ACD == ADB), then we are done.
//         // Otherwise, D is on the other side of AB, so CD intersects AB, and we need to see if we
//         // can flip CD.
//         QuarterEdge* qe_dual =
//             Tor(qe_a);  // The dual quarter edge that starts inside ACD and points across AC.
//         VertexData* c_data = (qe_dual->next->rot)->vertex;
//         VertexData* d_data = (qe_dual->next->next->rot)->vertex;
//         const common::Vec2f& c = c_data->vertex;
//         const common::Vec2f& d = d_data->vertex;

//         if (d_data == b_data) {
//             // We have B, so we have produced edge ab and are done.
//             return true;
//         } else {
//             // See if we can flip CD. We have to ensure that it is a convex quadrilateral.
//             // E is the vertex on the far side.
//             VertexData* e_data = Tor(Sym(qe_dual->next)->next)->vertex;
//             const common::Vec2f& e = e_data->vertex;
//             if (GetRightHandedness(a, c, e) > 0 && GetRightHandedness(a, e, d) > 0) {
//                 // Flip it
//                 QuarterEdge* qe_cd = qe_dual->next->rot;
//                 FlipEdge(qe_cd);
//             } else {
//                 // TODO: Handle this case. We need to progress past C and try the next triangle
//                 that
//                 //       could overlap.
//                 return false;
//             }
//         }
//     }

//     // AAAAAAH this should never happen.
//     return false;
// }

}  // namespace mesh