#include "delaunay_mesh.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace core {

// ------------------------------------------------------------------------------------------------
bool IsDualEdge(const QuarterEdge& qe) { return qe.vertex == nullptr; }

// ------------------------------------------------------------------------------------------------
bool IsPrimalEdge(const QuarterEdge& qe) { return qe.vertex != nullptr; }

// ------------------------------------------------------------------------------------------------
void SwapNexts(QuarterEdge* a, QuarterEdge* b) {
    auto tmp = a->next;
    a->next = b->next;
    b->next = tmp;
}

// ------------------------------------------------------------------------------------------------
DelaunayMesh::DelaunayMesh(float bounding_radius_, float min_dist_to_vertex,
                           float min_dist_to_edge) :
    bounding_radius_(bounding_radius_),
    min_dist_to_vertex_(min_dist_to_vertex),
    min_dist_to_edge_(min_dist_to_edge) {
    // The triangle radius is 2r + eps(), which guarantees that it is large enough.
    float r = 2 * bounding_radius_ + min_dist_to_edge + min_dist_to_vertex;

    common::Vec2f* a = AddVertex(r * std::cos(90 * M_PI / 180.0), r * std::sin(90 * M_PI / 180.0));
    common::Vec2f* b =
        AddVertex(r * std::cos(210 * M_PI / 180.0), r * std::sin(210 * M_PI / 180.0));
    common::Vec2f* c =
        AddVertex(r * std::cos(-30 * M_PI / 180.0), r * std::sin(-30 * M_PI / 180.0));

    QuarterEdge* ab = AddEdge(a, b);
    QuarterEdge* bc = AddEdge(b, c);
    QuarterEdge* ca = AddEdge(c, a);

    Splice(Sym(ab), bc);
    Splice(Sym(bc), ca);
    Splice(Sym(ca), ab);
}

// ------------------------------------------------------------------------------------------------
DelaunayMesh::~DelaunayMesh() {
    for (auto it = vertices_.begin(); it != vertices_.end(); ++it) {
        delete *it;
    }
    for (auto it = quarter_edges_.begin(); it != quarter_edges_.end(); ++it) {
        delete *it;
    }
}

// ------------------------------------------------------------------------------------------------
QuarterEdge* DelaunayMesh::Next(const QuarterEdge* qe) { return qe->next; }

// ------------------------------------------------------------------------------------------------
QuarterEdge* DelaunayMesh::Rot(const QuarterEdge* qe) { return qe->rot; }

// ------------------------------------------------------------------------------------------------
QuarterEdge* DelaunayMesh::Sym(const QuarterEdge* qe) { return qe->rot->rot; }

// ------------------------------------------------------------------------------------------------
QuarterEdge* DelaunayMesh::Tor(const QuarterEdge* qe) { return qe->rot->rot->rot; }

// ------------------------------------------------------------------------------------------------
QuarterEdge* DelaunayMesh::Prev(const QuarterEdge* qe) { return qe->rot->next->rot; }

// ------------------------------------------------------------------------------------------------
QuarterEdge* DelaunayMesh::Lnext(const QuarterEdge* qe) { return Tor(qe)->next->rot; }

// ------------------------------------------------------------------------------------------------
bool DelaunayMesh::IsBoundaryVertex(const common::Vec2f* v) {
    return v == vertices_[0] || v == vertices_[1] || v == vertices_[2];
}

// ------------------------------------------------------------------------------------------------
common::Vec2f* DelaunayMesh::AddVertex(float x, float y) {
    common::Vec2f* vertex = new common::Vec2f({.x = x, .y = y});
    vertices_.push_back(vertex);
    return vertex;
}

// ------------------------------------------------------------------------------------------------
QuarterEdge* DelaunayMesh::AddEdge(common::Vec2f* a, common::Vec2f* b) {
    int base_index = quarter_edges_.size();
    QuarterEdge* ab = new QuarterEdge({base_index, a, nullptr, nullptr});
    QuarterEdge* lr = new QuarterEdge({base_index + 1, nullptr, nullptr, nullptr});
    QuarterEdge* ba = new QuarterEdge({base_index + 2, b, nullptr, nullptr});
    QuarterEdge* rl = new QuarterEdge({base_index + 3, nullptr, nullptr, nullptr});

    ab->next = ab;
    ab->rot = lr;
    lr->next = rl;
    lr->rot = ba;
    ba->next = ba;
    ba->rot = rl;
    rl->next = lr;
    rl->rot = ab;

    quarter_edges_.push_back(ab);
    quarter_edges_.push_back(lr);
    quarter_edges_.push_back(ba);
    quarter_edges_.push_back(rl);

    return ab;
}

// ------------------------------------------------------------------------------------------------
void DelaunayMesh::Splice(QuarterEdge* a, QuarterEdge* b) {
    SwapNexts(a->next->rot, b->next->rot);
    SwapNexts(a, b);
}

// ------------------------------------------------------------------------------------------------
QuarterEdge* DelaunayMesh::GetEnclosingTriangle(const common::Vec2f& p, QuarterEdge* qe_dual) {
    // This should always terminate
    constexpr int kMaxIters = 100;
    for (int iter = 0; iter < kMaxIters; iter++) {
        QuarterEdge* qe_ab = qe_dual->rot;
        QuarterEdge* qe_bc = qe_dual->next->rot;
        QuarterEdge* qe_ca = qe_dual->next->next->rot;

        const auto& a = *(qe_ab->vertex);
        const auto& b = *(qe_bc->vertex);
        const auto& c = *(qe_ca->vertex);

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
QuarterEdge* DelaunayMesh::GetEnclosingTriangle(const common::Vec2f& p) {
    return GetEnclosingTriangle(p, Tor(quarter_edges_.front()));
}

// ------------------------------------------------------------------------------------------------
const common::Vec2f& DelaunayMesh::GetTriangleVertex1(const QuarterEdge* qe_dual) {
    return *((qe_dual->rot)->vertex);
}

// ------------------------------------------------------------------------------------------------
const common::Vec2f& DelaunayMesh::GetTriangleVertex2(const QuarterEdge* qe_dual) {
    return *((qe_dual->next->rot)->vertex);
}

// ------------------------------------------------------------------------------------------------
const common::Vec2f& DelaunayMesh::GetTriangleVertex3(const QuarterEdge* qe_dual) {
    return *((qe_dual->next->next->rot)->vertex);
}

// ------------------------------------------------------------------------------------------------
void DelaunayMesh::FlipEdge(QuarterEdge* qe) {
    QuarterEdge* qe_sym = Sym(qe);
    QuarterEdge* qe_a = Prev(qe);
    QuarterEdge* qe_b = Prev(qe_sym);

    Splice(qe, qe_a);
    Splice(qe_sym, qe_b);
    Splice(qe, Lnext(qe_a));
    Splice(qe_sym, Lnext(qe_b));

    qe->vertex = Sym(qe_a)->vertex;
    qe_sym->vertex = Sym(qe_b)->vertex;
}

// ------------------------------------------------------------------------------------------------
int DelaunayMesh::AddDelaunayVertex(const common::Vec2f& p) {
    // Ensure that the triangle is within bounds
    if (common::Norm(p) >= bounding_radius_) {
        return kInvalidIndex;
    }

    // Get the enclosing triangle
    QuarterEdge* qe_dual = GetEnclosingTriangle(p);
    QuarterEdge* qe_ab = qe_dual->rot;
    QuarterEdge* qe_bc = qe_dual->next->rot;
    QuarterEdge* qe_ca = qe_dual->next->next->rot;

    // Grab the vertices
    const common::Vec2f& a = *(qe_ab->vertex);
    const common::Vec2f& b = *(qe_bc->vertex);
    const common::Vec2f& c = *(qe_ca->vertex);

    // If we are too close to an existing vertex, do nothing.
    if (std::min({common::Norm(a - p), common::Norm(b - p), common::Norm(c - p)}) <
        min_dist_to_vertex_) {
        // Return an already existing vertex
        // common::Vec2f* v = nullptr;
        // if (common::Norm(a - p) < min_dist_to_vertex_) {
        //     v = qe_ab->vertex;
        // } else if (common::Norm(b - p) < min_dist_to_vertex_) {
        //     v = qe_bc->vertex;
        // } else {
        //     v = qe_ca->vertex;
        // }

        // for (int i = 0; i < (int)NumVertices(); i++) {
        //     if (vertices_[i] == v) {
        //         return i;
        //     }
        // }

        return kInvalidIndex;
    }

    // Add the vertex
    common::Vec2f* p_ptr = AddVertex(p.x, p.y);

    // Check whether we are (effectively) on an edge.
    // If we are, we delete the existing edge and add four instead.
    float dist_to_ab = common::GetDistanceToLine(p, a, b);
    float dist_to_bc = common::GetDistanceToLine(p, b, c);
    float dist_to_ca = common::GetDistanceToLine(p, c, a);
    QuarterEdge* qe_start = nullptr;
    if (std::min({dist_to_ab, dist_to_bc, dist_to_ca}) > min_dist_to_edge_) {
        // Normal case. We are not too close to an edge.

        // Add the new edges and splice them in
        QuarterEdge* ap = AddEdge(qe_ab->vertex, p_ptr);
        QuarterEdge* bp = AddEdge(qe_bc->vertex, p_ptr);
        QuarterEdge* cp = AddEdge(qe_ca->vertex, p_ptr);

        Splice(ap, qe_ab);
        Splice(bp, qe_bc);
        Splice(cp, qe_ca);

        // TODO: Figure out splice such that we get this desired effect.
        //       i.e. can we replace these next calls with three calls to splice?
        QuarterEdge* pa = Sym(ap);
        QuarterEdge* pb = Sym(bp);
        QuarterEdge* pc = Sym(cp);

        pa->next = bp->rot->rot;  // sym(bp)
        pb->next = cp->rot->rot;  // sym(bc)
        pc->next = ap->rot->rot;  // sym(ap)

        pa->rot->next = cp->rot;
        pb->rot->next = ap->rot;
        pc->rot->next = bp->rot;

        // Set our start quarter edge
        qe_start = pa;
    } else {
        // We are effectively on an edge.
        // Identify that edge, and cut it with the new vertex.

        // Let DE be the edge we are on, F be the far vertex qe, and G be the vertex qe across DE
        // from F.
        QuarterEdge* qe_dp = nullptr;
        QuarterEdge* qe_ef = nullptr;
        QuarterEdge* qe_fd = nullptr;
        if (dist_to_ab < min_dist_to_edge_) {
            qe_dp = qe_ab;
            qe_ef = qe_bc;
            qe_fd = qe_ca;
        } else if (dist_to_bc < min_dist_to_edge_) {
            qe_dp = qe_bc;
            qe_ef = qe_ca;
            qe_fd = qe_ab;
        } else {
            qe_dp = qe_ca;
            qe_ef = qe_ab;
            qe_fd = qe_bc;
        }
        QuarterEdge* qe_ge = Prev(Sym(Prev(qe_dp)));

        // Our edge may not be the boundary edge
        int num_boundary_vertices =
            IsBoundaryVertex(qe_dp->vertex) + IsBoundaryVertex(Sym(qe_dp)->vertex);
        if (num_boundary_vertices == 2) {
            return kInvalidIndex;
        }

        // Grab another quarter edge we will need
        QuarterEdge* qe_dg = Prev(qe_dp);

        // Reset the things that point to DE
        qe_dg->next = Rot(qe_fd)->rot;  // sym(fd)
        QuarterEdge* qe_ge_tor = Tor(qe_ge);
        qe_ge_tor->next = Sym(qe_ef)->rot;
        qe_ef->next = Rot(qe_ge)->rot;  // sym(ge)
        QuarterEdge* qe_fd_tor = Tor(qe_fd);
        qe_fd_tor->next = Sym(qe_dg)->rot;

        // Reset DP as a quarter edge, and change it to DP
        QuarterEdge* qe_pd = Sym(qe_dp);
        qe_pd->vertex = p_ptr;
        {
            qe_pd->next = qe_pd;
            qe_pd->rot->next = qe_dp->rot;
            qe_dp->next = qe_dp;
            qe_dp->rot->next = qe_pd->rot;
        }

        // Create three new edges EP, FP, and F->P
        QuarterEdge* qe_ep = AddEdge(qe_ef->vertex, p_ptr);
        QuarterEdge* qe_fp = AddEdge(qe_fd->vertex, p_ptr);
        QuarterEdge* qe_gp = AddEdge(qe_ge->vertex, p_ptr);

        // Splice them all correctly
        Splice(qe_dp, qe_dg);
        Splice(qe_gp, qe_ge);
        Splice(qe_ep, qe_ef);
        Splice(qe_fp, qe_fd);

        // TODO: Figure out splice such that we get this desired effect.
        //       i.e. can we replace these next calls with three calls to splice?
        qe_pd = Sym(qe_dp);
        QuarterEdge* qe_pe = Sym(qe_ep);
        QuarterEdge* qe_pf = Sym(qe_fp);
        QuarterEdge* qe_pg = Sym(qe_gp);

        qe_pd->next = Rot(qe_gp)->rot;  // sym(gp)
        qe_pe->next = Rot(qe_fp)->rot;  // sym(fp)
        qe_pf->next = Rot(qe_dp)->rot;  // sym(dp)
        qe_pg->next = Rot(qe_ep)->rot;  // sym(ep)

        Rot(qe_pd)->next = qe_fp->rot;
        Rot(qe_pe)->next = qe_gp->rot;
        Rot(qe_pf)->next = qe_ep->rot;
        Rot(qe_pg)->next = qe_dp->rot;

        // Set our start quarter edge
        qe_start = qe_pd;
    }

    // Check if we are locally Delaunay.
    // Walk around the outer edges and flip any edges that are not locally delaunay.
    // We can check an edge by seeing if the opposite vertex is within the inscribed circle
    // of the inner vertex + edge vertices.

    // We start at pa and rotate around it to get all edges that we have to check.
    // We need to walk the outer edges until we get back to pa.
    // Each outer edge is given by next(mesh, sym(mesh, qe))

    QuarterEdge* qe = qe_start;
    bool done = false;
    while (!done) {
        QuarterEdge* qe_outer_edge = Sym(qe)->next;

        // Advance
        qe = qe->next;
        done = qe == qe_start;

        // Only consider the edge if it is not a bounding edge
        common::Vec2f* src_ptr = qe_outer_edge->vertex;
        common::Vec2f* dst_ptr = Sym(qe_outer_edge)->vertex;
        int num_boundary_vertices = IsBoundaryVertex(src_ptr) + IsBoundaryVertex(dst_ptr);
        if (num_boundary_vertices == 2) {  // one edge being on the boundary is okay
            continue;
        }

        // Check the edge from qe_outer_edge to sym(qe_outer_edge)
        const common::Vec2f& src = *src_ptr;
        const common::Vec2f& dst = *dst_ptr;

        // Get the far vertex across the dividing edge
        common::Vec2f* far_ptr = Sym(qe_outer_edge->next)->vertex;
        const common::Vec2f& far = *far_ptr;

        // If the edge contains a boundary vertex, don't flip it if it would produce an inside-out
        // triangle.
        if (num_boundary_vertices == 1) {
            if (common::GetTriangleContainment(dst, far, p, src) >= 0 ||
                common::GetTriangleContainment(dst, far, src, p) >= 0 ||
                common::GetTriangleContainment(src, far, dst, p) >= 0 ||
                common::GetTriangleContainment(src, far, p, dst) >= 0) {
                continue;
            }
        }

        if (common::GetCircleContainment(p, src, dst, far) > 0 ||
            common::GetCircleContainment(far, p, dst, src) > 0) {
            // Either p is inside the circle passing through src, dst, and far, or
            //  far is inside the circle passing through p, src, and dst.
            // We have to flip the edge.
            FlipEdge(qe_outer_edge);

            // We flipped the edge, so qe_outer_edge has to be traversed again.
            // Back it up.
            qe = Prev(qe);
            if (qe != qe_start) {
                qe = Prev(qe);
            }
            done = false;
        }
    }

    // Return the index of the newly added vertex
    return vertices_.size() - 1;
}

// ------------------------------------------------------------------------------------------------
QuarterEdge* DelaunayMesh::GetQuarterEdge(int i, int j) {
    if (i == j) {
        return nullptr;  // Self edges do not exist
    }

    int n = NumVertices();
    if (i < 0 || i >= n || j < 0 || j >= n) {
        return nullptr;  // Invalid index
    }

    // If we already have an edge between these two vertices, we are done.
    common::Vec2f* a = vertices_.at(i);
    common::Vec2f* b = vertices_.at(j);
    for (QuarterEdge* qe : quarter_edges_) {
        if (qe->vertex == a) {
            if (Sym(qe)->vertex == b) {
                return qe;
            }
        }
    }

    return nullptr;
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMesh::HasEdge(int i, int j) {
    if (i == j) {
        return false;  // Self edges do not exist
    }

    int n = NumVertices();
    if (i < 0 || i >= n || j < 0 || j >= n) {
        return false;  // Invalid index
    }

    // If we already have an edge between these two vertices, we are done.
    common::Vec2f* a = vertices_.at(i);
    common::Vec2f* b = vertices_.at(j);
    for (QuarterEdge* qe : quarter_edges_) {
        if (qe->vertex == a) {
            if (Sym(qe)->vertex == b) {
                return true;
            }
        }
    }

    return false;
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMesh::ConstrainEdge(int i, int j) {
    // TODO: If we have vertices intersecting the edge, we should split our constrained edge
    //       by those vertices and call this for each sub-edge.

    // TODO: If we try to constrain an edge that overlaps with an already-constrained edge,
    //       then we need to produce an intersection and fix it that way.

    if (i == j) {
        return false;  // Cannot add self-edges
    }

    int n = NumVertices();
    if (i < 0 || i >= n || j < 0 || j >= n) {
        return false;  // Invalid index
    }

    // If we already have an edge between these two vertices, we are done.
    common::Vec2f* a = vertices_.at(i);
    common::Vec2f* b = vertices_.at(j);
    QuarterEdge* qe_a = nullptr;
    for (QuarterEdge* qe : quarter_edges_) {
        if (qe->vertex == a) {
            qe_a = qe;
            if (Sym(qe)->vertex == b) {
                return true;  // already exists
            }
        }
    }

    std::cout << "A: (" << a->x << ", " << a->y << ")" << std::endl;
    std::cout << "B: (" << b->x << ", " << b->y << ")" << std::endl;

    // --------------------------------------------------------------------------
    // The edge does not yet exist.
    // Walk around, from A, toward B and flip any offending edges.
    // Repeat until we have flipped our way to producing AB.

    for (int iter = 0; iter < 100; iter++) {
        std::cout << "Iter: " << iter << std::endl;

        // Rotate qe_a to the last coincident quarter edge that is CCW of the new segment.
        while (GetRightHandedness(*a, *(Sym(qe_a)->vertex), *b) <= 0.0) {
            qe_a = qe_a->next;
        }
        while (GetRightHandedness(*a, *(Sym(qe_a->next)->vertex), *b) > 0.0) {
            qe_a = qe_a->next;
        }

        std::cout << "Have qe_a (" << qe_a->vertex->x << ", " << qe_a->vertex->y << ") -> ("
                  << Sym(qe_a)->vertex->x << ", " << Sym(qe_a)->vertex->y << ")" << std::endl;

        // qe_a will then have a CCW triangle ACD.
        // If the far side of the triangle is B (ACD == ADB), then we are done.
        // Otherwise, D is on the other side of AB, so CD intersects AB, and we need to see if we
        // can flip CD.
        QuarterEdge* qe_dual =
            Tor(qe_a);  // The dual quarter edge that starts inside ACD and points across AC.
        common::Vec2f* c = (qe_dual->next->rot)->vertex;
        common::Vec2f* d = (qe_dual->next->next->rot)->vertex;

        std::cout << "C: (" << c->x << ", " << c->y << ")" << std::endl;
        std::cout << "D: (" << d->x << ", " << d->y << ")" << std::endl;

        if (d == b) {
            // We have B, so we have produced edge ab and are done.
            return true;
        } else {
            // See if we can flip CD. We have to ensure that it is a convex quadrilateral.
            // E is the vertex on the far side.
            common::Vec2f* e = Tor(Sym(qe_dual->next)->next)->vertex;
            std::cout << "E: (" << e->x << ", " << e->y << ")" << std::endl;
            if (GetRightHandedness(*a, *c, *e) > 0 && GetRightHandedness(*a, *e, *d) > 0) {
                // Flip it
                std::cout << "Flip it!" << std::endl;
                QuarterEdge* qe_cd = qe_dual->next->rot;
                FlipEdge(qe_cd);
            } else {
                // TODO: Handle this case. We need to progress past C and try the next triangle that
                //       could overlap.
                return false;
            }
        }
    }

    // AAAAAAH this should never happen.
    return false;
}

}  // namespace core