#pragma once

#include <cstddef>
#include <iterator>
#include <limits>
#include <tuple>
#include <vector>

#include "geometry_utils.hpp"
#include "typedefs.hpp"

namespace core {

constexpr usize kInvalidIndex = std::numeric_limits<usize>::max();

struct VertexIndex {
    usize i;

    bool operator==(const VertexIndex& rhs) const { return i == rhs.i; }
    bool operator!=(const VertexIndex& rhs) const { return i != rhs.i; }
    bool operator>=(const VertexIndex& rhs) const { return i >= rhs.i; }
    bool operator>(const VertexIndex& rhs) const { return i > rhs.i; }
    bool operator<=(const VertexIndex& rhs) const { return i <= rhs.i; }
    bool operator<(const VertexIndex& rhs) const { return i < rhs.i; }
};

bool IsValid(VertexIndex i);

struct QuarterEdgeIndex {
    usize i;

    bool operator==(const QuarterEdgeIndex& rhs) const { return i == rhs.i; }
    bool operator!=(const QuarterEdgeIndex& rhs) const { return i != rhs.i; }
    bool operator>=(const QuarterEdgeIndex& rhs) const { return i >= rhs.i; }
    bool operator>(const QuarterEdgeIndex& rhs) const { return i > rhs.i; }
    bool operator<=(const QuarterEdgeIndex& rhs) const { return i <= rhs.i; }
    bool operator<(const QuarterEdgeIndex& rhs) const { return i < rhs.i; }
};

bool IsValid(QuarterEdgeIndex i);

struct VertexData {
    common::Vec2f v;     // Typically set to (NaN,NaN) if free
    VertexIndex i_self;  // Index of itself
    VertexIndex i_next;  // Index of the next alive/free VertexData (typemax(size_t) if none)
    VertexIndex i_prev;  // Index of the previous alive/free VertexData (typemax(size_t) if none)
};

// Whether the given quarter edge is part of a constrained edge, which means we will not flip it.
constexpr u32 QE_FLAG_CONSTRAINED = 1;

// A quarter edge represents a directed edge in a QuadEdgeTree.
// Each undirected edge in the graph A <-> B with face L on the left of A->B and face R on the
// right of A->B has four quarter edges:
//    A->B, L->R, B->A, and R->L.
// If this quarter edge is primal, then the vertex index points to its vertex.
// If this quarter edge is dual (i.e., represents a face), then ivertex is set to
// typemax(size_t).
struct QuarterEdge {
    VertexIndex i_vertex;  // Non-null if this is a quarter edge originating at a vertex
                           // Set to typemax(size_t) if free or dual
    u32 flags;

    QuarterEdgeIndex i_nxt;   // The next right-hand (CCW) QuarterEdge with the same origin
    QuarterEdgeIndex i_rot;   // The next right-hand (CCW) QuarterEdge associated with the same
                              // undirected original edge.
    QuarterEdgeIndex i_self;  // Index of itself
    QuarterEdgeIndex i_next;  // Index of the next alive/free QuarterEdge (typemax(size_t) if none)
    QuarterEdgeIndex i_prev;  // Index of the previous alive/free QuarterEdge (if any)
};

// A dual edge connects two faces.
bool IsDualEdge(const QuarterEdge& qe);

// A primal edge connects two vertices.
bool IsPrimalEdge(const QuarterEdge& qe);

// A DelaunayMesh is a Quad Edge Mesh that can dynamically generate a Delaunay (or constrained
// Delaunay) triangulation. A Quad Edge Mesh's name comes from the fact that each edge in a
// "normal" mesh graph that connects two vertices is represented by 4 edges in the quad edge
// tree - one directed edge between each vertex and one directed edge between the two adjoined
// faces. Because of this, we always have a multiple of 4 quater edges.
//
// This class allows for the insertion and removal of items. As such, our vectors may change in
// length as memory reallocates. We use integers to index into our vertices and quarter edges to
// avoid memeory reallocation issues.
class DelaunayMesh {
  public:
    // Construct a QuadEdgeMesh that is merely a bounding triangle with the given radius.
    DelaunayMesh(float bounding_radius, float min_dist_to_vertex, float min_dist_to_edge);
    ~DelaunayMesh();

    size_t NumVertices() const { return n_vertices_; }
    size_t NumQuarterEdges() const { return n_quarter_edges_; }
    size_t NumEdges() const { return n_quarter_edges_ >> 2; }

    void Clear();
    bool LoadFromData(const u8* data);

    const common::Vec2f& GetVertex(VertexIndex i) const { return vertices_[i.i].v; }
    const common::Vec2f& GetVertex(QuarterEdgeIndex i) const {
        return GetVertex(quarter_edges_[i.i].i_vertex);
    }
    const VertexData& GetVertexData(VertexIndex i) const { return vertices_[i.i]; }
    const VertexData& GetVertexData(QuarterEdgeIndex i) const {
        return GetVertexData(quarter_edges_[i.i].i_vertex);
    }
    const QuarterEdge& GetQuarterEdge(QuarterEdgeIndex i) const { return quarter_edges_[i.i]; }

    // Get any primal quarter edge whose origin is the given vertex index, if one exists.
    // Note that the current implementation loops over all quarter edges.
    QuarterEdgeIndex GetQuarterEdge(VertexIndex a) const;

    // Get the quarter edge pointing from i to j, if it exists.
    // Note that the current implementation loops over all quarter edges.
    QuarterEdgeIndex GetQuarterEdge(VertexIndex a, VertexIndex b) const;

    // Get the vertex index for the given quarter edge, assuming it is a primal edge.
    VertexIndex GetVertexIndex(QuarterEdgeIndex i) const { return quarter_edges_[i.i].i_vertex; }

    // Get the next vertex or quarter edge index in the list, if it exists.
    // Returns an invalid index if it does not.
    VertexIndex GetNext(VertexIndex i) const;
    QuarterEdgeIndex GetNext(QuarterEdgeIndex i) const;

    // Get the first alive vertex or quarter edge index, if one exists.
    VertexIndex GetFirstVertexIndex() const { return i_vertex_alive_first_; };
    QuarterEdgeIndex GetFirstQuarterEdgeIndex() const { return i_qe_alive_first_; };

    // Get the next right-hand (CCW) QuarterEdge with the same origin.
    QuarterEdgeIndex Next(QuarterEdgeIndex qe) const;

    // Get the next right-hand (CCW) QuarterEdge associated with the same undirected original edge.
    QuarterEdgeIndex Rot(QuarterEdgeIndex qe) const;

    // Get the symmetric QuarterEdge - given QuarterEdge A->B, returns B->A.
    QuarterEdgeIndex Sym(QuarterEdgeIndex qe) const;

    // Get the next left-hand (CW) QuarterEdge associated with the same undirected original edge.
    // This is rot, but the other way.
    QuarterEdgeIndex Tor(QuarterEdgeIndex qe) const;

    // Get the previous right-hand (CCW) QuarterEdge with the same origin.
    // i.e. Get the next left-hand (CW) QuarterEdge with the same origin.
    QuarterEdgeIndex Prev(QuarterEdgeIndex qe) const;

    // Get the next quarter edge that rotates around the same triangle (CCW), in the CCW direction.
    QuarterEdgeIndex Lnext(QuarterEdgeIndex qe) const;

    // Get the quarter edge that shares the same origin as qe_src that has the smallest right-hand
    // (counter clockwise) angle to pointing toward the given point.
    // If there is a quarter edge from src to dst, then this will return that quarter edge.
    QuarterEdgeIndex GetQuarterEdgeRightHandClosestTo(QuarterEdgeIndex qe_src,
                                                      const common::Vec2f& dst) const;

    // A dual edge connects two faces.
    bool IsDual(QuarterEdgeIndex i) const;

    // A primal edge connects two vertices.
    bool IsPrimal(QuarterEdgeIndex i) const;

    // Whether the given quarter edge is constrained.
    bool IsConstrained(QuarterEdgeIndex i) const;

    // Whether the given vertex is one of the boundary vertices.
    // (The first three in vertices_)
    bool IsBoundaryVertex(VertexIndex i_vertex) const;
    bool IsBoundaryVertex(const VertexData& vertex_data) const;
    bool IsBoundaryVertex(QuarterEdgeIndex qe_primal) const;

    bool IsBoundaryEdge(QuarterEdgeIndex qe_primal) const;

    // // Returns true if an edge from i to j exists.
    // // Note that the current implementation loops over all quarter edges.
    // bool HasEdge(int i, int j) const;

    // Set the given edge to be constrained. Note that qe can be primal or dual.
    // This prevents this edge from being flipped in the future.
    void ConstrainEdge(QuarterEdgeIndex i);
    void UnconstrainEdge(QuarterEdgeIndex i);

    // Insert a new vertex into our mesh, returning an index to the new vertex.
    // There are several cases:
    //    1. The new vertex is inside an existing face.
    //       We insert 3 new edges to split the face.
    //       Return the index of the new vertex and an index of any quarter edge with this vertex as
    //       its source.
    //    2. The new vertex lies on an existing non-boundary edge.
    //       We split the edge and add 2 new edges to split the face on either end.
    //       Return the index of the new vertex and an index of a quarter edge with this vertex as
    //       its source that points along the split edge.
    //    3. The new vertex is coincident with an existing vertex.
    //       Return the index of the vertex we are coincident with.
    //    4. The new vertex lies outside the bounding circle (which lies inside the bounding
    //       triangle). Return invalid indices.
    enum class InsertVertexResultCategory { IN_FACE, ON_EDGE, COINCIDENT, OUT_OF_BOUNDS };
    struct InsertVertexResult {
        VertexIndex i_vertex;
        QuarterEdgeIndex i_qe;
        InsertVertexResultCategory category;
    };
    InsertVertexResult InsertVertex(const common::Vec2f& p);

    // Update the mesh to continue to be a Delaunay triangularization, subject to edge constraints.
    // We pass in a primal quarter edge (originating at the vertex to walk around) at which to start
    // the enforcement.
    void EnforceLocallyDelaunay(QuarterEdgeIndex qe_start);

    // Update the mesh such that there is an edge between vertices A and B, where A and B are given
    // by primal quarter edges. The method returns the quarter edge from a to b if the modification
    // was successful, and the invalid quarter edge otherwise.
    // The current implentation may partially modify the mesh even if unsuccessful.
    // This method will not succeed if there is a constrained edge that overlaps AB.
    QuarterEdgeIndex EnforceEdge(QuarterEdgeIndex qe_a, QuarterEdgeIndex qe_b);

    // // Find the triangle in our mesh that encloses the given point.
    // // Return a dual quarter edge originating from the triangle face.
    QuarterEdgeIndex GetEnclosingTriangle(const common::Vec2f& p, QuarterEdgeIndex qe_dual) const;
    QuarterEdgeIndex GetEnclosingTriangle(const common::Vec2f& p) const;

    // Given a dual quarter edge, return the quarter edge originating at the given vertex on that
    // face, pointed in a right-hand orientation.
    std::tuple<QuarterEdgeIndex, QuarterEdgeIndex, QuarterEdgeIndex> GetTriangleQuarterEdges(
        QuarterEdgeIndex qe_dual) const;

    // Move the given vertex (as a primal qe) toward pos. The vertex cannot be moved outside of its
    // surrounding polygon, nor can it be moved closer than min_dist_to_vertex_ to any vertex or
    // min_dist_to_edge_ to any edge.
    void MoveVertexToward(QuarterEdgeIndex qe_primal, const common::Vec2f& pos);

    // Flip the given edge if we can. Returns whether a flip occurred.
    bool MaybeFlipEdge(QuarterEdgeIndex qe_primal);

  private:
    // Private mutable getters
    VertexData& Get(VertexIndex i) { return vertices_[i.i]; }
    QuarterEdge& Get(QuarterEdgeIndex i) { return quarter_edges_[i.i]; }
    const VertexData& Get(VertexIndex i) const { return vertices_[i.i]; }
    const QuarterEdge& Get(QuarterEdgeIndex i) const { return quarter_edges_[i.i]; }

    // Get the index of the next free vertex or quarter edge, retrieving it from the alive list
    // if there are any available, and otherwise appending a new QuarterEdge to the list.
    VertexIndex LivenVertex();
    QuarterEdgeIndex LivenQuarterEdge();

    // Add a new vertex to the mesh, returning the index of the newly added vertex.
    VertexIndex AddVertex(float x, float y);

    // Add a new edge between the vertices at index a and b.
    // Create the quarter edges associated with the given undirected edge.
    // We always create quarter-edges in groups of four.
    // Returns the index of the quarter-edge from A to B.
    QuarterEdgeIndex AddEdge(VertexIndex a, VertexIndex b);

    // Move the given vertex or quarter edge (which must be in the alive list) into the free list.
    // We assume it is no longer referenced by any edges.
    void FreeVertex(VertexIndex i);
    void FreeQuarterEdge(QuarterEdgeIndex i);

    // A utility function used by Splice.
    void SwapNexts(QuarterEdgeIndex a, QuarterEdgeIndex b);

    // A utility function used to join quarter edges
    void Splice(QuarterEdgeIndex a, QuarterEdgeIndex b);

    // A utility method for flipping an edge within a quadrilateral.
    // Should only be called if the bounding quad is convex.
    void FlipEdgeImpl(QuarterEdgeIndex qe);

    struct EnforceEdgeInternalResult {
        bool progress;  // whether progress was made
        QuarterEdgeIndex qe_ab;
    };
    EnforceEdgeInternalResult EnforceEdgeInternal(QuarterEdgeIndex qe_a, VertexIndex i_vertex_b);

    // The maximum radius that a point can be from the origin
    float bounding_radius_;

    // Points within this distance are considered colinear
    float min_dist_to_vertex_;

    // A point this close to an edge is considered coincident
    float min_dist_to_edge_;

    // All of the stored vertices.
    // The vector only grows - it never shrinks.
    std::vector<VertexData> vertices_ = {};
    VertexIndex i_vertex_alive_first_ = {kInvalidIndex};  // Index of the first alive vertex, if any
    VertexIndex i_vertex_alive_last_ = {kInvalidIndex};   // Index of the last alive vertex, if any
    VertexIndex i_vertex_free_first_ = {kInvalidIndex};   // Index of the first free vertex, if any
    VertexIndex i_vertex_free_last_ = {kInvalidIndex};    // Index of the last free vertex, if any
    usize n_vertices_ = 0;

    // All of the stored quarter edges.
    // The vector only grows - it never shrinks.
    std::vector<QuarterEdge> quarter_edges_ = {};
    QuarterEdgeIndex i_qe_alive_first_ = {kInvalidIndex};  // Index of the first alive qe, if any
    QuarterEdgeIndex i_qe_alive_last_ = {kInvalidIndex};   // Index of the last alive qe, if any
    QuarterEdgeIndex i_qe_free_first_ = {kInvalidIndex};   // Index of the first free qe, if any
    QuarterEdgeIndex i_qe_free_last_ = {kInvalidIndex};    // Index of the last free qe, if any
    usize n_quarter_edges_ = 0;
};

}  // namespace core