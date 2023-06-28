#pragma once

#include <cstddef>
#include <limits>
#include <vector>

#include "geometry_utils.hpp"
#include "typedefs.hpp"

namespace mesh {

constexpr usize kInvalidIndex = std::numeric_limits<usize>::max();

struct VertexData {
    common::Vec2f v;  // Typically set to (NaN,NaN) if free
    usize i_self;     // Index of itself
    usize i_next;     // Index of the next alive/free VertexData (typemax(size_t) if none)
    usize i_prev;     // Index of the previous alive/free VertexData (typemax(size_t) if none)
};

// A quarter edge represents a directed edge in a QuadEdgeTree.
// Each undirected edge in the graph A <-> B with face L on the left of A->B and face R on the right
// of A->B has four quarter edges:
//    A->B, L->R, B->A, and R->L.
// If this quarter edge is primal, then the vertex index points to its vertex.
// If this quarter edge is dual (i.e., represents a face), then ivertex is set to typemax(size_t).
struct QuarterEdge {
    usize i_vertex;   // Non-null if this is a quarter edge originating at a vertex
                      // Set to typemax(size_t) if free or dual
    usize i_qe_next;  // The next right-hand (CCW) QuarterEdge with the same origin
    usize i_qe_rot;   // The next right-hand (CCW) QuarterEdge associated with the same
                      // undirected original edge.
    usize i_self;     // Index of itself
    usize i_next;     // Index of the next alive/free QuarterEdge (typemax(size_t) if none)
    usize i_prev;     // Index of the previous alive/free QuarterEdge (typemax(size_t) if none)
};

// A dual edge connects two faces.
bool IsDualEdge(const QuarterEdge& qe);

// A primal edge connects two vertices.
bool IsPrimalEdge(const QuarterEdge& qe);

// A DelaunayMesh is a Quad Edge Mesh that can dynamically generate a Delaunay (or constrained
// Delaunay) triangulation. A Quad Edge Mesh's name comes from the fact that each edge in a "normal"
// mesh graph that connects two vertices is represented by 4 edges in the quad edge tree - one
// directed edge between each vertex and one directed edge between the two adjoined faces. Because
// of this, we always have a multiple of 4 quater edges.
//
// This class allows for the insertion and removal of items. As such, our vectors may change in
// length as memory reallocates. We use integers to index into our vertices and quarter edges to
// avoid memeory reallocation issues.
class DelaunayMesh {
  public:
    // Construct a QuadEdgeMesh that is merely a bounding triangle with the given radius.
    DelaunayMesh(float bounding_radius, float min_dist_to_vertex = 1e-3,
                 float min_dist_to_edge = 1e-3);
    ~DelaunayMesh();

    // size_t NumVertices() const { return vertices_.size(); }
    // size_t NumQuarterEdges() const { return quarter_edges_.size(); }
    // size_t NumEdges() const { return NumQuarterEdges() / 4; }

    // void Clear();
    // bool LoadFromData(const u8* data);

    // const common::Vec2f& GetVertex(int i) const { return vertices_.at(i)->vertex; }
    // QuarterEdge* GetQuarterEdge(int i) const { return quarter_edges_.at(i); }

    // // Get the quarter edge pointing from i to j, if it exists.
    // // Note that the current implementation loops over all quarter edges.
    // QuarterEdge* GetQuarterEdge(int i, int j) const;

    // // Returns true if an edge from i to j exists.
    // // Note that the current implementation loops over all quarter edges.
    // bool HasEdge(int i, int j) const;

    // // Insert a new vertex into our mesh, and update the mesh to continue to be a Delaunay
    // // triangularization. Returns the index of the added vertex if a new point was added, and
    // // kInvalidIndex instead. (If the new point is coincident with an existing point or the
    // boundary
    // // edge, no new point is added.)
    // int AddDelaunayVertex(const common::Vec2f& p);

    // // Enforce an edge between the ith and jth vertices.
    // // Returns true if this operation was successful.
    // bool ConstrainEdge(int i, int j);

    // // Get the next right-hand (CCW) QuarterEdge with the same origin.
    // QuarterEdge* Next(const QuarterEdge* qe) const;

    // // Get the next right-hand (CCW) QuarterEdge associated with the same undirected original
    // edge. QuarterEdge* Rot(const QuarterEdge* qe) const;

    // // Get the symmetric QuarterEdge - given QuarterEdge A->B, returns B->A.
    // QuarterEdge* Sym(const QuarterEdge* qe) const;

    // // Get the next left-hand (CW) QuarterEdge associated with the same undirected original edge.
    // // This is rot, but the other way.
    // QuarterEdge* Tor(const QuarterEdge* qe) const;

    // // Get the previous right-hand (CCW) QuarterEdge with the same origin.
    // // i.e. Get the next left-hand (CW) QuarterEdge with the same origin.
    // QuarterEdge* Prev(const QuarterEdge* qe) const;

    // // Get the next quarter edge that rotates around the same triangle (CCW), in the CCW
    // direction. QuarterEdge* Lnext(const QuarterEdge* qe) const;

    // // Find the triangle in our mesh that encloses the given point.
    // // Return a dual quarter edge originating from the triangle face.
    // QuarterEdge* GetEnclosingTriangle(const common::Vec2f& p, QuarterEdge* qe_ab) const;
    // QuarterEdge* GetEnclosingTriangle(const common::Vec2f& p) const;

    // // Given a dual quarter edge, return the quarter edge originating at the given vertex on that
    // // face, pointed in a right-hand orientation.
    // QuarterEdge* GetTriangleQuarterEdge1(const QuarterEdge* qe_dual) const;
    // QuarterEdge* GetTriangleQuarterEdge2(const QuarterEdge* qe_dual) const;
    // QuarterEdge* GetTriangleQuarterEdge3(const QuarterEdge* qe_dual) const;

    // // Given a dual quarter edge, return the first vertex on that face.
    // const common::Vec2f& GetTriangleVertex1(const QuarterEdge* qe_dual) const;
    // const common::Vec2f& GetTriangleVertex2(const QuarterEdge* qe_dual) const;
    // const common::Vec2f& GetTriangleVertex3(const QuarterEdge* qe_dual) const;

    // // Whether the given vertex is one of the boundary vertices.
    // // (The first three in vertices_)
    // bool IsBoundaryVertex(const common::Vec2f* v) const;
    // bool IsBoundaryVertex(const VertexData& vertex_data) const;
    // bool IsBoundaryVertex(const VertexData* const& vertex_data) const;

  private:
    // Add a new vertex to the mesh, returning the index of the newly added vertex.
    usize AddVertex(float x, float y);

    // Get the index of the next free quarter edge, retrieving it from the alive list
    // if there are any available, and otherwise appending a new QuarterEdge to the list.
    usize PrepareFreeQuarterEdge();

    // Add a new edge between the vertices at index a and b.
    // Create the quarter edges associated with the given undirected edge.
    // We always create quarter-edges in groups of four.
    // Returns the index of the quarter-edge from A to B.
    usize AddEdge(usize i_vertex_a, usize i_vertex_b);

    // // A utility function used to join quarter edges
    // void Splice(QuarterEdge* a, QuarterEdge* b);

    // // A utility method for flipping an edge within a quadrilateral.
    // // Should only be called if the bounding quad is convex.
    // void FlipEdge(QuarterEdge* qe);

    // The maximum radius that a point can be from the origin
    float bounding_radius_;

    // Points within this distance are considered colinear
    float min_dist_to_vertex_;

    // A point this close to an edge is considered coincident
    float min_dist_to_edge_;

    // All of the stored vertices.
    // The vector only grows - it never shrinks.
    std::vector<VertexData> vertices_ = {};
    usize i_vertex_alive_first_ = kInvalidIndex;  // Index of the first alive vertex, if any
    usize i_vertex_alive_last_ = kInvalidIndex;   // Index of the last alive vertex, if any
    usize i_vertex_free_first_ = kInvalidIndex;   // Index of the first free vertex, if any
    usize i_vertex_free_last_ = kInvalidIndex;    // Index of the last free vertex, if any

    // All of the stored quarter edges.
    // The vector only grows - it never shrinks.
    std::vector<QuarterEdge> quarter_edges_ = {};
    usize i_qe_alive_first_ = kInvalidIndex;  // Index of the first alive quarter edge, if any
    usize i_qe_alive_last_ = kInvalidIndex;   // Index of the last alive quarter edge, if any
    usize i_qe_free_first_ = kInvalidIndex;   // Index of the first free quarter edge, if any
    usize i_qe_free_last_ = kInvalidIndex;    // Index of the last free quarter edge, if any
};

}  // namespace mesh