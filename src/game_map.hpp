#pragma once

#include <map>
#include <memory>
#include <optional>
#include <vector>

#include "assets_exporter.hpp"
#include "delaunay_mesh.hpp"
#include "typedefs.hpp"

#define MESH_BOUNDING_RADIUS 1000.0f
#define MESH_MIN_DIST_TO_VERTEX 0.1f
#define MESH_MIN_DIST_TO_EDGE 0.1f

namespace core {

const std::string kAssetEntryGeometryMesh = "geometry_mesh";
const std::string kAssetEntrySideInfos = "side_infos";

constexpr u16 kSideInfoFlag_DARK = 1;

// The information associated with one side of an edge between vertices in the map.
// If this represents the directed edge A->B, then it describes the edge viewed on the right
// side of A->B.
struct SideInfo {
    u16 flags;
    u16 texture_id;
    i16 x_offset;  // Texture x offset
    i16 y_offset;  // Texture y offset
    usize a_ind;   // The index of vertex A (the source vertex)
    usize b_ind;   // The index of vertex B (the dest vertex)
};

// Represents our game map
class GameMap {
  public:
    GameMap() = default;
    ~GameMap() = default;

    // Empty the game map
    void Clear();

    void InvalidateMesh() { mesh_.reset(); }

    const std::vector<common::Vec2f>& GetVertices() const { return vertices_; };
    const std::vector<SideInfo>& GetSideInfos() const { return side_infos_; };

    // Whether the given edge exists
    bool HasEdge(int a_ind, int b_ind) const;
    std::optional<usize> GetEdgeIndex(int a_ind, int b_ind) const;

    // Insert a directed edge into the game map, returning its index.
    usize AddDirectedEdge(int a_ind, int b_ind);

    // Remove a directed edge (side_info). This action invalidates the mesh.
    bool RemoveDirectedEdge(usize edge_index);

    // Write the GameMap entries into the exporter.
    bool Export(core::AssetsExporter* exporter) const;

    // Load the GameMap from the given file. This results the GameMap.
    bool Import(const core::AssetsExporter& exporter);

    // Finds a vertex (map index) near the given position, if there is one.
    // If there is a mesh, use the provided dual quarter edge (representing a face) to more quickly
    // locate the vertex. Otherwise, fall back to slow and simple iteration over all vertices.
    // Always returns a map vertex index.
    std::optional<usize> FindVertexNearPosition(const common::Vec2f& pos,
                                                core::QuarterEdge* qe_face,
                                                f32 tolerance = 0.3) const;

    // Finds an edge (map index -> map index) near the given position, if there is one.
    // If there is a mesh, use the provided dual quarter edge (representing a face) to more quickly
    // locate the vertex. Otherwise, fall back to slow and simple iteration over all edges.
    std::optional<std::pair<usize, usize>> FindEdgeNearPosition(const common::Vec2f& pos,
                                                                core::QuarterEdge* qe_face,
                                                                f32 tolerance = 0.2) const;

  private:
    AssetsExporterEntry ExportDelaunayMesh(const std::string& name) const;
    AssetsExporterEntry ExportSideInfos(const std::string& name) const;

    bool LoadDelaunayMesh(core::DelaunayMesh* mesh, const std::string& name,
                          const core::AssetsExporter& exporter);
    bool LoadSideInfos(const std::string& name, const core::AssetsExporter& exporter);

    // FindVertexNearPosition when a mesh is present.
    // The given quarter edge must be in the mesh.
    std::optional<usize> FindVertexNearPosition(core::DelaunayMesh* mesh, const common::Vec2f& pos,
                                                core::QuarterEdge* qe_face, f32 tolerance) const;

    // FindEdgeNearPosition when a mesh is present.
    // The given quarter edge must be in the mesh.
    std::optional<std::pair<usize, usize>> FindEdgeNearPosition(core::DelaunayMesh* mesh,
                                                                const common::Vec2f& pos,
                                                                core::QuarterEdge* qe_face,
                                                                f32 tolerance) const;

    // The editable map vertices.
    // This will exactly match the vertices in the DelaunayMesh.
    // (We error if adding any vertex to the DelaunayMesh fails (due to coincidence, for
    // example).)
    std::vector<common::Vec2f> vertices_;

    // All of the map-related side information.
    // If we have side information for an edge A -> B, then that edge must end up in the mesh.
    // Edges without side information may exist in the mesh. Such edges are assumed transparent.
    std::vector<SideInfo> side_infos_;

    // Map <a_ind, b_ind> to index in side_infos.
    std::map<std::tuple<usize, usize>, usize> side_to_info_;

    // The map geometry
    std::unique_ptr<core::DelaunayMesh> mesh_;
};

int MapVertexIndexToMeshVertexIndex(int ind);
int MeshVertexIndexToMapVertexIndex(int ind);

}  // namespace core