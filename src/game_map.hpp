#pragma once

#include <map>
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

    // The editable map vertices.
    // This will exactly match the vertices in the DelaunayMesh.
    // (We error if adding any vertex to the DelaunayMesh fails (due to coincidence, for example).)
    std::vector<common::Vec2f> vertices;

    // All of the map-related side information.
    // If we have side information for an edge A -> B, then that edge must end up in the mesh.
    // Edges without side information may exist in the mesh. Such edges are assumed transparent.
    std::vector<SideInfo> side_infos;

    // Map <a_ind, b_ind> to index in side_infos.
    std::map<std::tuple<usize, usize>, usize> side_to_info;

    // The map geometry
    core::DelaunayMesh mesh =
        core::DelaunayMesh(MESH_BOUNDING_RADIUS, MESH_MIN_DIST_TO_VERTEX, MESH_MIN_DIST_TO_EDGE);

    void Clear();

    bool HasEdge(int a_ind, int b_ind) const;

    // Insert a directed edge into the game map, returning its index.
    usize AddDirectedEdge(int a_ind, int b_ind);

    // Write the GameMap entries into the exporter.
    bool Export(core::AssetsExporter* exporter) const;

    // Load the GameMap from the given file. This results the GameMap.
    bool Import(const core::AssetsExporter& exporter);

  private:
    AssetsExporterEntry ExportDelaunayMesh(const std::string& name) const;
    AssetsExporterEntry ExportSideInfos(const std::string& name) const;

    bool LoadDelaunayMesh(core::DelaunayMesh* mesh, const std::string& name,
                          const core::AssetsExporter& exporter);
    bool LoadSideInfos(const std::string& name, const core::AssetsExporter& exporter);
};

int MapVertexIndexToMeshVertexIndex(int ind);
int MeshVertexIndexToMapVertexIndex(int ind);

}  // namespace core