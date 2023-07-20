#pragma once

#include <map>
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

struct TextureInfo {
    u16 texture_id;  // Index in the texture atlas
    i16 x_offset;    // Texture x offset
    i16 y_offset;    // Texture y offset
    // TODO: scale
};

constexpr u16 kSideInfoFlag_DARK = 1;
constexpr u16 kSideInfoFlag_PASSABLE = 2;

// The information associated with one side of an edge between vertices in the map.
// If this represents the directed edge A->B, then it describes the edge viewed on the right
// side of A->B.

// A solid sideinfo does not set its lower or upper texture, and thus only has the middle
// texture, z_floor, and z_ceil.

// A floor height change sideinfo has no middle texture, and thus has a lower and/or upper texture.

struct SideInfo {
    u16 flags;
    TextureInfo texture_info_lower;   // Texture displayed if the floor height increases
    TextureInfo texture_info_middle;  // Texture displayed if the wall is solid
    TextureInfo texture_info_upper;   // Texture displayed if the floor ceiling decreases
    f32 z_floor;  // Height of the floor on the side the edge is viewed from (right of A->B)
    f32 z_lower;  // Height of the boundary between the lower and middle texture
    f32 z_upper;  // Height of the boundary between the middle and upper texture
    f32 z_ceil;   // Height of the ceiling on the side the edge is viewed from (right of A-B)
    QuarterEdgeIndex qe;  // The primary quarter edge that represents A->B
};

// Represents our game map
class GameMap {
  public:
    GameMap();
    ~GameMap() = default;

    // Empty the game map
    void Clear();

    const std::map<QuarterEdgeIndex, SideInfo>& GetSideInfos() const { return side_infos_; };

    // Returns a pointer to the corresponding side info, or nullptr if it does not exist.
    SideInfo* GetEditableSideInfo(QuarterEdgeIndex qe_primal);

    // // Whether the given edge exists
    // bool HasEdge(int a_ind, int b_ind) const;
    // std::optional<usize> GetEdgeIndex(int a_ind, int b_ind) const;

    const core::DelaunayMesh& GetMesh() const { return mesh_; }

    // Insert a new vertex into the mesh.
    VertexIndex AddVertex(const common::Vec2f& pos);

    // // Insert a directed edge into the game map, returning its index.
    // usize AddDirectedEdge(int a_ind, int b_ind);

    // // Remove a vertex (side_info). Can only be called if there are no ajoining sides.
    // // This action invalidates the mesh.
    // bool RemoveVertex(usize edge_index);

    // // Remove a directed edge (side_info). This action does not invalidate the mesh.
    // bool RemoveDirectedEdge(usize edge_index);

    // Have the GameMap create a new FaceInfo for the given face.
    bool AddFaceInfo(QuarterEdgeIndex qe_dual);

    // Move the vertex given by the primal quarter edge toward pos, within its surrounding polygon.
    void MoveVertexToward(QuarterEdgeIndex qe_primal, const common::Vec2f& pos);

    // Flip the given edge if we can. Returns whether a flip occurred.
    bool MaybeFlipEdge(QuarterEdgeIndex qe_primal);

    // // Write the GameMap entries into the exporter.
    bool Export(core::AssetsExporter* exporter) const;

    // Load the GameMap from the given file. This results the GameMap.
    bool Import(const core::AssetsExporter& exporter);

    // Finds a vertex (as a primal quarter edge) near the given position, if there is one, using the
    // provided dual quarter edge representing the face containing pos to more quickly locate the
    // vertex.
    QuarterEdgeIndex FindVertexNearPosition(const common::Vec2f& pos, QuarterEdgeIndex qe_dual,
                                            f32 tolerance = 0.3) const;

    // Finds an edge near the given position, if there is one.
    QuarterEdgeIndex FindEdgeNearPosition(const common::Vec2f& pos, QuarterEdgeIndex qe_dual,
                                          f32 tolerance = 0.2) const;

  private:
    AssetsExporterEntry ExportDelaunayMesh(const std::string& name) const;
    AssetsExporterEntry ExportSideInfos(const std::string& name) const;

    bool LoadDelaunayMesh(core::DelaunayMesh* mesh, const std::string& name,
                          const core::AssetsExporter& exporter);
    bool LoadSideInfos(const std::string& name, const core::AssetsExporter& exporter);

    // The map geometry
    core::DelaunayMesh mesh_;

    // All of the map-related side information.
    // Only primal quarter edges should ever be associated with side_infos.
    // Any edges that do not have side infos are simply transparent.
    std::map<QuarterEdgeIndex, SideInfo> side_infos_;
};

}  // namespace core