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
const std::string kAssetEntrySectors = "sectors";

struct TextureInfo {
    u16 texture_id;  // Index in the texture atlas. Id 0 indicates not to draw anything.
    i16 x_offset;    // Texture x offset
    i16 y_offset;    // Texture y offset
    // TODO: scale
};

constexpr u16 kSideInfoFlag_PASSABLE = 2;  // The side info can be traversed

// The information associated with one side of an edge between vertices in the map.
// If this represents the directed edge A->B, then it describes the edge viewed on the right
// side of A->B.

// A solid sideinfo does not set its lower or upper texture, and thus only has the middle
// texture.

// A passable side info's middle texture is not displayed.

struct SideInfo {
    u16 flags;
    u16 sector_id;                    // Index into the sector list
    TextureInfo texture_info_lower;   // Texture displayed if the floor height increases
    TextureInfo texture_info_middle;  // Texture displayed if the wall is solid
    TextureInfo texture_info_upper;   // Texture displayed if the floor ceiling decreases
    QuarterEdgeIndex qe;              // The primary quarter edge that represents A->B

    QuarterEdgeIndex qe_portal;       // The sideinfo for the other side of a portal
};

// All side infos refer to a sector to contain the corresponding floor/ceiling information.
struct Sector {
    u16 flags;
    f32 z_floor;
    f32 z_ceil;
};

// Represents our game map
class GameMap {
  public:
    GameMap();
    ~GameMap() = default;

    // Empty the game map
    void Clear();

    u16 GetMaxSectorIndex() const { return (u16)(sectors_.size() - 1); }

    const std::map<QuarterEdgeIndex, SideInfo>& GetSideInfos() const { return side_infos_; };

    // Returns a pointer to the corresponding side info, or nullptr if it does not exist.
    const SideInfo* GetSideInfo(QuarterEdgeIndex qe_primal) const;
    SideInfo* GetEditableSideInfo(QuarterEdgeIndex qe_primal);

    // Returns a pointer to the corresponding sector, or nullptr if it does not exist.
    const Sector* GetSector(u16 sector_index) const;
    Sector* GetEditableSector(u16 sector_index);

    const core::DelaunayMesh& GetMesh() const { return mesh_; }

    // Insert a new vertex into the mesh.
    VertexIndex AddVertex(const common::Vec2f& pos);

    // Insert a sideinfo for a directed edge into the game map, returning whether it was successful.
    bool AddSideInfo(QuarterEdgeIndex qe_primal);

    // // Remove a vertex (side_info). Can only be called if there are no ajoining sides.
    // // This action invalidates the mesh.
    // bool RemoveVertex(usize edge_index);

    // // Remove a directed edge (side_info). This action does not invalidate the mesh.
    // bool RemoveDirectedEdge(usize edge_index);

    // Have the GameMap create a new sector, and return its index.
    u16 AddSector();

    // Move the vertex given by the primal quarter edge toward pos, within its surrounding polygon.
    void MoveVertexToward(QuarterEdgeIndex qe_primal, const common::Vec2f& pos);

    // Flip the given edge if we can. Returns whether a flip occurred.
    bool MaybeFlipEdge(QuarterEdgeIndex qe_primal);

    bool MaybeEnforceEdge(QuarterEdgeIndex qe_a, QuarterEdgeIndex qe_b) {
        return mesh_.EnforceEdge(qe_a, qe_b);
    }

    // // Write the GameMap entries into the exporter.
    bool Export(core::AssetsExporter* exporter) const;

    // Load the GameMap from the given file. This results the GameMap.
    bool Import(const core::AssetsExporter& exporter);

    // Load a game map from DOOM data.
    bool LoadFromDoomData(const u8* vertex_data, u32 vertex_data_size, const u8* segs_data,
                          u32 segs_data_size, const u8* subsectors_data, u32 subsectors_data_size,
                          const u8* linedefs_data, u32 linedefs_data_size);

    // Finds a vertex (as a primal quarter edge) near the given position, if there is one, using the
    // provided dual quarter edge representing the face containing pos to more quickly locate
    // the vertex.
    QuarterEdgeIndex FindVertexNearPosition(const common::Vec2f& pos, QuarterEdgeIndex qe_dual,
                                            f32 tolerance = 0.3) const;

    // Finds an edge near the given position, if there is one.
    QuarterEdgeIndex FindEdgeNearPosition(const common::Vec2f& pos, QuarterEdgeIndex qe_dual,
                                          f32 tolerance = 0.2) const;

  private:
    AssetsExporterEntry ExportDelaunayMesh(const std::string& name) const;
    AssetsExporterEntry ExportSideInfos(const std::string& name) const;
    AssetsExporterEntry ExportSectors(const std::string& name) const;

    bool LoadDelaunayMesh(core::DelaunayMesh* mesh, const std::string& name,
                          const core::AssetsExporter& exporter);
    bool LoadSideInfos(const std::string& name, const core::AssetsExporter& exporter);
    bool LoadSectors(const std::string& name, const core::AssetsExporter& exporter);

    // The map geometry
    core::DelaunayMesh mesh_;

    // All of the map-related side information.
    // Only primal quarter edges should ever be associated with side_infos.
    // Any edges that do not have side infos are simply transparent.
    std::map<QuarterEdgeIndex, SideInfo> side_infos_;

    // Provides height information for a set of faces.
    std::vector<Sector> sectors_;
};

}  // namespace core