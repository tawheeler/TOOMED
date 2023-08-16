#include "game_map.hpp"

#include <cstring>
#include <iostream>
#include <sstream>

#include "math_utils.hpp"

namespace core {

struct ExportedSideInfo {
    u16 flags;
    u16 sector_id;                    // Index into the sector list
    TextureInfo texture_info_lower;   // Texture displayed if the floor height increases
    TextureInfo texture_info_middle;  // Texture displayed if the wall is solid
    TextureInfo texture_info_upper;   // Texture displayed if the floor ceiling decreases
};

// ------------------------------------------------------------------------------------------------
GameMap::GameMap() : mesh_(MESH_BOUNDING_RADIUS, MESH_MIN_DIST_TO_VERTEX, MESH_MIN_DIST_TO_EDGE) {
    // Start off with a default sector (id 0)
    Sector sector = {};
    sector.z_floor = 0.0;
    sector.z_ceil = 1.0;
    sectors_.push_back(sector);
}

// ------------------------------------------------------------------------------------------------
void GameMap::Clear() {
    mesh_.Clear();
    side_infos_.clear();
}

// ------------------------------------------------------------------------------------------------
const SideInfo* GameMap::GetSideInfo(QuarterEdgeIndex qe_primal) const {
    auto it = side_infos_.find(qe_primal);
    if (it == side_infos_.end()) {
        return nullptr;
    }
    return &(it->second);
}

// ------------------------------------------------------------------------------------------------
SideInfo* GameMap::GetEditableSideInfo(QuarterEdgeIndex qe_primal) {
    auto it = side_infos_.find(qe_primal);
    if (it == side_infos_.end()) {
        return nullptr;
    }
    return &(it->second);
}

// ------------------------------------------------------------------------------------------------
const Sector* GameMap::GetSector(u16 sector_index) const {
    if (sectors_.size() <= sector_index) {
        return nullptr;
    }
    return &(sectors_[sector_index]);
}

// ------------------------------------------------------------------------------------------------
Sector* GameMap::GetEditableSector(u16 sector_index) {
    if (sectors_.size() <= sector_index) {
        return nullptr;
    }
    return &(sectors_[sector_index]);
}

// ------------------------------------------------------------------------------------------------
VertexIndex GameMap::AddVertex(const common::Vec2f& pos) {
    DelaunayMesh::InsertVertexResult result = mesh_.InsertVertex(pos);

    // TODO: Handle side_info data based on how the vertex was inserted.

    if (IsValid(result.i_vertex)) {
        mesh_.EnforceLocallyDelaunay(result.i_qe);
    }

    return result.i_vertex;
}

// ------------------------------------------------------------------------------------------------
bool GameMap::AddSideInfo(QuarterEdgeIndex qe_primal) {
    if (!mesh_.IsPrimal(qe_primal)) {
        return false;  // Only add side infos for primal quarter edges
    }
    if (side_infos_.find(qe_primal) != side_infos_.end()) {
        return false;  // There already is one
    }

    // Ensure that the side is constrained.
    mesh_.ConstrainEdge(qe_primal);

    // Add a new sideinfo.
    SideInfo side_info = {};
    side_info.qe = qe_primal;
    side_infos_[side_info.qe] = side_info;

    return true;
}

// ------------------------------------------------------------------------------------------------
u16 GameMap::AddSector() {
    u16 idx = sectors_.size();
    Sector sector = {};
    sector.z_floor = 0.0;
    sector.z_ceil = 1.0;
    sectors_.push_back(sector);
    return idx;
}

// //
// ------------------------------------------------------------------------------------------------
// bool GameMap::RemoveVertex(usize vertex_index) {
//     if (vertex_index >= vertices_.size()) {
//         return false;
//     }

//     // As a sanity check, verify that no sideinfos share this edge.
//     // @efficiency This is expensive.
//     for (const auto& side_info : side_infos_) {
//         if (side_info.a_ind == vertex_index || side_info.b_ind == vertex_index) {
//             return false;
//         }
//     }

//     // Remove the vertex
//     vertices_.erase(vertices_.begin() + vertex_index);

//     return true;
// }

// //
// ------------------------------------------------------------------------------------------------
// bool GameMap::RemoveDirectedEdge(usize edge_index) {
//     if (edge_index >= side_infos_.size()) {
//         return false;
//     }

//     // Remove the edge from side_to_info
//     const SideInfo& side_info = side_infos_[edge_index];
//     side_to_info_.erase(std::make_pair(side_info.a_ind, side_info.b_ind));

//     // TODO: Update side-to-info because we need to update all of the indices for things higher
//     than
//     // what we removed. AAAAARG.
//     //       Maybe just have an ordered map of side infos instead of a list plus map?
//     //       I also need each vertex to have a uid, so that if I delete a vertex I don't have to
//     //       update all side_infos.
//     //       ????

//     // Remove the item for side_infos
//     side_infos_.erase(side_infos_.begin() + edge_index);

//     return true;
// }

// ------------------------------------------------------------------------------------------------
void GameMap::MoveVertexToward(QuarterEdgeIndex qe_primal, const common::Vec2f& pos) {
    mesh_.MoveVertexToward(qe_primal, pos);
}

// ------------------------------------------------------------------------------------------------
bool GameMap::MaybeFlipEdge(QuarterEdgeIndex qe_primal) { return mesh_.MaybeFlipEdge(qe_primal); }

// ------------------------------------------------------------------------------------------------
AssetsExporterEntry GameMap::ExportDelaunayMesh(const std::string& name) const {
    AssetsExporterEntry entry;
    entry.SetName(name);

    // Serialize the mesh
    std::stringstream buffer;

    // Serialize the counts
    u32 n_vertices = mesh_.NumVertices();
    buffer.write(reinterpret_cast<const char*>(&n_vertices), sizeof(u32));
    u32 n_quarter_edges = mesh_.NumQuarterEdges();
    buffer.write(reinterpret_cast<const char*>(&n_quarter_edges), sizeof(u32));

    // Serialize the vertices

    std::map<VertexIndex, u32> vertex_2_serialization_index;  // Remember the order

    VertexIndex i_vertex = mesh_.GetFirstVertexIndex();
    while (IsValid(i_vertex)) {
        const common::Vec2f& v = mesh_.GetVertex(i_vertex);
        buffer.write(reinterpret_cast<const char*>(&v.x), sizeof(f32));
        buffer.write(reinterpret_cast<const char*>(&v.y), sizeof(f32));
        vertex_2_serialization_index[i_vertex] = (u32)(size(vertex_2_serialization_index));
        i_vertex = mesh_.GetNext(i_vertex);
    }

    // Figure out the quarter edge serialization order
    std::map<QuarterEdgeIndex, u32> qe_2_serialization_index;  // Remember the order
    QuarterEdgeIndex qe = mesh_.GetFirstQuarterEdgeIndex();
    while (IsValid(qe)) {
        qe_2_serialization_index[qe] = (u32)(size(qe_2_serialization_index));
        qe = mesh_.GetNext(qe);
    }

    // Serialize the quarter edges
    qe = mesh_.GetFirstQuarterEdgeIndex();
    while (IsValid(qe)) {
        u32 vertex_index = std::numeric_limits<u32>::max();
        const QuarterEdge& qe_data = mesh_.GetQuarterEdge(qe);
        if (IsValid(qe_data.i_vertex)) {
            vertex_index = vertex_2_serialization_index[qe_data.i_vertex];
        }

        u32 nxt = qe_2_serialization_index[qe_data.i_nxt];
        u32 rot = qe_2_serialization_index[qe_data.i_rot];

        buffer.write(reinterpret_cast<const char*>(&vertex_index), sizeof(u32));
        buffer.write(reinterpret_cast<const char*>(&nxt), sizeof(u32));
        buffer.write(reinterpret_cast<const char*>(&rot), sizeof(u32));

        qe = mesh_.GetNext(qe);
    }

    entry.SetData(buffer.str());

    return entry;
}

// ------------------------------------------------------------------------------------------------
AssetsExporterEntry GameMap::ExportSideInfos(const std::string& name) const {
    AssetsExporterEntry entry;
    entry.SetName(name);

    // --------------------------------------
    // Serialize the content
    std::stringstream buffer;

    // Serialize the counts
    u32 n_side_infos = side_infos_.size();
    buffer.write(reinterpret_cast<const char*>(&n_side_infos), sizeof(u32));

    // Serialize the side infos

    std::map<QuarterEdgeIndex, u16> sideinfo_2_serialization_index;  // Remember the order

    for (const auto& it : side_infos_) {
        const SideInfo& side_info = it.second;

        ExportedSideInfo exported_side_info = {.flags = side_info.flags,
                                               .sector_id = side_info.sector_id,
                                               .texture_info_lower = side_info.texture_info_lower,
                                               .texture_info_middle = side_info.texture_info_middle,
                                               .texture_info_upper = side_info.texture_info_upper};

        buffer.write(reinterpret_cast<const char*>(&exported_side_info),
                     sizeof(exported_side_info));

        sideinfo_2_serialization_index[side_info.qe] = (u16)(size(sideinfo_2_serialization_index));
    }

    // Serialize a vector[u16] of quarter edge index -> side info index.
    // Uses 0xFFFF to indicate no side info.
    QuarterEdgeIndex qe = mesh_.GetFirstQuarterEdgeIndex();
    while (IsValid(qe)) {
        u16 side_info_index = std::numeric_limits<u16>::max();
        auto it_side = sideinfo_2_serialization_index.find(qe);
        if (it_side != sideinfo_2_serialization_index.end()) {
            side_info_index = it_side->second;
        }

        buffer.write(reinterpret_cast<const char*>(&side_info_index), sizeof(u16));

        qe = mesh_.GetNext(qe);
    }

    // --------------------------------------
    entry.SetData(buffer.str());

    return entry;
}

// ------------------------------------------------------------------------------------------------
AssetsExporterEntry GameMap::ExportSectors(const std::string& name) const {
    AssetsExporterEntry entry;
    entry.SetName(name);

    // --------------------------------------
    // Serialize the content
    std::stringstream buffer;

    // Serialize the counts
    u32 n_sectors = sectors_.size();
    buffer.write(reinterpret_cast<const char*>(&n_sectors), sizeof(u32));

    // Serialize the sectors
    for (const Sector& sector : sectors_) {
        buffer.write(reinterpret_cast<const char*>(&sector), sizeof(sector));
    }

    // --------------------------------------
    entry.SetData(buffer.str());

    return entry;
}

// ------------------------------------------------------------------------------------------------
bool GameMap::Export(AssetsExporter* exporter) const {
    exporter->AddEntry(ExportDelaunayMesh(kAssetEntryGeometryMesh));
    exporter->AddEntry(ExportSideInfos(kAssetEntrySideInfos));
    exporter->AddEntry(ExportSectors(kAssetEntrySectors));
    return true;
}

// ------------------------------------------------------------------------------------------------
bool GameMap::LoadDelaunayMesh(DelaunayMesh* mesh, const std::string& name,
                               const AssetsExporter& exporter) {
    const AssetsExporterEntry* entry = exporter.FindEntry(name);
    if (entry == nullptr) {
        std::cout << "Failed to find entry " << name << std::endl;
        return false;
    }

    const u8* data = entry->data.data();
    return mesh->LoadFromData(data);
}

// ------------------------------------------------------------------------------------------------
bool GameMap::LoadSideInfos(const std::string& name, const AssetsExporter& exporter) {
    const AssetsExporterEntry* entry = exporter.FindEntry(name);
    if (entry == nullptr) {
        std::cout << "Failed to find entry " << name << std::endl;
        return false;
    }

    const u8* data = entry->data.data();
    u32 offset = 0;

    u32 n_side_infos = *(u32*)(data + offset);
    offset += sizeof(u32);

    std::vector<SideInfo> side_infos;
    side_infos.reserve(n_side_infos);
    for (u32 i = 0; i < n_side_infos; i++) {
        ExportedSideInfo exported_side_info = *(ExportedSideInfo*)(data + offset);
        offset += sizeof(ExportedSideInfo);

        SideInfo side_info = {};
        side_info.flags = exported_side_info.flags;
        side_info.sector_id = exported_side_info.sector_id;
        side_info.texture_info_lower = exported_side_info.texture_info_lower;
        side_info.texture_info_middle = exported_side_info.texture_info_middle;
        side_info.texture_info_upper = exported_side_info.texture_info_upper;
        side_info.qe_portal = {kInvalidIndex};
        // NOTE: quarter edge index not yet set.

        side_infos.emplace_back(side_info);
    }

    // Set source and dest from side to info list and store in map.
    side_infos_.clear();
    for (size_t i = 0; i < mesh_.NumQuarterEdges(); i++) {
        u16 side_info_index = *(u16*)(data + offset);
        offset += sizeof(u16);
        if (side_info_index != std::numeric_limits<u16>::max()) {
            SideInfo& side_info = side_infos[side_info_index];
            side_info.qe = {i};
            side_infos_[side_info.qe] = side_info;

            // Ensure that the quarter edges associated with any side info are solid.
            mesh_.ConstrainEdge(side_info.qe);
        }

        // u16 portal_side_info_index = *(u16*)(data + offset);
        // offset += sizeof(u16);
        // if (side_info_index != std::numeric_limits<u16>::max()) {
        //     SideInfo& side_info = side_infos[side_info_index];
        //     side_info.qe = {i};
        //     side_infos_[side_info.qe] = side_info;

        //     // Ensure that the quarter edges associated with any side info are solid.
        //     mesh_.ConstrainEdge(side_info.qe);
        // }
    }

    return true;
}

// ------------------------------------------------------------------------------------------------
bool GameMap::LoadSectors(const std::string& name, const AssetsExporter& exporter) {
    const AssetsExporterEntry* entry = exporter.FindEntry(name);
    if (entry == nullptr) {
        std::cout << "Failed to find entry " << name << std::endl;
        return false;
    }

    const u8* data = entry->data.data();
    u32 offset = 0;

    u32 n_sectors = *(u32*)(data + offset);
    offset += sizeof(u32);

    sectors_.clear();
    sectors_.reserve(n_sectors);
    for (u32 i = 0; i < n_sectors; i++) {
        sectors_.push_back(*(Sector*)(data + offset));
        offset += sizeof(Sector);
    }

    return true;
}

// ------------------------------------------------------------------------------------------------
bool GameMap::Import(const AssetsExporter& exporter) {
    // Clear the game map
    Clear();

    // First, load the geometry mesh.
    bool success = true;
    if (exporter.HasEntry(kAssetEntryGeometryMesh)) {
        success &= LoadDelaunayMesh(&mesh_, kAssetEntryGeometryMesh, exporter);
    } else {
        success = false;
        std::cout << "Failed to load geometry mesh!" << std::endl;
    }

    if (success) {
        success &= LoadSideInfos(kAssetEntrySideInfos, exporter);
    }
    if (success) {
        success &= LoadSectors(kAssetEntrySectors, exporter);
    }

    return success;
}

// ------------------------------------------------------------------------------------------------
bool GameMap::LoadFromDoomData(const u8* linedefs_data, u32 linedefs_data_size,
                               const u8* sidedefs_data, u32 sidedefs_data_size,
                               const u8* vertex_data, u32 vertex_data_size, const u8* sectors_data,
                               u32 sectors_data_size, const core::RenderAssets& render_assets) {
    // Empty the game map.
    Clear();

    // Reset the mesh with one that is guaranteed to be large enough to hold the Doom level.
    constexpr f32 kDoomUnitsPerWorldUnit = 64.0f;  // DOOM units are per-pixel.
    constexpr f32 kMaxDoomX = std::numeric_limits<i16>::max() / kDoomUnitsPerWorldUnit;
    mesh_ = DelaunayMesh(kMaxDoomX * 1.41 + 50.0, MESH_MIN_DIST_TO_VERTEX, MESH_MIN_DIST_TO_EDGE);

    // Load the vertices into a temp vector.
    u32 sizeof_vertex = 4;
    u32 n_vertices = vertex_data_size / sizeof_vertex;
    std::vector<common::Vec2f> doom_vertices;
    doom_vertices.reserve(n_vertices);

    u32 vertices_offset = 0;
    for (u32 i = 0; i < n_vertices; i++) {
        i16 x = *(i16*)(vertex_data + vertices_offset);
        vertices_offset += sizeof(i16);
        i16 y = *(i16*)(vertex_data + vertices_offset);
        vertices_offset += sizeof(i16);

        doom_vertices.emplace_back(x / kDoomUnitsPerWorldUnit, y / kDoomUnitsPerWorldUnit);
    }

    // Sidedefs
    struct Sidedef {
        i16 offset_x;
        i16 offset_y;
        u8 texture_name_upper[8];
        u8 texture_name_lower[8];
        u8 texture_name_middle[8];
        i16 sector_index;
    };

    // Keep track of the vertex indices assigned to each doom vertex.
    // Unlike quarter edge indices, these are not going to change as we add more lines.
    std::vector<VertexIndex> vertex_indices(n_vertices, {kInvalidIndex});

    // Add the linedefs one at a time.
    // Introduce the vertices. If there already is a vertex, that is okay.
    // Ensure that we flip edges such that the new linedefs have corresponding edges.
    // Constrain those edges.

    // Linedefs
    struct Linedef {
        i16 i_vertex_start;
        i16 i_vertex_end;
        i16 flags;
        i16 special_type;
        i16 sector_tag;
        i16 i_sidedef_front;
        i16 i_sidedef_back;
    };
    u32 n_linedefs = linedefs_data_size / sizeof(Linedef);

    u32 linedefs_offset = 0;
    for (u32 i_linedef = 0; i_linedef < n_linedefs; i_linedef++) {
        Linedef linedef = *(Linedef*)(linedefs_data + linedefs_offset);
        linedefs_offset += sizeof(Linedef);

        // Insert vertex a and b.
        // The vertex should either be in a face or coincident with an existing point.
        // TODO: It may be possible to be on an edge, but the edge should not be constrained.
        for (i16 i_vertex : {linedef.i_vertex_start, linedef.i_vertex_end}) {
            const common::Vec2f& v = doom_vertices[i_vertex];

            DelaunayMesh::InsertVertexResult res = mesh_.InsertVertex(v);
            if (res.category == DelaunayMesh::InsertVertexResultCategory::OUT_OF_BOUNDS) {
                return false;
            } else if (res.category == DelaunayMesh::InsertVertexResultCategory::IN_FACE ||
                       res.category == DelaunayMesh::InsertVertexResultCategory::ON_EDGE) {
                mesh_.EnforceLocallyDelaunay(res.i_qe);
                vertex_indices[i_vertex] = mesh_.GetVertexIndex(res.i_qe);
            }
        }

        // Now ensure that there is a line between A and B.
        // If there is not, flip edges until there is.
        QuarterEdgeIndex qe_a = mesh_.GetQuarterEdge(vertex_indices[linedef.i_vertex_start]);
        QuarterEdgeIndex qe_b = mesh_.GetQuarterEdge(vertex_indices[linedef.i_vertex_end]);
        if (!IsValid(qe_a) || !IsValid(qe_b)) {
            return false;
        }

        QuarterEdgeIndex qe_ab = mesh_.EnforceEdge(qe_a, qe_b);
        if (!IsValid(qe_ab)) {
            return false;
        }

        mesh_.ConstrainEdge(qe_ab);

        // Next, assign sidedefs as appropriate.
        qe_a = qe_ab;
        qe_b = mesh_.Sym(qe_a);

        // Always asign the front sidedef.
        // Note that DOOM convention has the opposite handedness, so we add qe_b.
        if (!AddSideInfo(qe_b)) {
            return false;
        }
        SideInfo* side_info = GetEditableSideInfo(qe_b);
        Sidedef sidedef = *(Sidedef*)(sidedefs_data + linedef.i_sidedef_front * sizeof(Sidedef));

        side_info->sector_id = sidedef.sector_index;
        side_info->texture_info_lower.texture_id =
            FindTextureIdForDoomTextureName(sidedef.texture_name_lower, render_assets).value_or(0);
        side_info->texture_info_lower.x_offset = sidedef.offset_x;
        side_info->texture_info_lower.y_offset = sidedef.offset_y;
        side_info->texture_info_upper.texture_id =
            FindTextureIdForDoomTextureName(sidedef.texture_name_upper, render_assets).value_or(0);
        side_info->texture_info_upper.x_offset = sidedef.offset_x;
        side_info->texture_info_upper.y_offset = sidedef.offset_y;
        side_info->texture_info_middle.texture_id =
            FindTextureIdForDoomTextureName(sidedef.texture_name_middle, render_assets).value_or(0);
        side_info->texture_info_middle.x_offset = sidedef.offset_x;
        side_info->texture_info_middle.y_offset = sidedef.offset_y;
        side_info->qe_portal = {kInvalidIndex};

        // Consider the sidedef passable if its middle texture is "AASTINKY" or "-".
        bool is_passable = strncmp((const char*)sidedef.texture_name_middle, "AASTINKY", 8) == 0 ||
                           strncmp((const char*)sidedef.texture_name_middle, "-", 8) == 0;
        if (is_passable) {
            side_info->flags |= core::kSideInfoFlag_PASSABLE;
        }

        // If there is a sidedef on the other side, add that.
        constexpr u16 LINEDEF_FLAG_TWO_SIDED = 0x0004;
        if ((linedef.flags & LINEDEF_FLAG_TWO_SIDED) > 0) {
            // @efficiency - reuse the previous code for inserting a sideinfo.
            if (!AddSideInfo(qe_a)) {
                return false;
            }
            SideInfo* side_info = GetEditableSideInfo(qe_a);
            Sidedef sidedef = *(Sidedef*)(sidedefs_data + linedef.i_sidedef_back * sizeof(Sidedef));

            side_info->sector_id = sidedef.sector_index;
            side_info->texture_info_lower.texture_id =
                FindTextureIdForDoomTextureName(sidedef.texture_name_lower, render_assets)
                    .value_or(0);
            side_info->texture_info_lower.x_offset = sidedef.offset_x;
            side_info->texture_info_lower.y_offset = sidedef.offset_y;
            side_info->texture_info_upper.texture_id =
                FindTextureIdForDoomTextureName(sidedef.texture_name_upper, render_assets)
                    .value_or(0);
            side_info->texture_info_upper.x_offset = sidedef.offset_x;
            side_info->texture_info_upper.y_offset = sidedef.offset_y;
            side_info->texture_info_middle.texture_id =
                FindTextureIdForDoomTextureName(sidedef.texture_name_middle, render_assets)
                    .value_or(0);
            side_info->texture_info_middle.x_offset = sidedef.offset_x;
            side_info->texture_info_middle.y_offset = sidedef.offset_y;
            side_info->qe_portal = {kInvalidIndex};

            if (is_passable) {
                side_info->flags |= core::kSideInfoFlag_PASSABLE;
            }
        }
    }

    // Load sectors
    struct DoomSector {
        i16 z_floor;
        i16 z_ceil;
        u8 texture_name_floor[8];
        u8 texture_name_ceil[8];
        i16 light_level;
        i16 special_type;
        i16 tag_number;
    };
    u32 n_sectors = sectors_data_size / sizeof(DoomSector);
    u32 sectors_data_offset = 0;
    sectors_.resize(n_sectors);
    for (u32 i_sector = 0; i_sector < n_sectors; i_sector++) {
        Sector& sector = sectors_[i_sector];
        DoomSector doom_sector = *(DoomSector*)(sectors_data + sectors_data_offset);
        sectors_data_offset += sizeof(DoomSector);

        sector.flags = 0x0000;  // TODO
        sector.z_floor = doom_sector.z_floor / kDoomUnitsPerWorldUnit;
        sector.z_ceil = doom_sector.z_ceil / kDoomUnitsPerWorldUnit;
    }

    // Segs and Subsectors
    // struct Seg {
    //     i16 i_vertex_start;
    //     i16 i_vertex_end;
    //     i16 discrete_angle;
    //     i16 i_linedef;
    //     i16 is_along_linedef;
    //     i16 offset;  // distance along linedef to start of seg
    // };
    // u32 n_segs = segs_data_size / sizeof(Seg);

    // struct SubsectorEntry {
    //     i16 n_segs;
    //     i16 i_seg_start;
    // };
    // u32 n_subsectors = subsectors_data_size / sizeof(SubsectorEntry);

    // std::vector<common::Vec2f> sector_vs;
    // SubsectorEntry subsector = *(SubsectorEntry*)(subsectors_data + 0);
    // for (i16 i_seg = subsector.i_seg_start; i_seg < subsector.i_seg_start + subsector.n_segs;
    //      i_seg++) {
    //     Seg seg = *(Seg*)(segs_data + sizeof(Seg) * i_seg);
    //     common::Vec2f v_start = vertices_[seg.i_vertex_start].v;
    //     common::Vec2f b = vertices_[seg.i_vertex_end].v;
    //     common::Vec2f dir = common::Normalize(b - v_start);
    //     common::Vec2f a = v_start + dir * seg.offset;
    //     sector_vs.push_back(a);
    //     sector_vs.push_back(b);
    // }

    return true;
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex GameMap::FindVertexNearPosition(const common::Vec2f& pos, QuarterEdgeIndex qe_dual,
                                                 f32 tolerance) const {
    const auto [qe_ab, qe_bc, qe_ca] = mesh_.GetTriangleQuarterEdges(qe_dual);
    const common::Vec2f& a = mesh_.GetVertex(qe_ab);
    const common::Vec2f& b = mesh_.GetVertex(qe_bc);
    const common::Vec2f& c = mesh_.GetVertex(qe_ca);

    f32 dist_a = common::Norm(a - pos);
    f32 dist_b = common::Norm(b - pos);
    f32 dist_c = common::Norm(c - pos);

    QuarterEdgeIndex selected_index = {kInvalidIndex};
    f32 min_distance = tolerance;
    if (dist_a < min_distance) {
        selected_index = qe_ab;
        min_distance = dist_a;
    }
    if (dist_b < min_distance) {
        selected_index = qe_bc;
        min_distance = dist_b;
    }
    if (dist_c < min_distance) {
        selected_index = qe_ca;
        // min_distance = dist_c;
    }
    return selected_index;
}

// ------------------------------------------------------------------------------------------------
QuarterEdgeIndex GameMap::FindEdgeNearPosition(const common::Vec2f& pos, QuarterEdgeIndex qe_dual,
                                               f32 tolerance) const {
    const auto [qe_ab, qe_bc, qe_ca] = mesh_.GetTriangleQuarterEdges(qe_dual);
    const common::Vec2f& a = mesh_.GetVertex(qe_ab);
    const common::Vec2f& b = mesh_.GetVertex(qe_bc);
    const common::Vec2f& c = mesh_.GetVertex(qe_ca);

    f32 dist_ab = common::GetDistanceToLine(pos, a, b);
    f32 dist_bc = common::GetDistanceToLine(pos, b, c);
    f32 dist_ca = common::GetDistanceToLine(pos, c, a);

    QuarterEdgeIndex qe_selected = {kInvalidIndex};
    f32 min_distance = tolerance;
    if (dist_ab < min_distance) {
        qe_selected = qe_ab;
        min_distance = dist_ab;
    }
    if (dist_bc < min_distance) {
        qe_selected = qe_bc;
        min_distance = dist_bc;
    }
    if (dist_ca < min_distance) {
        qe_selected = qe_ca;
        min_distance = dist_ca;
    }

    // If this edge does not have a side_info, but the reversed edge does, select that instead.
    if (IsValid(qe_selected) && side_infos_.find(qe_selected) == side_infos_.end()) {
        QuarterEdgeIndex qe_rev = mesh_.Sym(qe_selected);
        if (side_infos_.find(qe_rev) != side_infos_.end()) {
            qe_selected = qe_rev;
        }
    }

    return qe_selected;
}

}  // namespace core