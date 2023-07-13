#include "game_map.hpp"

#include <iostream>
#include <sstream>

#include "math_utils.hpp"

namespace core {

struct ExportedSideInfo {
    u16 flags;
    u16 texture_id;
    i16 x_offset;
    i16 y_offset;
};

// ------------------------------------------------------------------------------------------------
GameMap::GameMap() : mesh_(MESH_BOUNDING_RADIUS, MESH_MIN_DIST_TO_VERTEX, MESH_MIN_DIST_TO_EDGE) {}

// ------------------------------------------------------------------------------------------------
void GameMap::Clear() {
    mesh_.Clear();
    side_infos_.clear();
}

// //
// ------------------------------------------------------------------------------------------------
// bool GameMap::HasEdge(int a_ind, int b_ind) const {
//     auto tup = std::make_pair(a_ind, b_ind);
//     return side_to_info_.find(tup) != side_to_info_.end();
// }

// //
// ------------------------------------------------------------------------------------------------
// bool GameMap::HasMesh() const {
//     if (mesh_) {
//         return true;
//     }
//     return false;
// };

// //
// ------------------------------------------------------------------------------------------------
// std::optional<usize> GameMap::GetEdgeIndex(int a_ind, int b_ind) const {
//     auto tup = std::make_pair(a_ind, b_ind);
//     auto it = side_to_info_.find(tup);
//     if (it != side_to_info_.end()) {
//         return it->second;
//     }
//     return std::nullopt;
// }

// ------------------------------------------------------------------------------------------------
SideInfo* GameMap::GetEditableSideInfo(QuarterEdgeIndex qe_primal) {
    auto it = side_infos_.find(qe_primal);
    if (it == side_infos_.end()) {
        return nullptr;
    }
    return &(it->second);
}
// ------------------------------------------------------------------------------------------------
FaceInfo* GameMap::GetEditableFaceInfo(QuarterEdgeIndex qe_dual) {
    auto it = face_infos_.find(qe_dual);
    if (it == face_infos_.end()) {
        return nullptr;
    }
    return &(it->second);
}

// ------------------------------------------------------------------------------------------------
VertexIndex GameMap::AddVertex(const common::Vec2f& pos) {
    DelaunayMesh::InsertVertexResult result = mesh_.InsertVertex(pos);

    // TODO: Handle side_info data based on how the vertex was inserted.
    // TODO: Handle solid triangle info based on how the vertex was inserted.
    // TODO: Probably want to automatically enforce Delaunay-ness.

    if (IsValid(result.i_vertex)) {
        mesh_.EnforceLocallyDelaunay(result.i_qe);
    }

    return result.i_vertex;
}

// //
// ------------------------------------------------------------------------------------------------
// usize GameMap::AddDirectedEdge(int a_ind, int b_ind) {
//     usize side_info_index = side_infos_.size();
//     SideInfo side_info;
//     side_info.a_ind = a_ind;
//     side_info.b_ind = b_ind;
//     side_infos_.emplace_back(side_info);
//     side_to_info_[std::make_pair(a_ind, b_ind)] = side_info_index;

//     // Invalidate the mesh if this is not aleady an edge in our mesh.
//     if (HasMesh()) {
//         int mesh_ind_a = MapVertexIndexToMeshVertexIndex(a_ind);
//         int mesh_ind_b = MapVertexIndexToMeshVertexIndex(b_ind);
//         if (mesh_->GetQuarterEdge(mesh_ind_a, mesh_ind_b) == nullptr) {
//             InvalidateMesh();
//         }
//     }

//     // InvalidateMesh();  // We will need to re-generate the mesh
//     return side_info_index;
// }

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
bool GameMap::AddFaceInfo(QuarterEdgeIndex qe_dual) {
    if (!mesh_.IsDual(qe_dual)) {
        return false;
    }

    auto it = face_infos_.find(qe_dual);
    if (it != face_infos_.end()) {
        return false;
    }

    FaceInfo face_info = {};
    face_info.qe = qe_dual;

    face_infos_[qe_dual] = (face_info);
    QuarterEdgeIndex qe = mesh_.Next(qe_dual);
    while (qe != qe_dual) {
        face_info.qe = qe;
        face_infos_[qe] = (face_info);
        qe = mesh_.Next(qe);
    }

    return true;
}

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
                                               .texture_id = side_info.texture_id,
                                               .x_offset = side_info.x_offset,
                                               .y_offset = side_info.y_offset};

        // // TODO: Remove this calculation
        // // Calculate the angle of the side info, and set it to shaded depending.
        // common::Vec2f a = vertices_[side_info.a_ind];
        // common::Vec2f b = vertices_[side_info.b_ind];
        // float angle = atan2(b.y - a.y, b.x - a.x);
        // float angledist =
        //     std::min(common::AngleDist(angle, 0.0f), common::AngleDist(angle, 3.14159265f));
        // if (angledist > M_PI / 4) {
        //     // shaded
        //     exported_side_info.flags |= kSideInfoFlag_DARK;
        // }

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
bool GameMap::Export(AssetsExporter* exporter) const {
    // First, write out the geometry_mesh entry for the delaunay mesh
    exporter->AddEntry(ExportDelaunayMesh(kAssetEntryGeometryMesh));

    // Then write out all side info definitions
    exporter->AddEntry(ExportSideInfos(kAssetEntrySideInfos));

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
        side_info.texture_id = exported_side_info.texture_id;
        side_info.x_offset = exported_side_info.x_offset;
        side_info.y_offset = exported_side_info.y_offset;
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

    return success;
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