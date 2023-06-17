#pragma once

#include <map>
#include <vector>

#include "delaunay_mesh.hpp"
#include "typedefs.hpp"

#define MESH_BOUNDING_RADIUS 1000.0f
#define MESH_MIN_DIST_TO_VERTEX 0.1f
#define MESH_MIN_DIST_TO_EDGE 0.1f

namespace core {

// The information associated with one side of an edge between vertices in the map.
// If this represents the directed edge A->B, then it describes the edge viewed on the right side of
// A->B.
struct SideInfo {
    u16 flags;
    u16 texture_id;
    i16 x_offset;  // Texture x offset
    i16 y_offset;  // Texture y offset
    usize a_ind;   // The index of vertex A (the source vertex)
    usize b_ind;   // The index of vertex B (the dest vertex)
};

// Represents our game map
struct GameMap {
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

    void Clear() {
        vertices.clear();
        side_infos.clear();
        side_to_info.clear();
        mesh.Clear();
    }

    bool HasEdge(int a_ind, int b_ind) {
        auto tup = std::make_pair(a_ind, b_ind);
        return side_to_info.find(tup) != side_to_info.end();
    }

    usize AddEdge(int a_ind, int b_ind) {
        usize side_info_index = side_infos.size();
        SideInfo side_info;
        side_info.a_ind = a_ind;
        side_info.b_ind = b_ind;
        side_infos.emplace_back(side_info);
        side_to_info[std::make_pair(a_ind, b_ind)] = side_info_index;
        return side_info_index;
    }
};

int MapVertexIndexToMeshVertexIndex(int ind);
int MeshVertexIndexToMapVertexIndex(int ind);

}  // namespace core