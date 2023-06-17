#include "assets_exporter.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

namespace core {

struct ExportedTocEntry {
    u32 offset;
    u32 size;
    char name[ASSET_NAME_BYTE_COUNT];
};

// --------------------------------------
void SetEntryName(AssetsExporterEntry& entry, const std::string& name) {
    memset(entry.name, 0, ASSET_NAME_BYTE_COUNT);
    memcpy(entry.name, name.data(), std::min<size_t>(ASSET_NAME_BYTE_COUNT, name.size()));
}

// ------------------------------------------------------------------------------------------------
AssetsExporter::AssetsExporter() {}

// ------------------------------------------------------------------------------------------------
int AssetsExporter::FindEntry(const std::string& entry_name) {
    for (int i = 0; i < (int)(entries_.size()); i++) {
        if (strcmp(entries_[i].name, entry_name.c_str()) == 0) {
            return i;
        }
    }
    return -1;
}

// ------------------------------------------------------------------------------------------------
bool AssetsExporter::HasEntry(const std::string& entry_name) { return FindEntry(entry_name) >= 0; }

// ------------------------------------------------------------------------------------------------
void AssetsExporter::Clear() { entries_.clear(); }

// ------------------------------------------------------------------------------------------------
void AssetsExporter::AddEntry(const AssetsExporterEntry& entry) { entries_.emplace_back(entry); }

// ------------------------------------------------------------------------------------------------
void AssetsExporter::AddMeshEntry(const DelaunayMesh& mesh, const std::string& name) {
    AssetsExporterEntry entry;
    SetEntryName(entry, name);

    // --------------------------------------
    // Serialize the mesh
    std::stringstream buffer;

    // Serialize the counts
    u32 n_vertices = mesh.NumVertices();
    buffer.write(reinterpret_cast<const char*>(&n_vertices), sizeof(u32));
    u32 n_quarter_edges = mesh.NumQuarterEdges();
    buffer.write(reinterpret_cast<const char*>(&n_quarter_edges), sizeof(u32));

    // Serialize the vertices
    for (size_t i = 0; i < n_vertices; i++) {
        const common::Vec2f& v = mesh.GetVertex(i);
        buffer.write(reinterpret_cast<const char*>(&v.x), sizeof(f32));
        buffer.write(reinterpret_cast<const char*>(&v.y), sizeof(f32));
    }

    // Serialize the quarter edges
    for (size_t i = 0; i < n_quarter_edges; i++) {
        const core::QuarterEdge* qe = mesh.GetQuarterEdge(i);

        u32 vertex_index = std::numeric_limits<u32>::max();
        if (qe->vertex != nullptr) {
            vertex_index = qe->vertex->index;
        }
        buffer.write(reinterpret_cast<const char*>(&vertex_index), sizeof(u32));
        buffer.write(reinterpret_cast<const char*>(&qe->next->index), sizeof(u32));
        buffer.write(reinterpret_cast<const char*>(&qe->rot->index), sizeof(u32));
    }

    // --------------------------------------
    // Set the entry data from the buffer.
    std::string buffer_str = buffer.str();
    entry.data.reserve(buffer_str.size());
    const u8* buffer_data = reinterpret_cast<const u8*>(buffer_str.data());
    for (size_t i = 0; i < buffer_str.size(); i++) {
        entry.data.push_back(buffer_data[i]);
    }

    AddEntry(entry);
}

// ------------------------------------------------------------------------------------------------
bool AssetsExporter::LoadMeshEntry(DelaunayMesh* mesh, const std::string& name) {
    int entry_index = FindEntry(name);
    if (entry_index < 0) {
        std::cout << "Failed to find entry " << name << std::endl;
        return false;
    }

    const u8* data = entries_[entry_index].data.data();
    return mesh->LoadFromData(data);
}

// ------------------------------------------------------------------------------------------------
void AssetsExporter::AddSideInfos(const GameMap& game_map, const std::string& name) {
    AssetsExporterEntry entry;
    SetEntryName(entry, name);

    // --------------------------------------
    // Serialize the content
    std::stringstream buffer;

    // Serialize the counts
    u32 n_side_infos = game_map.side_infos.size();
    buffer.write(reinterpret_cast<const char*>(&n_side_infos), sizeof(u32));

    // Serialize the side infos
    struct ExportedSideInfo {
        u16 flags;
        u16 texture_id;
        i16 x_offset;
        i16 y_offset;
    };

    for (const auto& side_info : game_map.side_infos) {
        ExportedSideInfo exported_side_info = {.flags = side_info.flags,
                                               .texture_id = side_info.texture_id,
                                               .x_offset = side_info.x_offset,
                                               .y_offset = side_info.y_offset};
        buffer.write(reinterpret_cast<const char*>(&exported_side_info),
                     sizeof(exported_side_info));
    }

    // Serialize a vector[u16] of quarter edge index -> side info index.
    // Uses 0xFFFF to indicate no quarter edge.
    for (size_t i = 0; i < game_map.mesh.NumQuarterEdges(); i++) {
        const QuarterEdge* qe = game_map.mesh.GetQuarterEdge(i);
        u16 side_info_index = std::numeric_limits<u16>::max();
        if (IsPrimalEdge(*qe) && !game_map.mesh.IsBoundaryVertex(qe->vertex) &&
            !game_map.mesh.IsBoundaryVertex(game_map.mesh.Sym(qe)->vertex)) {
            size_t a_ind = MeshVertexIndexToMapVertexIndex(qe->vertex->index);
            size_t b_ind = MeshVertexIndexToMapVertexIndex(game_map.mesh.Sym(qe)->vertex->index);
            auto it_side = game_map.side_to_info.find(std::make_pair(a_ind, b_ind));
            if (it_side != game_map.side_to_info.end()) {
                side_info_index = (u16)(it_side->second);
            }
        }
        buffer.write(reinterpret_cast<const char*>(&side_info_index), sizeof(u16));
    }

    // --------------------------------------
    // Set the entry data from the buffer.
    std::string buffer_str = buffer.str();
    entry.data.reserve(buffer_str.size());
    const u8* buffer_data = reinterpret_cast<const u8*>(buffer_str.data());
    for (size_t i = 0; i < buffer_str.size(); i++) {
        entry.data.push_back(buffer_data[i]);
    }

    AddEntry(entry);
}

// ------------------------------------------------------------------------------------------------
bool AssetsExporter::WriteToFile(const std::string& filepath) {
    std::ofstream fout;
    fout.open(filepath, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!fout.is_open()) {
        return false;
    }

    char header[4] = {'T', 'O', 'O', 'M'};
    fout.write(reinterpret_cast<const char*>(&header), 4);

    // Write out the number of entries
    u32 n_entries = NumEntries();
    fout.write(reinterpret_cast<const char*>(&n_entries), sizeof(u32));

    // Temporarily skip the table of contents location
    auto toc_location = fout.tellp();
    u32 toc_offset = 0;
    fout.write(reinterpret_cast<const char*>(&toc_offset), sizeof(u32));

    // Write out all entry data, keeping track of its offset in the file
    std::vector<u32> entry_offsets(n_entries);
    for (u32 i = 0; i < n_entries; i++) {
        entry_offsets[i] = fout.tellp();
        const AssetsExporterEntry& entry = entries_[i];
        fout.write(reinterpret_cast<const char*>(entry.data.data()), entry.data.size());
    }

    // Write the table of contents
    toc_offset = fout.tellp();
    for (u32 i = 0; i < n_entries; i++) {
        ExportedTocEntry output_entry;
        output_entry.offset = entry_offsets[i];
        output_entry.size =
            (i < n_entries - 1 ? entry_offsets[i + 1] : toc_offset) - output_entry.offset;
        memcpy(output_entry.name, entries_[i].name, ASSET_NAME_BYTE_COUNT);
        fout.write(reinterpret_cast<const char*>(&output_entry), sizeof(output_entry));
    }

    // Set the table of contents location
    fout.seekp(toc_location);
    fout.write(reinterpret_cast<const char*>(&toc_offset), sizeof(u32));

    fout.close();
    return true;
}

// ------------------------------------------------------------------------------------------------
bool AssetsExporter::LoadAssetsFile(const std::string& filepath) {
    Clear();

    FILE* fileptr = fopen(filepath.c_str(), "rb");
    if (fileptr == nullptr) {
        std::cout << "Error opening assets file [" << filepath << "]" << std::endl;
        return false;
    }

    // Count the number of bytes
    fseek(fileptr, 0, SEEK_END);
    u32 n_bytes_in_file = ftell(fileptr);

    // Read in the entire file
    fseek(fileptr, 0, SEEK_SET);
    u8* bytes = (u8*)malloc(n_bytes_in_file);
    if (bytes == nullptr) {
        std::cout << "Failed to allocate memory for the assets file [" << filepath << "]"
                  << std::endl;
        return false;
    }
    u32 bytes_read = fread(bytes, sizeof(u8), n_bytes_in_file, fileptr);
    if (bytes_read != n_bytes_in_file) {
        std::cout << "Failed to read in assets file [" << filepath << "]!" << std::endl;
        std::cout << "Read " << bytes_read << " bytes but expected " << n_bytes_in_file
                  << std::endl;
    }

    u32 n_entries = *(u32*)(bytes + 0x04);
    u32 toc_offset = *(u32*)(bytes + 0x08);

    // Process the table of contents
    u32 offset = toc_offset;
    for (u32 toc_index = 0; toc_index < n_entries; toc_index++) {
        ExportedTocEntry* imported_entry = (ExportedTocEntry*)(bytes + offset);
        printf("Entry %d: %.16s at offset %d with size %d\n", toc_index, imported_entry->name,
               imported_entry->offset, imported_entry->size);
        offset += sizeof(ExportedTocEntry);

        // Pull its data
        AssetsExporterEntry assets_entry;
        memcpy(assets_entry.name, imported_entry->name, ASSET_NAME_BYTE_COUNT);
        assets_entry.data.resize(imported_entry->size);
        memcpy(assets_entry.data.data(), bytes + imported_entry->offset, imported_entry->size);
        AddEntry(assets_entry);
    }

    fclose(fileptr);
    return true;
}

}  // namespace core