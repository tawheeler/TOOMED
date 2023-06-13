#include "assets_exporter.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

namespace core {

// ------------------------------------------------------------------------------------------------
AssetsExporter::AssetsExporter() {}

// ------------------------------------------------------------------------------------------------
void AssetsExporter::AddEntry(const AssetsExporterEntry& entry) { entries_.emplace_back(entry); }

// ------------------------------------------------------------------------------------------------
void AssetsExporter::AddMeshEntry(const DelaunayMesh& mesh, std::string name) {
    AssetsExporterEntry entry;

    // Store the name
    memset(entry.name, 0, 16);
    memcpy(entry.name, name.data(), std::min<size_t>(16, name.size()));

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
bool AssetsExporter::WriteToFile(const std::string& filename) {
    std::ofstream fout;
    fout.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
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
    struct TocEntry {
        u32 offset;
        u32 size;
        char name[16];
    };

    toc_offset = fout.tellp();
    for (u32 i = 0; i < n_entries; i++) {
        TocEntry output_entry;
        output_entry.offset = entry_offsets[i];
        output_entry.size =
            (i < n_entries - 1 ? entry_offsets[i + 1] : toc_offset) - output_entry.offset;
        memcpy(output_entry.name, entries_[i].name, 16);
        fout.write(reinterpret_cast<const char*>(&output_entry), sizeof(output_entry));
    }

    // Set the table of contents location
    fout.seekp(toc_location);
    fout.write(reinterpret_cast<const char*>(&toc_offset), sizeof(u32));

    fout.close();
    return true;
}

}  // namespace core