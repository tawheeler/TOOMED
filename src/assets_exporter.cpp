#include "assets_exporter.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>

namespace core {

struct ExportedTocEntry {
    u32 offset;
    u32 size;
    char name[ASSET_NAME_BYTE_COUNT];
};

// --------------------------------------
void AssetsExporterEntry::SetName(const std::string& newname) {
    memset(name, 0, ASSET_NAME_BYTE_COUNT);
    memcpy(name, newname.data(), std::min<size_t>(ASSET_NAME_BYTE_COUNT, newname.size()));
}

// --------------------------------------
void AssetsExporterEntry::SetData(const std::string& buffer_str) {
    const u8* buffer_data = reinterpret_cast<const u8*>(buffer_str.data());
    data.resize(buffer_str.size());
    memcpy(data.data(), buffer_data, buffer_str.size());
}

// ------------------------------------------------------------------------------------------------
AssetsExporter::AssetsExporter() {}

// ------------------------------------------------------------------------------------------------
const AssetsExporterEntry* AssetsExporter::FindEntry(const std::string& entry_name) const {
    for (int i = 0; i < (int)(entries_.size()); i++) {
        if (strcmp(entries_[i].name, entry_name.c_str()) == 0) {
            return &(entries_[i]);
        }
    }
    return nullptr;
}

// ------------------------------------------------------------------------------------------------
bool AssetsExporter::HasEntry(const std::string& entry_name) const {
    return FindEntry(entry_name) != nullptr;
}

// ------------------------------------------------------------------------------------------------
void AssetsExporter::Clear() { entries_.clear(); }

// ------------------------------------------------------------------------------------------------
void AssetsExporter::AddEntry(const AssetsExporterEntry& entry) { entries_.emplace_back(entry); }

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