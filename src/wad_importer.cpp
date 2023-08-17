#include "wad_importer.hpp"

#include <cstdio>
#include <cstring>

namespace core {

// ------------------------------------------------------------------------------------------------
std::unique_ptr<WadImporter> WadImporter::LoadFromFile(std::string filepath) {
    std::unique_ptr<WadImporter> ptr = {};

    FILE* fileptr = fopen(filepath.c_str(), "rb");
    if (fileptr == nullptr) {
        return ptr;
    }

    // Count the number of bytes
    fseek(fileptr, 0, SEEK_END);
    u32 n_bytes_in_file = ftell(fileptr);

    if (n_bytes_in_file < 8) {
        return ptr;  // too few
    }

    // Read in the binary assets as a single blob, including the header
    fseek(fileptr, 0, SEEK_SET);
    u8* blob = (u8*)malloc(n_bytes_in_file);
    if (blob == nullptr) {
        return ptr;  // failed to allocate blob
    }

    usize n_bytes_read = fread(blob, sizeof(u8), n_bytes_in_file, fileptr);
    if (n_bytes_read != n_bytes_in_file) {
        free(blob);
        return ptr;
    }

    return std::make_unique<WadImporter>(blob, n_bytes_in_file);
}

// ------------------------------------------------------------------------------------------------
WadImporter::WadImporter(u8* blob, u32 n_bytes) : blob_(blob), blob_size_(n_bytes) {
    // Process the loaded assets from the loaded binary blob

    // Read the number of table of content entries
    n_entries_ = *(u32*)(blob_ + 0x04);
    u32 dir_loc = *(u32*)(blob_ + 0x08);

    // Set the blob entries
    entries_ = (WadDirectoryEntry*)(blob_ + dir_loc);
}

// ------------------------------------------------------------------------------------------------
WadImporter::~WadImporter() { free(blob_); }

// ------------------------------------------------------------------------------------------------
std::optional<const u8*> WadImporter::FindEntryData(const std::string& entry_name) const {
    for (u32 i = 0; i < n_entries_; i++) {
        WadDirectoryEntry* entry = entries_ + i;
        if (strncmp(entry->name, entry_name.c_str(), 8) == 0) {
            return blob_ + entry->byte_offset;
        }
    }
    return std::nullopt;
}

// ------------------------------------------------------------------------------------------------
std::optional<int> WadImporter::FindEntryDataIndex(const std::string& entry_name,
                                                   int index_start) const {
    for (int i = index_start; i < (int)n_entries_; i++) {
        WadDirectoryEntry* entry = entries_ + i;
        if (strncmp(entry->name, entry_name.c_str(), 8) == 0) {
            return i;
        }
    }
    return std::nullopt;
}

// ------------------------------------------------------------------------------------------------
std::optional<std::pair<const u8*, u32>> WadImporter::GetEntryData(int index) const {
    if (0 <= index && index < (int)n_entries_) {
        WadDirectoryEntry* entry = entries_ + index;
        const u8* ptr = blob_ + entry->byte_offset;
        return std::make_pair(ptr, entry->size);
    }
    return std::nullopt;
}

}  // namespace core