#include "assets_importer.hpp"

#include <cstdio>
#include <cstring>

namespace core {

// ------------------------------------------------------------------------------------------------
int OldStyleBitmap::GetColumnMajorPixelIndex(int x, int y) const {
    return y + x * n_pixels_per_column;
}

// ------------------------------------------------------------------------------------------------
int OldStyleBitmap::GetRowMajorPixelIndex(int x, int y) const { return x + y * n_pixels_per_row; }

// ------------------------------------------------------------------------------------------------
u32 OldStyleBitmap::GetColumnMajorPixelAt(int x, int y) const {
    return abgr[y + x * n_pixels_per_column];
}

// ------------------------------------------------------------------------------------------------
u32 OldStyleBitmap::GetRowMajorPixelAt(int x, int y) const {
    return abgr[x + y * n_pixels_per_row];
}

// ------------------------------------------------------------------------------------------------
OldStyleBitmap LoadBitmap(u8* data) {
    OldStyleBitmap bitmap = {};

    bitmap.n_pixels = *(u32*)(data);
    data += sizeof(u32);
    bitmap.n_pixels_per_column = *(u32*)(data);
    data += sizeof(u32);
    bitmap.column_major = *data;
    data += sizeof(u8);
    bitmap.abgr = (u32*)(data);
    bitmap.n_pixels_per_row = bitmap.n_pixels / bitmap.n_pixels_per_column;

    return bitmap;
}

// ------------------------------------------------------------------------------------------------
std::unique_ptr<AssetsImporter> AssetsImporter::LoadFromFile(std::string filepath) {
    std::unique_ptr<AssetsImporter> ptr = {};

    FILE* fileptr = fopen(filepath.c_str(), "rb");
    if (fileptr == nullptr) {
        return ptr;
    }

    // Count the number of bytes
    fseek(fileptr, 0, SEEK_END);
    u32 n_bytes_in_file = ftell(fileptr);
    u32 n_bytes_past_header = n_bytes_in_file - 4;

    if (n_bytes_in_file < 8) {
        return ptr;  // too few
    }

    // Read in the binary assets as a single blob
    fseek(fileptr, 4, SEEK_SET);  // Skip the header bytes
    u8* blob = (u8*)malloc(n_bytes_past_header);
    if (blob == nullptr) {
        return ptr;  // failed to allocate blob
    }

    usize n_bytes_read = fread(blob, sizeof(u8), n_bytes_past_header, fileptr);
    if (n_bytes_read != n_bytes_past_header) {
        free(blob);
        return ptr;
    }

    return std::make_unique<AssetsImporter>(blob, n_bytes_past_header);
}

// ------------------------------------------------------------------------------------------------
AssetsImporter::AssetsImporter(u8* blob, u32 n_bytes) : blob_(blob), blob_size_(n_bytes) {
    // Process the loaded assets from the loaded binary blob

    // Read the number of table of content entries
    u32 byte_index = n_bytes - sizeof(u32);
    n_entries_ = *(u32*)(blob_ + byte_index);

    // Scan through them in reverse order
    for (int i = n_entries_; i > 0; i--) {
        byte_index -= sizeof(AssetsImporterEntry);
        entries_ = (AssetsImporterEntry*)(blob_ + byte_index);
        // Make the name null-terminated just in case.
        entries_->name[ASSET_NAME_BYTE_COUNT - 1] = '\0';
        printf("Entry %d: %s at offset %d\n", i, entries_->name, entries_->byte_offset);
    }
}

// ------------------------------------------------------------------------------------------------
AssetsImporter::~AssetsImporter() { free(blob_); }

// ------------------------------------------------------------------------------------------------
std::optional<u8*> AssetsImporter::FindEntryData(const std::string& entry_name) const {
    for (u32 i = 0; i < n_entries_; i++) {
        AssetsImporterEntry* entry = entries_ + i;
        if (strcmp(entry->name, entry_name.c_str()) == 0) {
            return blob_ + entry->byte_offset;
        }
    }
    return std::nullopt;
}

}  // namespace core