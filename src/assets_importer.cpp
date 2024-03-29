#include "assets_importer.hpp"

#include <cstdio>
#include <cstring>
#include <sstream>

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
std::vector<doom::Patch> ExtractPatches(const OldStyleBitmap& bitmap, const std::string& prefix,
                                        const Palette& palette) {
    std::vector<doom::Patch> patches;

    std::vector<LAB> lab_palette = GetLabPalette(palette);

    constexpr u32 kTextureSize = 64;

    u32 n_patches = bitmap.n_pixels_per_column / kTextureSize;
    patches.reserve(n_patches);
    for (u32 i_patch; i_patch < n_patches; i_patch++) {
        doom::Patch patch;

        std::stringstream ss;
        ss << prefix;
        if (i_patch < 10) {
            ss << "0";
        }
        ss << i_patch;

        patch.name = ss.str();
        patch.size_x = kTextureSize;
        patch.size_y = kTextureSize;

        // Just default the origin to the bottom center
        patch.origin_x = patch.size_x / 2;
        patch.origin_y = patch.size_y;

        patch.column_offsets.resize(patch.size_x, 0);

        u32 idx = 0;
        for (u16 x = 0; x < patch.size_x; x++) {
            // Create the posts for column x.

            // Just add fresh columns.
            // TODO: @efficiency: reuse repeated columns.
            patch.column_offsets[x] = (u32)(patch.post_data.size());

            bool has_active_post = false;
            u32 byte_offset_for_active_post = 0;
            u16 y_delta = 0;
            u16 length = 0;

            for (u16 y = 0; y < patch.size_y; y++) {
                idx++;
                if (!has_active_post) {
                    // We started an active post.
                    has_active_post = true;
                    length = 1;

                    // start a new active post.
                    byte_offset_for_active_post = (u32)(patch.post_data.size());
                    patch.post_data.push_back(y_delta);
                    patch.post_data.push_back(length);  // will need to overwrite later

                    u32 abgr = bitmap.GetColumnMajorPixelAt(x, y + i_patch * kTextureSize);
                    patch.post_data.push_back(GetClosestColorInPalette(abgr, lab_palette));
                    y_delta = 0;
                } else {
                    length++;
                    u32 abgr = bitmap.GetColumnMajorPixelAt(x, y + i_patch * kTextureSize);
                    patch.post_data.push_back(GetClosestColorInPalette(abgr, lab_palette));
                }
            }

            // Write out the length if we need to terminate the final post
            if (has_active_post) {
                patch.post_data[byte_offset_for_active_post + 1] = length;
            }

            // Terminate each column
            patch.post_data.push_back(0xFF);
        }

        patches.push_back(patch);
    }

    return patches;
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