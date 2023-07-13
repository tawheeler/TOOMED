#pragma once

#include <memory>
#include <optional>
#include <string>

#include "typedefs.hpp"

#ifndef ASSET_NAME_BYTE_COUNT
#define ASSET_NAME_BYTE_COUNT 16
#endif

namespace core {

// TODO: Eventually move this elsewhere.
//       We want to move to a system where we index into the palette rather than store full ABGR
//       values.
struct OldStyleBitmap {
    u32 n_pixels;
    u32 n_pixels_per_column;
    u32 n_pixels_per_row;
    bool column_major;
    u32* abgr;

    int GetColumnMajorPixelIndex(int x, int y) const;
    int GetRowMajorPixelIndex(int x, int y) const;
    u32 GetColumnMajorPixelAt(int x, int y) const;
    u32 GetRowMajorPixelAt(int x, int y) const;
};

OldStyleBitmap LoadBitmap(u8* data);

struct AssetsImporterEntry {
    u32 byte_offset;
    char name[ASSET_NAME_BYTE_COUNT];
};

class AssetsImporter {
  public:
    static std::unique_ptr<AssetsImporter> LoadFromFile(std::string filepath);
    AssetsImporter(u8* blob, u32 n_bytes);
    ~AssetsImporter();

    u32 NumEntries() const { return n_entries_; }

    // If an entry with the given name exists, returns a pointer (into blob_)
    // of where that entry's data begins.
    std::optional<u8*> FindEntryData(const std::string& entry_name) const;

  private:
    u8* blob_;
    u32 blob_size_;

    // These are in the blob
    AssetsImporterEntry* entries_;
    u32 n_entries_;
};

}  // namespace core