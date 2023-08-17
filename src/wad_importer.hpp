#pragma once

#include <memory>
#include <optional>
#include <string>

#include "typedefs.hpp"

namespace core {

struct WadDirectoryEntry {
    u32 byte_offset;
    u32 size;
    char name[8];
};

class WadImporter {
  public:
    static std::unique_ptr<WadImporter> LoadFromFile(std::string filepath);
    WadImporter(u8* blob, u32 n_bytes);
    ~WadImporter();

    u32 NumEntries() const { return n_entries_; }

    // If an entry with the given name exists, returns a pointer (into blob_)
    // of where that entry's data begins.
    std::optional<const u8*> FindEntryData(const std::string& entry_name) const;

    std::optional<int> FindEntryDataIndex(const std::string& entry_name, int index_start = 0) const;

    std::optional<std::pair<const u8*, u32>> GetEntryData(int index) const;

  private:
    u8* blob_;
    u32 blob_size_;

    // These are in the blob
    WadDirectoryEntry* entries_;
    u32 n_entries_;
};

}  // namespace core