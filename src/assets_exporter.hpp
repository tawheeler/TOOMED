#pragma once

#include <optional>
#include <string>
#include <vector>

// #include "delaunay_mesh.hpp"
// #include "game_map.hpp"
#include "typedefs.hpp"

#define ASSET_NAME_BYTE_COUNT 16

namespace core {

struct AssetsExporterEntry {
    char name[ASSET_NAME_BYTE_COUNT];
    std::vector<u8> data;

    void SetName(const std::string& newname);
    void SetData(const std::string& buffer_str);
};

class AssetsExporter {
  public:
    AssetsExporter();
    ~AssetsExporter() = default;

    size_t NumEntries() const { return entries_.size(); }

    bool HasEntry(const std::string& entry_name) const;
    const AssetsExporterEntry* FindEntry(const std::string& entry_name) const;

    void Clear();

    void AddEntry(const AssetsExporterEntry& entry);

    bool WriteToFile(const std::string& filepath);

    bool LoadAssetsFile(const std::string& filepath);

  private:
    std::vector<AssetsExporterEntry> entries_ = {};
};

}  // namespace core