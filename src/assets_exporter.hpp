#pragma once

#include <string>
#include <vector>

#include "delaunay_mesh.hpp"
#include "game_map.hpp"
#include "typedefs.hpp"

#define ASSET_NAME_BYTE_COUNT 16

namespace core {

struct AssetsExporterEntry {
    char name[ASSET_NAME_BYTE_COUNT];
    std::vector<u8> data;
};

class AssetsExporter {
  public:
    AssetsExporter();
    ~AssetsExporter() = default;

    size_t NumEntries() const { return entries_.size(); }

    bool HasEntry(const std::string& entry_name);

    void Clear();

    void AddEntry(const AssetsExporterEntry& entry);

    void AddMeshEntry(const DelaunayMesh& mesh, const std::string& name);
    bool LoadMeshEntry(DelaunayMesh* mesh, const std::string& name);

    void AddSideInfos(const GameMap& game_map, const std::string& name);

    bool WriteToFile(const std::string& filepath);

    bool LoadAssetsFile(const std::string& filepath);

  private:
    int FindEntry(const std::string& entry_name);

    std::vector<AssetsExporterEntry> entries_ = {};
};

}  // namespace core