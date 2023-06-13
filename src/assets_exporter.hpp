#pragma once

#include <string>
#include <vector>

#include "delaunay_mesh.hpp"
#include "typedefs.hpp"

namespace core {

struct AssetsExporterEntry {
    char name[16];
    std::vector<u8> data;
};

class AssetsExporter {
  public:
    AssetsExporter();
    ~AssetsExporter() = default;

    size_t NumEntries() const { return entries_.size(); }

    void AddEntry(const AssetsExporterEntry& entry);

    void AddMeshEntry(const DelaunayMesh& mesh, std::string name);

    bool WriteToFile(const std::string& filename);

  private:
    std::vector<AssetsExporterEntry> entries_ = {};
};

}  // namespace core