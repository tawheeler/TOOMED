#include "palette.hpp"

#include <cstring>
#include <iostream>
#include <sstream>

namespace core {

// ------------------------------------------------------------------------------------------------
AssetsExporterEntry ExportPalettes(const std::vector<Palette>& palettes) {
    AssetsExporterEntry entry;
    entry.SetName(kAssetEntryPalettes);

    // Serialize the data
    std::stringstream buffer;

    u32 n_palettes = palettes.size();
    buffer.write(reinterpret_cast<const char*>(&n_palettes), sizeof(u32));

    for (const Palette& palette : palettes) {
        buffer.write(reinterpret_cast<const char*>(&palette.rgbs), sizeof(palette.rgbs));
    }
    entry.SetData(buffer.str());

    return entry;
}

// ------------------------------------------------------------------------------------------------
AssetsExporterEntry ExportColormaps(const std::vector<Colormap>& colormaps) {
    AssetsExporterEntry entry;
    entry.SetName(kAssetEntryColormaps);

    // Serialize the data
    std::stringstream buffer;

    u32 n_colormaps = colormaps.size();
    buffer.write(reinterpret_cast<const char*>(&n_colormaps), sizeof(u32));

    for (const Colormap& colormap : colormaps) {
        buffer.write(reinterpret_cast<const char*>(&colormap.map), sizeof(colormap.map));
    }
    entry.SetData(buffer.str());

    return entry;
}

}  // namespace core

namespace doom {

// ------------------------------------------------------------------------------------------------
std::vector<core::Palette> ParseDoomPalettes(const std::unique_ptr<core::WadImporter>& importer) {
    std::vector<core::Palette> palettes;

    auto data_opt = importer->FindEntryData("PLAYPAL");
    if (!data_opt) {
        return palettes;
    }

    const u8* data = *data_opt;
    u32 data_offset = 0;

    // Pull out the palettes.
    palettes.resize(kNumDoomPalettes);
    for (i32 i = 0; i < kNumDoomPalettes; i++) {
        memcpy(palettes[i].rgbs, data + data_offset, 256 * 3);
        data_offset += 256 * 3;
    }

    return palettes;
}

// ------------------------------------------------------------------------------------------------
std::vector<core::Colormap> ParseDoomColormaps(const std::unique_ptr<core::WadImporter>& importer) {
    std::vector<core::Colormap> colormaps;

    auto data_opt = importer->FindEntryData("COLORMAP");
    if (!data_opt) {
        return colormaps;
    }

    const u8* data = *data_opt;
    u32 data_offset = 0;

    // Pull out the palettes.
    colormaps.resize(kNumDoomColormaps);
    for (i32 i = 0; i < kNumDoomColormaps; i++) {
        memcpy(colormaps[i].map, data + data_offset, 256);
        data_offset += 256;
    }

    return colormaps;
}

}  // namespace doom