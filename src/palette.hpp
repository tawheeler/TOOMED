#pragma once

#include <vector>

#include "assets_exporter.hpp"
#include "typedefs.hpp"
#include "wad_importer.hpp"

namespace core {

// A DOOM color palette.
// Each color entry is 3 u8's of RGB.
struct Palette {
    u8 rgbs[768];
};

// The name of the exported asset entry.
const std::string kAssetEntryPalettes = "palettes";

// Export the palettes.
AssetsExporterEntry ExportPalettes(const std::vector<Palette>& palettes);

// The DOOM colormaps.
// Each colormap maps palette indices (u8's) to new palette indices.
// DOOM contains 34 colormaps, the first 32 of which have decreasing brightness,
// the next one being used for invulnerability, and the last for nothing.
struct Colormap {
    u8 map[256];
};

// The name of the exported asset entry.
const std::string kAssetEntryColormaps = "colormaps";

// Export the colormaps.
AssetsExporterEntry ExportColormaps(const std::vector<Colormap>& colormaps);

}  // namespace core

// ------------------------------------------------------------------------------------------------
namespace doom {

// Parse the DOOM PLAYPAL. https://doomwiki.org/wiki/PLAYPAL
constexpr i32 kNumDoomPalettes = 14;
std::vector<core::Palette> ParseDoomPalettes(const std::unique_ptr<core::WadImporter>& importer);

// Parse the DOOM COLORMAPS. https://doomwiki.org/wiki/COLORMAP
constexpr i32 kNumDoomColormaps = 34;
std::vector<core::Colormap> ParseDoomColormaps(const std::unique_ptr<core::WadImporter>& importer);

}  // namespace doom