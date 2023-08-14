#pragma once

#include <optional>
#include <vector>

#include "palette.hpp"
#include "texture.hpp"
#include "typedefs.hpp"

namespace core {

// All of the assets used during rendering
struct RenderAssets {
    std::vector<core::Palette> palettes;
    std::vector<core::Colormap> colormaps;
    std::vector<doom::Patch> patches;
};

// Find the id of the texture with the given name (only 8 chars), if it exists.
std::optional<u16> FindTextureIdForDoomTextureName(u8 name[8], const RenderAssets& render_assets);

}  // namespace core