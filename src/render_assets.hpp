#pragma once

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

}  // namespace core