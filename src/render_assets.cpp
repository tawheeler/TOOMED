#include "render_assets.hpp"

#include <cstring>

namespace core {

std::optional<u16> FindTextureIdForDoomTextureName(u8 name[8], const RenderAssets& render_assets) {
    for (u16 texture_id = 0; texture_id < render_assets.patches.size(); texture_id++) {
        const doom::Patch& patch = render_assets.patches[texture_id];
        if (strncmp((const char*)name, patch.name.c_str(), 8) == 0) {
            return texture_id;
        }
    }

    return std::nullopt;
}

std::optional<u16> FindFlatIdForDoomFlatName(u8 name[8], const RenderAssets& render_assets) {
    for (u16 flat_id = 0; flat_id < render_assets.flats.size(); flat_id++) {
        const doom::Flat& flat = render_assets.flats[flat_id];
        if (strncmp((const char*)name, flat.name.c_str(), 8) == 0) {
            return flat_id;
        }
    }

    return std::nullopt;
}

}  // namespace core