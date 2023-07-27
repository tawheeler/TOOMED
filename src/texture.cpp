#include "texture.hpp"

namespace doom {

// ------------------------------------------------------------------------------------------------
Texture::Texture(const std::string& name, u16 size_x, u16 size_y) :
    name_(name), size_x_(size_x), size_y_(size_y) {
    data_.resize(size_x_ * size_y_);
}

// ------------------------------------------------------------------------------------------------
void Texture::SetPixel(u32 idx, u8 palette_index) { data_[idx] = palette_index; }
void Texture::SetPixel(u16 x, u16 y, u8 palette_index) {
    u32 idx = (u32)(size_y_)*x + y;
    data_[idx] = palette_index;
}

// ------------------------------------------------------------------------------------------------
u8 Texture::GetPixel(u32 idx) const { return data_[idx]; }
u8 Texture::GetPixel(u16 x, u16 y) const {
    u32 idx = (u32)(size_y_)*x + y;
    return data_[idx];
}

// ------------------------------------------------------------------------------------------------
std::vector<Texture> ParseDoomTextures(const std::unique_ptr<core::WadImporter>& importer) {
    std::vector<Texture> retval;

    auto texture_data_opt = importer->FindEntryData("TEXTURE1");
    auto pnames_data_opt = importer->FindEntryData("PNAMES");  // all the names for wall patches

    if (!texture_data_opt || !pnames_data_opt) {
        return retval;
    }

    // First, parse PNAMES (https://doomwiki.org/wiki/PNAMES)
    const u8* pnames_data = *pnames_data_opt;
    i32 n_map_patches = *(i32*)(pnames_data);  // not really used

    // Parse the texture data (https://doomwiki.org/wiki/TEXTURE1_and_TEXTURE2)
    const u8* texture_data = *texture_data_opt;

    i32 n_textures = *(i32*)(texture_data);

    retval.reserve(n_textures);
    for (i32 i_texture = 0; i_texture < n_textures; i_texture++) {
        // Offset of the map texture in `data`
        i32 byte_offset = *(i32*)(texture_data + 4 * (i_texture + 1));

        char* name = (char*)(texture_data + byte_offset);  // 8 bytes long
        byte_offset += 8;

        // skip masked
        byte_offset += 4;

        // Total width of the map texture
        u16 size_x = *(u16*)(texture_data + byte_offset);
        byte_offset += sizeof(u16);

        // Total height of the map texture
        u16 size_y = *(u16*)(texture_data + byte_offset);
        byte_offset += sizeof(u16);

        // skip columndirectory
        byte_offset += 4;

        // Number of patches that make up this texture
        u16 n_patches = *(u16*)(texture_data + byte_offset);
        byte_offset += sizeof(u16);

        // Create a new texture
        Texture texture(std::string(name, 8), size_x, size_y);

        // Array of map patches
        struct Subpatch {
            i16 origin_x;  // 	A short int defining the horizontal offset of the patch relative to
                           // the upper-left of the texture.
            i16 origin_y;  // A short int defining the vertical offset of the patch relative to the
                           // upper-left of the texture.
            i16 i_patch;   // A short int defining the patch number (as listed in PNAMES) to draw.
            i16 step_dir;  // A short int possibly intended to define if the patch was to be drawn
                           // normally or mirrored. Unused.
            i16 i_colormap;  // A short int possibly intended to define a special colormap to draw
                             // the patch with, like a brightmap. Unused.
        };
        for (u16 i_patch = 0; i_patch < n_patches; i_patch++) {
            Subpatch subpatch = *(Subpatch*)(texture_data + byte_offset);
            byte_offset += sizeof(Subpatch);

            // Get the patch name (max 8 chars)
            std::string subpatch_name((char*)(pnames_data + 4 + 8 * subpatch.i_patch), 8);

            auto patch_opt = importer->FindEntryData(subpatch_name);
            if (!patch_opt) {
                continue;
            }

            // Render the patch
            const u8* patch_data = *patch_opt;
            u16 patch_size_x = *(u16*)(patch_data);
            u16 patch_size_y = *(u16*)(patch_data + 2);
            // i16 offset_x = *(u16*)(patch_data + 4); // unused
            // i16 offset_y = *(u16*)(patch_data + 6); // unused
            for (u16 patch_x = 0; patch_x < patch_size_x; patch_x++) {
                i32 texture_x = patch_x + subpatch.origin_x;
                if (texture_x < 0 || texture_x >= texture.SizeX()) {
                    continue;
                }

                u32 column_offset = *(u32*)(patch_data + 8 + 4 * patch_x);

                u8 y_from_top = *(u8*)(patch_data + column_offset);
                while (y_from_top != 0xFF) {
                    i32 texture_y = y_from_top + subpatch.origin_y;
                    u8 n_pix_in_post = *(u8*)(patch_data + column_offset + 1);
                    for (u8 i_pix = 0; i_pix < n_pix_in_post; i_pix++) {
                        u8 palette_index = *(u8*)(patch_data + column_offset + 3 + i_pix);
                        if (texture_y >= 0 && texture_y < texture.SizeY()) {
                            texture.SetPixel(texture_x, texture_y, palette_index);
                        }
                        texture_y += 1;
                    }

                    column_offset += 4 + n_pix_in_post;
                    y_from_top = *(u8*)(patch_data + column_offset);
                }
            }
        }

        // Store the texture
        retval.push_back(texture);
    }

    return retval;
}

}  // namespace doom