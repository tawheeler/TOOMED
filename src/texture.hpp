#pragma once

#include <string>
#include <vector>

#include "typedefs.hpp"
#include "wad_importer.hpp"

namespace doom {

class Texture {
  public:
    Texture(const std::string& name, u16 size_x, u16 size_y);

    u16 SizeX() const { return size_x_; }
    u16 SizeY() const { return size_y_; }

    void SetPixel(u32 idx, u8 palette_index);
    void SetPixel(u16 x, u16 y, u8 palette_index);

    u8 GetPixel(u32 idx) const;
    u8 GetPixel(u16 x, u16 y) const;

  private:
    std::string name_;
    u16 size_x_;            // number of bytes wide
    u16 size_y_;            // number of bytes tall
    std::vector<u8> data_;  // indices into our palette, indexed by data[x*size_y + y]
    // The texture data is stored top-to-bottom, left to right:
    //  o ------> x
    //  |
    //  |
    //  v y
};

std::vector<Texture> ParseDoomTextures(const std::unique_ptr<core::WadImporter>& importer);

}  // namespace doom