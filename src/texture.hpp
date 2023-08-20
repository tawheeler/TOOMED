#pragma once

#include <string>
#include <vector>

#include "assets_exporter.hpp"
#include "typedefs.hpp"
#include "wad_importer.hpp"

namespace doom {

// Based on the Patch used by DOOM.
// https://doomwiki.org/wiki/Picture_format
// Patches are columnar images where each pixel column is represented by a list of posts.
// This allows for skipping transparent tiles and for avoiding duplicate columns.
// Note that this format does not use Posts exactly, but skips the unneeded padding bytes.
struct Patch {
    std::string name;
    u16 size_x;    // number of bytes wide
    u16 size_y;    // number of bytes tall
    u16 origin_x;  // the number of pixels from the left edge that this patch's origin is
    u16 origin_y;  // the number of pixels from the top edge that this patch's origin is

    // Array of indices into `post_data` for each pixel column (i.e. of length `size_x`).
    // Note that this is different than the column_offsets used in the DOOM format, which
    // store the number of bytes past the start of the mp, and thus include the Patch header.
    std::vector<u32> column_offsets;

    // Each post consists of:
    //     u8   y_delta: The number of transparent pixels to skip prior to this post
    //     u8   length: The number of pixels in this post
    //     u8[] palette_indices: The array of palette indices (of length `length`)
    // Posts are terminated by 0xFF
    std::vector<u8> post_data;
};

class Texture {
  public:
    Texture(const std::string& name, u16 size_x, u16 size_y);

    u16 SizeX() const { return size_x_; }
    u16 SizeY() const { return size_y_; }

    void SetPixel(u32 idx, u8 palette_index);
    void SetPixel(u16 x, u16 y, u8 palette_index);

    u8 GetPixel(u32 idx) const;
    u8 GetPixel(u16 x, u16 y) const;

    u8 IsTransparent(u32 idx) const;
    u8 IsTransparent(u16 x, u16 y) const;

    Patch ToPatch();

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
    std::vector<bool> transparent_;  // Track whether pixels have been written to
};

std::vector<Patch> ParseDoomTextures(const std::unique_ptr<core::WadImporter>& importer);

// The name of the exported asset entry.
const std::string kAssetEntryPatches = "patches";
core::AssetsExporterEntry ExportPatches(const std::vector<Patch>& patches);
bool ImportPatches(std::vector<doom::Patch>* patches, const core::AssetsExporter& exporter);

// A flat in DOOM is a 64x64 set of pixel indices.
struct Flat {
    std::string name;
    u8 data[4096];
};

std::vector<Flat> ParseDoomFlats(const std::unique_ptr<core::WadImporter>& importer);

// The name of the exported asset entry.
const std::string kAssetEntryFlats = "flats";
core::AssetsExporterEntry ExportFlats(const std::vector<Flat>& flats);
bool ImportFlats(std::vector<doom::Flat>* flats, const core::AssetsExporter& exporter);

}  // namespace doom