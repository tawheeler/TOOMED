#include "palette.hpp"

#include <cmath>
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
bool ImportPalettes(std::vector<Palette>* palettes, const AssetsExporter& exporter) {
    const AssetsExporterEntry* entry = exporter.FindEntry(kAssetEntryPalettes);
    if (entry == nullptr) {
        std::cout << "Failed to find entry " << kAssetEntryPalettes << std::endl;
        return false;
    }

    u32 offset = 0;
    const u8* data = entry->data.data();

    palettes->clear();

    u32 n_palettes = *(u32*)(data + offset);
    offset += sizeof(u32);

    palettes->resize(n_palettes);
    for (u32 i = 0; i < n_palettes; i++) {
        memcpy(palettes->at(i).rgbs, data + offset, 256 * 3);
        offset += 256 * 3;
    }

    return true;
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

// ------------------------------------------------------------------------------------------------
bool ImportColormaps(std::vector<Colormap>* colormaps, const AssetsExporter& exporter) {
    const AssetsExporterEntry* entry = exporter.FindEntry(kAssetEntryColormaps);
    if (entry == nullptr) {
        std::cout << "Failed to find entry " << kAssetEntryColormaps << std::endl;
        return false;
    }

    u32 offset = 0;
    const u8* data = entry->data.data();

    colormaps->clear();

    u32 n_colormaps = *(u32*)(data + offset);
    offset += sizeof(u32);

    colormaps->resize(n_colormaps);
    for (u32 i = 0; i < n_colormaps; i++) {
        memcpy(colormaps->at(i).map, data + offset, sizeof(colormaps->at(i).map));
        offset += sizeof(colormaps->at(i).map);
    }

    return true;
}

f32 PivotRGB(f32 c) {
    if (c > 0.04045) {
        return std::pow((c + 0.055) / 1.055, 2.4) * 100;
    } else {
        return c / 12.92 * 100;
    }
}

f32 PivotXYZ(f32 c) {
    f32 i = std::cbrt(c);
    if (c > 0.008856) {
        return i;
    } else {
        return 7.787 * c + 16.0 / 116.0;
    }
}

LAB RGB2LAB(u8 r, u8 g, u8 b) {
    f32 rp = PivotRGB(r / 255.0);
    f32 gp = PivotRGB(b / 255.0);
    f32 bp = PivotRGB(g / 255.0);

    f32 x = 0.4124 * rp + 0.3576 * gp + 0.1805 * bp;
    f32 y = 0.2126 * rp + 0.7152 * gp + 0.0722 * bp;
    f32 z = 0.0193 * rp + 0.1192 * gp + 0.9505 * bp;

    // Observer = deg, Illuminant = D65
    f32 xp = PivotXYZ(x / 95.047);
    f32 yp = PivotXYZ(y / 100.00);
    f32 zp = PivotXYZ(z / 108.883);

    LAB retval;
    retval.l = 116.0 * yp - 16.0;
    retval.a = 500.0 * (xp - yp);
    retval.b = 200.0 * (yp - zp);
    return retval;
}

// https://en.wikipedia.org/wiki/Color_difference#CIE94
f32 CalcColorDistanceCie94(LAB A, LAB B) {
    constexpr f32 kl = 1.0;
    constexpr f32 kc = 1.0;
    constexpr f32 kH = 1.0;
    constexpr f32 K1 = 0.045;
    constexpr f32 K2 = 0.015;

    f32 C1 = std::hypot(A.a, A.b);
    f32 C2 = std::hypot(B.a, B.b);
    f32 dC = C1 - C2;
    f32 dL = A.l - B.l;
    f32 da = A.a - B.a;
    f32 db = A.b - B.b;
    f32 dH = std::sqrt(std::max(da * da + db * db - dC * dC, 0.0f));  // TODO
    f32 SL = 1.0;
    f32 SC = 1.0 + K1 * C1;
    f32 SH = 1.0 + K2 * C1;
    return sqrt(std::pow(dL / (kl * SL), 2) + std::pow(dC / (kc * SC), 2) +
                std::pow(dH / (kH * SH), 2));
}

u8 GetClosestColorInPalette(u32 abgr, std::vector<LAB> palette) {
    u8 r = abgr & 0xFF;
    u8 g = (abgr >> 8) & 0xFF;
    u8 b = (abgr >> 16) & 0xFF;
    LAB lab = RGB2LAB(r, g, b);

    f32 min_dist = INFINITY;
    u8 best_index = 0;
    for (usize i = 0; i < palette.size(); i++) {
        f32 dist = CalcColorDistanceCie94(lab, palette[i]);
        if (dist < min_dist) {
            min_dist = dist;
            best_index = (u8)i;
        }
    }

    return best_index;
}

std::vector<LAB> GetLabPalette(const Palette& palette) {
    std::vector<LAB> retval;
    retval.reserve(256);

    u32 offset = 0;
    for (u32 i = 0; i < 256; i++) {
        u8 r = *(palette.rgbs + offset);
        offset++;
        u8 g = *(palette.rgbs + offset);
        offset++;
        u8 b = *(palette.rgbs + offset);
        offset++;

        retval.push_back(RGB2LAB(r, g, b));
    }

    return retval;
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