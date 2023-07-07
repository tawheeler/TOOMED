#include "palette.hpp"

#include <iostream>
#include <sstream>

namespace core {

// ------------------------------------------------------------------------------------------------
AssetsExporterEntry ExportColorPalette() {
    AssetsExporterEntry entry;
    entry.SetName(kAssetEntryColorPalette);

    // Serialize the data
    std::stringstream buffer;
    buffer.write(reinterpret_cast<const char*>(&COLOR_PALETTE), sizeof(COLOR_PALETTE));
    entry.SetData(buffer.str());

    return entry;
}

}