#include "game_map.hpp"

namespace core {

int MapVertexIndexToMeshVertexIndex(int ind) { return ind + 3; }
int MeshVertexIndexToMapVertexIndex(int ind) { return ind - 3; }

}