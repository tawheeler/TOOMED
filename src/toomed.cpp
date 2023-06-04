#include <SDL2/SDL.h>

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "delaunay_mesh.hpp"
#include "geometry_utils.hpp"
#include "math_utils.hpp"

#define SCREEN_SIZE_X 640
#define SCREEN_SIZE_Y 360

#define MESH_BOUNDING_RADIUS 1000.0f
#define MESH_MIN_DIST_TO_VERTEX 0.1f
#define MESH_MIN_DIST_TO_EDGE 0.1f

#define ASSERT(_e, ...)               \
    if (!(_e)) {                      \
        fprintf(stderr, __VA_ARGS__); \
        exit(1);                      \
    }

using namespace std;

typedef float f32;
typedef double f64;
typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;
typedef size_t usize;

// The information associated with one side of an edge between vertices in the map.
// If this represents the directed edge A->B, then it describes the edge viewed on the right side of
// A->B.
struct SideInfo {
    u16 flags;
    u16 texture_id;
    i16 x_offset;  // Texture x offset
    i16 y_offset;  // Texture y offset
    usize a_ind;   // The index of vertex A (the source vertex)
    usize b_ind;   // The index of vertex B (the dest vertex)
};

// Represents our game map
struct Map {
    // The editable map vertices.
    // This will exactly match the vertices in the DelaunayMesh.
    // (We error if adding any vertex to the DelaunayMesh fails (due to coincidence, for example).)
    std::vector<common::Vec2f> vertices;

    // All of the map-related side information.
    // If we have side information for an edge A -> B, then that edge must end up in the mesh.
    // Edges without side information may exist in the mesh. Such edges are assumed transparent.
    std::vector<SideInfo> side_infos;

    // Map <a_ind, b_ind> to index in side_infos.
    std::map<std::tuple<usize, usize>, usize> side_to_info;

    // The map geometry
    core::DelaunayMesh mesh =
        core::DelaunayMesh(MESH_BOUNDING_RADIUS, MESH_MIN_DIST_TO_VERTEX, MESH_MIN_DIST_TO_EDGE);

    // face data (is solid, height, is door, etc.)

    usize AddEdge(int a_ind, int b_ind) {
        SideInfo side_info;
        side_info.a_ind = a_ind;
        side_info.b_ind = b_ind;
        side_infos.emplace_back(side_info);
        usize edge_index = side_infos.size();
        side_to_info[std::make_pair(a_ind, b_ind)] = edge_index;
        return edge_index;
    }
};

int main() {
    // Initialize SDL
    ASSERT(SDL_Init(SDL_INIT_VIDEO) == 0, "SDL initialization failed: %s\n", SDL_GetError());

    // Create a window
    SDL_Window* window = SDL_CreateWindow("TOOM EDITOR", SDL_WINDOWPOS_CENTERED_DISPLAY(1),
                                          SDL_WINDOWPOS_CENTERED_DISPLAY(1), SCREEN_SIZE_X,
                                          SCREEN_SIZE_Y, SDL_WINDOW_ALLOW_HIGHDPI);
    ASSERT(window, "Error creating SDL window: %s\n", SDL_GetError());

    // Create a renderer
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC);
    ASSERT(renderer, "Error creating SDL renderer: %s\n", SDL_GetError());

    // Create a texture
    SDL_Texture* texture =
        SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STREAMING,
                          SCREEN_SIZE_X, SCREEN_SIZE_Y);
    ASSERT(texture, "Error creating SDL texture: %s\n", SDL_GetError());

    // Create our map
    Map map;
    map.vertices.emplace_back(50.0, 50.0);
    map.vertices.emplace_back(100.0, 50.0);
    map.vertices.emplace_back(100.0, 100.0);
    map.vertices.emplace_back(50.0, 100.0);
    map.AddEdge(0, 1);
    map.AddEdge(1, 2);
    map.AddEdge(2, 3);
    map.AddEdge(3, 0);

    // Construct the mesh
    {
        // Add all vertices
        for (const auto& v : map.vertices) {
            int new_vertex = map.mesh.AddDelaunayVertex(v);
            ASSERT(new_vertex != core::kInvalidIndex, "Failed to add vertex into mesh\n");
        }

        // Constrain all edges
        for (const auto& side_info : map.side_infos) {
            ASSERT(map.mesh.ConstrainEdge(side_info.a_ind, side_info.b_ind),
                   "Failed to constraint edge %d -> %d\n", (int)(side_info.a_ind),
                   (int)(side_info.b_ind));
        }
    }

    bool continue_running = true;
    while (continue_running) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                continue_running = false;
                break;
            }
        }

        // Clear screen
        SDL_SetRenderDrawColor(renderer, 0x41, 0x41, 0x41, 0xFF);
        SDL_RenderClear(renderer);

        // Render the texture to the screen.
        // SDL_UpdateTexture(texture, NULL, pixels, SCREEN_SIZE_X * 4);
        // SDL_RenderCopyEx(renderer, texture, NULL, NULL, 0.0, NULL, SDL_FLIP_VERTICAL);

        {  // Render the mesh
            SDL_SetRenderDrawColor(renderer, 0xFF, 0x48, 0xCF, 0xFF);

            for (size_t qe_index = 0; qe_index < map.mesh.NumQuarterEdges(); qe_index++) {
                core::QuarterEdge* qe = map.mesh.GetQuarterEdge(qe_index);

                if (core::IsPrimalEdge(*qe) && !map.mesh.IsBoundaryVertex(qe->vertex)) {
                    // Get its opposite side.
                    core::QuarterEdge* qe_sym = map.mesh.Sym(qe);

                    const common::Vec2f* a = qe->vertex;
                    const common::Vec2f* b = qe_sym->vertex;
                    if (a > b && !map.mesh.IsBoundaryVertex(b)) {  // Avoid rendering edges twice
                        int ax = (int)(a->x);
                        int ay = SCREEN_SIZE_Y - (int)(a->y);
                        int bx = (int)(b->x);
                        int by = SCREEN_SIZE_Y - (int)(b->y);
                        SDL_RenderDrawLine(renderer, ax, ay, bx, by);
                    }
                }
            }
        }

        {  // Render all side_infos
            SDL_SetRenderDrawColor(renderer, 0x80, 0x80, 0x80, 0xFF);

            for (const auto& side_info : map.side_infos) {
                common::Vec2f a = map.vertices[side_info.a_ind];
                common::Vec2f b = map.vertices[side_info.b_ind];

                int ax = (int)(a.x);
                int ay = SCREEN_SIZE_Y - (int)(a.y);
                int bx = (int)(b.x);
                int by = SCREEN_SIZE_Y - (int)(b.y);
                SDL_RenderDrawLine(renderer, ax, ay, bx, by);
            }
        }

        {  // Render all vertices

            SDL_SetRenderDrawColor(renderer, 0x90, 0x90, 0x90, 0xFF);

            for (const auto& v : map.vertices) {
                int x = (int)(v.x);
                int y = SCREEN_SIZE_Y - (int)(v.y);

                SDL_Rect rect;

                // Outline with darker color
                SDL_SetRenderDrawColor(renderer, 0x41, 0x41, 0x41, 0xFF);
                rect.x = (int)(x - 2);
                rect.y = (int)(y - 2);
                rect.h = 5;
                rect.w = 5;
                SDL_RenderFillRect(renderer, &rect);

                // Fill with lighter color
                SDL_SetRenderDrawColor(renderer, 0x90, 0x90, 0x90, 0xFF);
                rect.x = (int)(x - 1);
                rect.y = (int)(y - 1);
                rect.h = 3;
                rect.w = 3;
                SDL_RenderFillRect(renderer, &rect);
            }
        }

        // SDL_RENDERER_PRESENTVSYNC means this is syncronized with the monitor refresh rate. (30Hz)
        SDL_RenderPresent(renderer);
    }

    // Destroy our window
    SDL_DestroyWindow(window);
}