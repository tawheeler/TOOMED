#include <SDL2/SDL.h>

#include <cstdio>
#include <iostream>
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

// Represents our game map
struct Map {
    // The map geometry
    core::DelaunayMesh mesh =
        core::DelaunayMesh(MESH_BOUNDING_RADIUS, MESH_MIN_DIST_TO_VERTEX, MESH_MIN_DIST_TO_EDGE);

    // edge data (texture, is opaque, etc.)
    // face data (is solid, height, is door, etc.)
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
    map.mesh.AddDelaunayVertex({0.0, 0.0});
    map.mesh.AddDelaunayVertex({0.0, 50.0});
    map.mesh.AddDelaunayVertex({50.0, 50.0});
    map.mesh.AddDelaunayVertex({50.0, 0.0});
    map.mesh.AddDelaunayVertex({100.0, 100.0});

    bool continue_running = true;
    while (continue_running) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                continue_running = false;
                break;
            }
        }

        // Render the texture to the screen.
        // SDL_UpdateTexture(texture, NULL, pixels, SCREEN_SIZE_X * 4);
        SDL_RenderCopyEx(renderer, texture, NULL, NULL, 0.0, NULL, SDL_FLIP_VERTICAL);

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
                        int ax = a->x;
                        int ay = SCREEN_SIZE_Y - a->y;
                        int bx = b->x;
                        int by = SCREEN_SIZE_Y - b->y;
                        SDL_RenderDrawLine(renderer, ax, ay, bx, by);
                    }
                }
            }
        }

        // SDL_RENDERER_PRESENTVSYNC means this is syncronized with the monitor refresh rate. (30Hz)
        SDL_RenderPresent(renderer);
    }

    // Destroy our window
    SDL_DestroyWindow(window);
}