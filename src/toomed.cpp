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

#define SCREEN_SIZE_X 1280
#define SCREEN_SIZE_Y 720

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

common::Vec2f GlobalToCamera(common::Vec2f g, common::Vec2f camera_pos, f32 camera_zoom) {
    common::Vec2f c = (g - camera_pos) * camera_zoom;
    return common::Vec2f(SCREEN_SIZE_X / 2 + c.x, SCREEN_SIZE_Y / 2 - c.y);
}

common::Vec2f CameraToGlobal(common::Vec2f c, common::Vec2f camera_pos, f32 camera_zoom) {
    common::Vec2f g_offset = common::Vec2f(c.x - SCREEN_SIZE_X / 2, SCREEN_SIZE_Y / 2 - c.y);
    return g_offset / camera_zoom + camera_pos;
}

int main() {
    SDL_version ver;
    SDL_GetVersion(&ver);
    fprintf(stdout, "Running with SDL2 version %d.%d.%d\n", ver.major, ver.minor, ver.patch);

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
    map.vertices.emplace_back(1.0, 1.0);
    map.vertices.emplace_back(7.0, 1.0);
    map.vertices.emplace_back(7.0, 3.0);
    map.vertices.emplace_back(6.0, 3.0);
    map.vertices.emplace_back(9.0, 3.0);
    map.vertices.emplace_back(9.0, 1.0);
    map.vertices.emplace_back(12.0, 1.0);
    map.vertices.emplace_back(12.0, 7.0);
    map.vertices.emplace_back(9.0, 6.0);
    map.vertices.emplace_back(9.0, 5.0);
    map.vertices.emplace_back(6.0, 5.0);
    map.vertices.emplace_back(6.0, 6.0);
    map.vertices.emplace_back(7.0, 6.0);
    map.vertices.emplace_back(7.0, 7.0);
    map.vertices.emplace_back(1.0, 7.0);
    map.vertices.emplace_back(0.0, 0.0);
    map.vertices.emplace_back(13.0, 0.0);
    map.vertices.emplace_back(13.0, 8.0);
    map.vertices.emplace_back(0.0, 8.0);

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
                   "Failed to constrain edge %d -> %d\n", (int)(side_info.a_ind),
                   (int)(side_info.b_ind));
        }
    }

    // Camera parameters
    common::Vec2f camera_pos = {2.0, 2.0};
    f32 camera_zoom = 50.0;

    // Mouse press params
    bool mouse_is_pressed = false;
    common::Vec2f mouse_pos = {0.0, 0.0};
    common::Vec2f mouse_click_pos = {0.0, 0.0};
    common::Vec2f camera_pos_at_mouse_click = {0.0, 0.0};
    int selected_vertex_index = -1;
    core::QuarterEdge* qe_mouse_face = map.mesh.GetEnclosingTriangle(mouse_pos);

    bool continue_running = true;
    while (continue_running) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                continue_running = false;
                break;
            } else if (event.type == SDL_WINDOWEVENT &&
                       event.window.event == SDL_WINDOWEVENT_CLOSE &&
                       event.window.windowID == SDL_GetWindowID(window)) {
                continue_running = false;
                break;
            } else if (event.type == SDL_MOUSEWHEEL) {
                if (event.wheel.preciseY > 0) {
                    camera_zoom *= 1.1;
                } else {
                    camera_zoom /= 1.1;
                }
            } else if (event.type == SDL_MOUSEBUTTONDOWN) {
                if (!mouse_is_pressed) {
                    // New press
                    mouse_is_pressed = true;
                    mouse_click_pos = CameraToGlobal(common::Vec2f(event.button.x, event.button.y),
                                                     camera_pos, camera_zoom);
                    camera_pos_at_mouse_click = camera_pos;

                    // Check for a selected vertex near the current mouse position.
                    f32 min_dist = 999.0;
                    core::QuarterEdge* qe_a = map.mesh.GetTriangleQuarterEdge1(qe_mouse_face);
                    core::QuarterEdge* qe_b = map.mesh.GetTriangleQuarterEdge2(qe_mouse_face);
                    core::QuarterEdge* qe_c = map.mesh.GetTriangleQuarterEdge3(qe_mouse_face);
                    const common::Vec2f& a = qe_a->vertex->vertex;
                    const common::Vec2f& b = qe_b->vertex->vertex;
                    const common::Vec2f& c = qe_c->vertex->vertex;

                    f32 dist_a = common::Norm(a - mouse_click_pos);
                    f32 dist_b = common::Norm(b - mouse_click_pos);
                    f32 dist_c = common::Norm(c - mouse_click_pos);

                    selected_vertex_index = -1;
                    if (dist_a < min_dist) {
                        selected_vertex_index = qe_a->vertex->index;
                        min_dist = dist_a;
                    }
                    if (dist_b < min_dist) {
                        selected_vertex_index = qe_b->vertex->index;
                        min_dist = dist_b;
                    }
                    if (dist_c < min_dist) {
                        selected_vertex_index = qe_c->vertex->index;
                        min_dist = dist_c;
                    }
                }
            } else if (event.type == SDL_MOUSEBUTTONUP) {
                if (mouse_is_pressed) {
                    // New release
                    mouse_is_pressed = false;
                }
            } else if (event.type == SDL_MOUSEMOTION) {
                // Move the mouse
                mouse_pos = CameraToGlobal(common::Vec2f(event.motion.x, event.motion.y),
                                           camera_pos, camera_zoom);
                qe_mouse_face = map.mesh.GetEnclosingTriangle(mouse_pos, qe_mouse_face);

                if (mouse_is_pressed) {
                    // Move the camera
                    camera_pos = camera_pos_at_mouse_click + mouse_click_pos - mouse_pos;
                }
            }
        }

        // Clear screen
        SDL_SetRenderDrawColor(renderer, 0x41, 0x41, 0x41, 0xFF);
        SDL_RenderClear(renderer);

        // Render the texture to the screen.
        // SDL_UpdateTexture(texture, NULL, pixels, SCREEN_SIZE_X * 4);
        // SDL_RenderCopyEx(renderer, texture, NULL, NULL, 0.0, NULL, SDL_FLIP_VERTICAL);

        {
            // Fill enclosing triangle
            const common::Vec2f& a = map.mesh.GetTriangleVertex1(qe_mouse_face);
            const common::Vec2f& b = map.mesh.GetTriangleVertex2(qe_mouse_face);
            const common::Vec2f& c = map.mesh.GetTriangleVertex3(qe_mouse_face);

            auto a_cam = GlobalToCamera(a, camera_pos, camera_zoom);
            auto b_cam = GlobalToCamera(b, camera_pos, camera_zoom);
            auto c_cam = GlobalToCamera(c, camera_pos, camera_zoom);

            SDL_Vertex triangle[3];
            triangle[0] = {
                SDL_FPoint{a_cam.x, a_cam.y},
                SDL_Color{55, 55, 55, 255},
                SDL_FPoint{0},
            };
            triangle[1] = {
                SDL_FPoint{b_cam.x, b_cam.y},
                SDL_Color{55, 55, 55, 255},
                SDL_FPoint{0},
            };
            triangle[2] = {
                SDL_FPoint{c_cam.x, c_cam.y},
                SDL_Color{55, 55, 55, 255},
                SDL_FPoint{0},
            };
            SDL_RenderGeometry(renderer, nullptr, triangle, 3, nullptr, 0);
        }

        {  // Render the mesh
            SDL_SetRenderDrawColor(renderer, 0xFF, 0x48, 0xCF, 0xFF);

            for (size_t qe_index = 0; qe_index < map.mesh.NumQuarterEdges(); qe_index++) {
                core::QuarterEdge* qe = map.mesh.GetQuarterEdge(qe_index);

                if (core::IsPrimalEdge(*qe) && !map.mesh.IsBoundaryVertex(qe->vertex)) {
                    // Get its opposite side.
                    core::QuarterEdge* qe_sym = map.mesh.Sym(qe);

                    const core::VertexData* a = qe->vertex;
                    const core::VertexData* b = qe_sym->vertex;
                    if (a > b && !map.mesh.IsBoundaryVertex(b)) {  // Avoid rendering edges twice
                        auto a_cam = GlobalToCamera(a->vertex, camera_pos, camera_zoom);
                        auto b_cam = GlobalToCamera(b->vertex, camera_pos, camera_zoom);
                        SDL_RenderDrawLine(renderer, (int)(a_cam.x), (int)(a_cam.y), (int)(b_cam.x),
                                           (int)(b_cam.y));
                    }
                }
            }
        }

        {  // Render all side_infos
            SDL_SetRenderDrawColor(renderer, 0x80, 0x80, 0x80, 0xFF);

            for (const auto& side_info : map.side_infos) {
                common::Vec2f a = map.vertices[side_info.a_ind];
                common::Vec2f b = map.vertices[side_info.b_ind];

                auto a_cam = GlobalToCamera(a, camera_pos, camera_zoom);
                auto b_cam = GlobalToCamera(b, camera_pos, camera_zoom);
                SDL_RenderDrawLine(renderer, (int)(a_cam.x), (int)(a_cam.y), (int)(b_cam.x),
                                   (int)(b_cam.y));
            }
        }

        {  // Render all vertices

            SDL_SetRenderDrawColor(renderer, 0x90, 0x90, 0x90, 0xFF);

            for (const auto& v : map.vertices) {
                auto v_cam = GlobalToCamera(v, camera_pos, camera_zoom);

                SDL_Rect rect;

                // Outline with darker color
                SDL_SetRenderDrawColor(renderer, 0x41, 0x41, 0x41, 0xFF);
                rect.x = (int)(v_cam.x - 2);
                rect.y = (int)(v_cam.y - 2);
                rect.h = 5;
                rect.w = 5;
                SDL_RenderFillRect(renderer, &rect);

                // Fill with lighter color
                SDL_SetRenderDrawColor(renderer, 0x90, 0x90, 0x90, 0xFF);
                rect.x = (int)(v_cam.x - 1);
                rect.y = (int)(v_cam.y - 1);
                rect.h = 3;
                rect.w = 3;
                SDL_RenderFillRect(renderer, &rect);
            }
        }

        if (selected_vertex_index != -1) {
            // Render our selected vertex (Offset by 3 to skip the enclosing vertices).
            const auto& v = map.vertices[selected_vertex_index - 3];

            auto v_cam = GlobalToCamera(v, camera_pos, camera_zoom);

            SDL_Rect rect;

            // Outline with darker color
            SDL_SetRenderDrawColor(renderer, 0x41, 0x41, 0x41, 0xFF);
            rect.x = (int)(v_cam.x - 3);
            rect.y = (int)(v_cam.y - 3);
            rect.h = 7;
            rect.w = 7;
            SDL_RenderFillRect(renderer, &rect);

            // Fill with lighter color
            SDL_SetRenderDrawColor(renderer, 0xFF, 0xFF, 0xFF, 0xFF);
            rect.x = (int)(v_cam.x - 2);
            rect.y = (int)(v_cam.y - 2);
            rect.h = 5;
            rect.w = 5;
            SDL_RenderFillRect(renderer, &rect);

            // Render the mouse pos
            v_cam = GlobalToCamera(mouse_click_pos, camera_pos, camera_zoom);

            // Outline with darker color
            SDL_SetRenderDrawColor(renderer, 0x41, 0x41, 0x41, 0xFF);
            rect.x = (int)(v_cam.x - 3);
            rect.y = (int)(v_cam.y - 3);
            rect.h = 7;
            rect.w = 7;
            SDL_RenderFillRect(renderer, &rect);

            // Fill with lighter color
            SDL_SetRenderDrawColor(renderer, 0xFF, 0xFF, 0xFF, 0xFF);
            rect.x = (int)(v_cam.x - 2);
            rect.y = (int)(v_cam.y - 2);
            rect.h = 5;
            rect.w = 5;
            SDL_RenderFillRect(renderer, &rect);
        }

        // SDL_RENDERER_PRESENTVSYNC means this is syncronized with the monitor refresh rate.
        // (30Hz)
        SDL_RenderPresent(renderer);
    }

    // Destroy our window
    SDL_DestroyWindow(window);
}