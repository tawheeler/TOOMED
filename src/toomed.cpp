#include <SDL2/SDL.h>

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "assets_exporter.hpp"
#include "delaunay_mesh.hpp"
#include "game_map.hpp"
#include "geometry_utils.hpp"
#include "math_utils.hpp"
#include "typedefs.hpp"

#define SCREEN_SIZE_X 1280
#define SCREEN_SIZE_Y 720

#define ASSERT(_e, ...)               \
    if (!(_e)) {                      \
        fprintf(stderr, __VA_ARGS__); \
        exit(1);                      \
    }

common::Vec2f GlobalToCamera(const common::Vec2f& g, const common::Vec2f& camera_pos,
                             f32 camera_zoom) {
    common::Vec2f c = (g - camera_pos) * camera_zoom;
    return common::Vec2f(SCREEN_SIZE_X / 2 + c.x, SCREEN_SIZE_Y / 2 - c.y);
}

common::Vec2f CameraToGlobal(const common::Vec2f& c, const common::Vec2f& camera_pos,
                             f32 camera_zoom) {
    common::Vec2f g_offset = common::Vec2f(c.x - SCREEN_SIZE_X / 2, SCREEN_SIZE_Y / 2 - c.y);
    return g_offset / camera_zoom + camera_pos;
}

void SetColor(SDL_Renderer* renderer, u32 rgba) {
    u8 a = rgba & 0xFF;
    rgba >>= 8;
    u8 b = rgba & 0xFF;
    rgba >>= 8;
    u8 g = rgba & 0xFF;
    rgba >>= 8;
    u8 r = rgba & 0xFF;
    SDL_SetRenderDrawColor(renderer, r, g, b, a);
}

void RenderGrid(SDL_Renderer* renderer, u32 rgba, f32 line_spacing, const common::Vec2f& camera_pos,
                f32 camera_zoom) {
    SetColor(renderer, rgba);

    f32 line_spacing_camera = line_spacing * camera_zoom;
    f32 x_global = camera_pos.x - fmod(camera_pos.x, line_spacing);
    f32 x_camera = SCREEN_SIZE_X / 2 + (x_global - camera_pos.x) * camera_zoom;
    x_camera = fmod(x_camera, line_spacing_camera);
    while (x_camera < SCREEN_SIZE_X) {
        x_camera += line_spacing_camera;
        SDL_RenderDrawLine(renderer, (int)(x_camera), (int)(0), (int)(x_camera),
                           (int)(SCREEN_SIZE_Y));
    }

    f32 y_global = camera_pos.y - fmod(camera_pos.y, line_spacing);
    f32 y_camera = SCREEN_SIZE_Y / 2 - (y_global - camera_pos.y) * camera_zoom;
    y_camera = fmod(y_camera, line_spacing_camera);
    while (y_camera < SCREEN_SIZE_Y) {
        y_camera += line_spacing_camera;
        SDL_RenderDrawLine(renderer, (int)(0), (int)(y_camera), (int)(SCREEN_SIZE_X),
                           (int)(y_camera));
    }
}

// void ImportGameData(core::GameMap* map) {
//     std::cout << "--------------------------------------" << std::endl;
//     std::cout << "Importing game data" << std::endl;

//     core::AssetsExporter exporter;
//     exporter.LoadAssetsFile("../toom/assets/toomed.bin");
//     std::cout << "Num entries:" << exporter.NumEntries() << std::endl;

//     bool succeeded = map->Import(exporter);
//     if (!succeeded) {
//         std::cout << "Failed to load game map! Clearing possibly corrupted map." << std::endl;
//         map->Clear();
//     }

//     std::cout << "DONE" << std::endl;
//     std::cout << "--------------------------------------" << std::endl;
// }

// void ExportGameData(const core::GameMap& map) {
//     std::cout << "--------------------------------------" << std::endl;
//     std::cout << "Exporting game data" << std::endl;

//     core::AssetsExporter exporter;
//     bool succeeded = map.Export(&exporter);
//     if (!succeeded) {
//         std::cout << "Failure while exporting game map! Not written to file." << std::endl;
//         return;
//     }

//     std::cout << "Num entries:" << exporter.NumEntries() << std::endl;

//     exporter.WriteToFile("../toom/assets/toomed.bin");

//     std::cout << "DONE" << std::endl;
//     std::cout << "--------------------------------------" << std::endl;
// }

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
    core::GameMap map;

    // Camera parameters
    common::Vec2f camera_pos = {2.0, 2.0};
    f32 camera_zoom = 50.0;

    // Mouse press params
    bool mouse_is_pressed = false;
    common::Vec2f mouse_pos = {0.0, 0.0};
    common::Vec2f mouse_click_pos = {0.0, 0.0};
    common::Vec2f camera_pos_at_mouse_click = {0.0, 0.0};
    core::VertexIndex selected_vertex_index = {core::kInvalidIndex};
    core::QuarterEdgeIndex selected_edge_index = {core::kInvalidIndex};
    core::QuarterEdgeIndex qe_mouse_face = map.GetMesh().GetEnclosingTriangle(mouse_pos);

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

                    // Check for a selected vertex near the current mouse position
                    selected_vertex_index =
                        map.FindVertexNearPosition(mouse_click_pos, qe_mouse_face);

                    // Check for a selected edge near the current mouse position if we did not
                    // select a vertex
                    selected_edge_index = {core::kInvalidIndex};
                    if (!IsValid(selected_vertex_index)) {
                        selected_edge_index =
                            map.FindEdgeNearPosition(mouse_click_pos, qe_mouse_face);
                    }
                }
            } else if (event.type == SDL_MOUSEBUTTONUP) {
                if (mouse_is_pressed) {
                    // New release
                    mouse_is_pressed = false;

                    if (core::IsValid(selected_vertex_index)) {
                        // Check for a selected vertex near the released mouse position
                        // core::VertexIndex released_vertex_index =
                        //     map.FindVertexNearPosition(mouse_pos, qe_mouse_face);
                        //     if (released_vertex_index.has_value() &&
                        //         released_vertex_index != selected_vertex_index) {
                        //         // Join those edges
                        //         int src = *selected_vertex_index;
                        //         int dst = *released_vertex_index;
                        //         if (!map.HasEdge(src, dst)) {
                        //             map.AddDirectedEdge(src, dst);
                        //         }
                        //     }
                    } else {
                        // We do not have a selected vertex index
                        const bool holding_v = SDL_GetKeyboardState(nullptr)[SDL_SCANCODE_V];
                        const bool did_not_drag = common::Norm(mouse_pos - mouse_click_pos) < 0.2;
                        if (holding_v && did_not_drag) {
                            // Add a new vertex.
                            map.AddVertex(mouse_click_pos);
                        }
                    }
                }
            } else if (event.type == SDL_MOUSEMOTION) {
                // Move the mouse
                mouse_pos = CameraToGlobal(
                    common::Vec2f(event.motion.x, event.motion.y),
                    mouse_is_pressed ? camera_pos_at_mouse_click : camera_pos, camera_zoom);
                // Update the mouse face
                qe_mouse_face = map.GetMesh().GetEnclosingTriangle(mouse_pos, qe_mouse_face);

                // Pan the camera
                if (mouse_is_pressed && !core::IsValid(selected_vertex_index)) {
                    camera_pos = camera_pos_at_mouse_click + mouse_click_pos - mouse_pos;
                }

            } else if (event.type == SDL_KEYDOWN) {
                if (event.key.keysym.sym == SDLK_e) {
                    // ExportGameData(map); // TODO
                } else if (event.key.keysym.sym == SDLK_i) {
                    // ImportGameData(&map); // TODO
                } else if (event.key.keysym.sym == SDLK_DELETE) {
                    // Delete key pressed!
                    if (core::IsValid(selected_vertex_index)) {
                        std::cout << "Delete vertex! [unimplemented]" << std::endl;
                        // map.RemoveVertex(*selected_vertex_index);
                        selected_vertex_index = {core::kInvalidIndex};
                    }
                    if (core::IsValid(selected_edge_index)) {
                        // Delete that edge.
                        std::cout << "Delete edge! [unimplemented]" << std::endl;
                        // map.RemoveDirectedEdge(*selected_edge_index);
                        selected_edge_index = {core::kInvalidIndex};
                    }
                }
            }
        }

        u32 color_white = 0xFFFFFFFF;
        u32 color_background = 0x414141FF;
        u32 color_light_background = 0x454545FF;
        u32 color_light_gray = 0x909090FF;
        u32 color_qe_constrained = 0xAAAACFFF;
        u32 color_qe_normal = 0xFF48CFFF;
        u32 color_active = 0xFFA0A0FF;

        // Clear screen
        SetColor(renderer, color_background);
        SDL_RenderClear(renderer);

        if (core::IsValid(qe_mouse_face)) {
            const core::DelaunayMesh& mesh = map.GetMesh();

            // Fill enclosing triangle
            const auto [qe_ab, qe_bc, qe_ca] = mesh.GetTriangleQuarterEdges(qe_mouse_face);
            const common::Vec2f& a = mesh.GetVertex(qe_ab);
            const common::Vec2f& b = mesh.GetVertex(qe_bc);
            const common::Vec2f& c = mesh.GetVertex(qe_ca);

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

        {
            // Draw major vertical lines
            f32 major_line_spacing = 1.0;
            f32 minor_line_spacing = major_line_spacing / 8.0;

            if (camera_zoom > 75.0) {
                RenderGrid(renderer, 0x404055FF, minor_line_spacing, camera_pos, camera_zoom);
            }

            RenderGrid(renderer, color_light_background, major_line_spacing, camera_pos,
                       camera_zoom);
        }

        {
            // Render the mesh
            const core::DelaunayMesh& mesh = map.GetMesh();

            core::QuarterEdgeIndex qe = mesh.GetFirstQuarterEdgeIndex();
            while (core::IsValid(qe)) {
                if (mesh.IsPrimal(qe)) {
                    // Get its opposite side.
                    core::QuarterEdgeIndex qe_sym = mesh.Sym(qe);

                    const core::VertexData& a = mesh.GetVertexData(qe);
                    const core::VertexData& b = mesh.GetVertexData(qe_sym);
                    if (a.i_self > b.i_self) {  // Avoid rendering edges twice

                        // Set the color
                        if (mesh.IsConstrained(qe)) {
                            SetColor(renderer, color_qe_constrained);
                        } else {
                            SetColor(renderer, color_qe_normal);
                        }

                        auto a_cam = GlobalToCamera(a.v, camera_pos, camera_zoom);
                        auto b_cam = GlobalToCamera(b.v, camera_pos, camera_zoom);
                        SDL_RenderDrawLine(renderer, (int)(a_cam.x), (int)(a_cam.y), (int)(b_cam.x),
                                           (int)(b_cam.y));
                    }
                }
                // Get the next one
                qe = mesh.GetNext(qe);
            }
        }

        // {  // Render all side_infos (these are directed)
        //     SDL_SetRenderDrawColor(renderer, 0x80, 0x80, 0x80, 0xFF);

        //     for (const auto& side_info : map.GetSideInfos()) {
        //         common::Vec2f a = map.GetVertices()[side_info.a_ind];
        //         common::Vec2f b = map.GetVertices()[side_info.b_ind];

        //         auto a_cam = GlobalToCamera(a, camera_pos, camera_zoom);
        //         auto b_cam = GlobalToCamera(b, camera_pos, camera_zoom);
        //         SDL_RenderDrawLine(renderer, (int)(a_cam.x), (int)(a_cam.y), (int)(b_cam.x),
        //                            (int)(b_cam.y));

        //         common::Vec2f c = (a + b) / 2.0;
        //         common::Vec2f d = c + Rotr(Normalize(b - a)) * 0.2;
        //         auto c_cam = GlobalToCamera(c, camera_pos, camera_zoom);
        //         auto d_cam = GlobalToCamera(d, camera_pos, camera_zoom);
        //         SDL_RenderDrawLine(renderer, (int)(c_cam.x), (int)(c_cam.y), (int)(d_cam.x),
        //                            (int)(d_cam.y));
        //     }
        // }

        {  // Render all vertices
            const core::DelaunayMesh& mesh = map.GetMesh();

            core::VertexIndex i_vertex = mesh.GetFirstVertexIndex();
            while (core::IsValid(i_vertex)) {
                auto v_cam = GlobalToCamera(mesh.GetVertex(i_vertex), camera_pos, camera_zoom);

                SDL_Rect rect;

                // Outline with darker color
                SetColor(renderer, color_background);
                rect.x = (int)(v_cam.x - 2);
                rect.y = (int)(v_cam.y - 2);
                rect.h = 5;
                rect.w = 5;
                SDL_RenderFillRect(renderer, &rect);

                // Fill with lighter color
                SetColor(renderer, color_light_gray);
                rect.x = (int)(v_cam.x - 1);
                rect.y = (int)(v_cam.y - 1);
                rect.h = 3;
                rect.w = 3;
                SDL_RenderFillRect(renderer, &rect);

                // Get the next one
                i_vertex = mesh.GetNext(i_vertex);
            }
        }

        if (IsValid(selected_edge_index)) {
            // Render our selected edge
            SetColor(renderer, color_white);

            const core::DelaunayMesh& mesh = map.GetMesh();
            common::Vec2f a = mesh.GetVertex(selected_edge_index);
            common::Vec2f b = mesh.GetVertex(mesh.Sym(selected_edge_index));

            auto a_cam = GlobalToCamera(a, camera_pos, camera_zoom);
            auto b_cam = GlobalToCamera(b, camera_pos, camera_zoom);
            SDL_RenderDrawLine(renderer, (int)(a_cam.x), (int)(a_cam.y), (int)(b_cam.x),
                               (int)(b_cam.y));

            // Draw a little extra tick mark for directed edges
            // common::Vec2f c = (a + b) / 2.0;
            // common::Vec2f d = c + Rotr(Normalize(b - a)) * 0.2;
            // auto c_cam = GlobalToCamera(c, camera_pos, camera_zoom);
            // auto d_cam = GlobalToCamera(d, camera_pos, camera_zoom);
            // SDL_RenderDrawLine(renderer, (int)(c_cam.x), (int)(c_cam.y), (int)(d_cam.x),
            //                    (int)(d_cam.y));
        }

        if (core::IsValid(selected_vertex_index)) {
            // Render our selected vertex
            const core::DelaunayMesh& mesh = map.GetMesh();
            const auto& v = mesh.GetVertex(selected_vertex_index);
            auto v_cam = GlobalToCamera(v, camera_pos, camera_zoom);

            SDL_Rect rect;

            // Outline with darker color
            SetColor(renderer, color_background);
            rect.x = (int)(v_cam.x - 3);
            rect.y = (int)(v_cam.y - 3);
            rect.h = 7;
            rect.w = 7;
            SDL_RenderFillRect(renderer, &rect);

            // Fill with lighter color
            if (mouse_is_pressed) {
                SetColor(renderer, color_active);
            } else {
                SetColor(renderer, color_white);
            }
            rect.x = (int)(v_cam.x - 2);
            rect.y = (int)(v_cam.y - 2);
            rect.h = 5;
            rect.w = 5;
            SDL_RenderFillRect(renderer, &rect);

            if (mouse_is_pressed) {
                // Render a line to the mouse position.
                auto b_cam = GlobalToCamera(mouse_pos, camera_pos, camera_zoom);
                SDL_RenderDrawLine(renderer, (int)(v_cam.x), (int)(v_cam.y), (int)(b_cam.x),
                                   (int)(b_cam.y));
            }
        }

        // SDL_RENDERER_PRESENTVSYNC means this is syncronized with the monitor
        // refresh rate. (30Hz)
        SDL_RenderPresent(renderer);
    }

    // Destroy our window
    SDL_DestroyWindow(window);
}