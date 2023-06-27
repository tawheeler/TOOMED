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

void ImportGameData(core::GameMap* map) {
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Importing game data" << std::endl;

    core::AssetsExporter exporter;
    exporter.LoadAssetsFile("../toom/assets/toomed.bin");
    std::cout << "Num entries:" << exporter.NumEntries() << std::endl;

    bool succeeded = map->Import(exporter);
    if (!succeeded) {
        std::cout << "Failed to load game map! Clearing possibly corrupted map." << std::endl;
        map->Clear();
    }

    std::cout << "DONE" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
}

void ExportGameData(const core::GameMap& map) {
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Exporting game data" << std::endl;

    core::AssetsExporter exporter;
    bool succeeded = map.Export(&exporter);
    if (!succeeded) {
        std::cout << "Failure while exporting game map! Not written to file." << std::endl;
        return;
    }

    std::cout << "Num entries:" << exporter.NumEntries() << std::endl;

    exporter.WriteToFile("../toom/assets/toomed.bin");

    std::cout << "DONE" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
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
    core::GameMap map;

    // TODO
    // map.vertices.emplace_back(1.0, 1.0);
    // map.vertices.emplace_back(7.0, 1.0);
    // map.vertices.emplace_back(7.0, 3.0);
    // map.vertices.emplace_back(6.0, 3.0);
    // map.vertices.emplace_back(6.0, 4.0);
    // map.vertices.emplace_back(9.0, 4.0);
    // map.vertices.emplace_back(9.0, 3.0);
    // map.vertices.emplace_back(9.0, 1.0);
    // map.vertices.emplace_back(12.0, 1.0);
    // map.vertices.emplace_back(12.0, 7.0);
    // map.vertices.emplace_back(9.0, 6.0);
    // map.vertices.emplace_back(9.0, 5.0);
    // map.vertices.emplace_back(6.0, 5.0);
    // map.vertices.emplace_back(6.0, 6.0);
    // map.vertices.emplace_back(7.0, 6.0);
    // map.vertices.emplace_back(7.0, 7.0);
    // map.vertices.emplace_back(1.0, 7.0);
    // map.vertices.emplace_back(0.0, 0.0);
    // map.vertices.emplace_back(13.0, 0.0);
    // map.vertices.emplace_back(13.0, 8.0);
    // map.vertices.emplace_back(0.0, 8.0);

    // map.vertices.emplace_back(2.0, 2.0);
    // map.vertices.emplace_back(3.0, 2.0);
    // map.vertices.emplace_back(3.0, 3.0);
    // map.vertices.emplace_back(2.0, 3.0);

    // map.vertices.emplace_back(3.0, 5.0);
    // map.vertices.emplace_back(4.0, 5.0);
    // map.vertices.emplace_back(4.0, 6.0);
    // map.vertices.emplace_back(3.0, 6.0);

    // // Add all vertices
    // for (const auto& v : map.vertices) {
    //     int new_vertex = map.mesh.AddDelaunayVertex(v);
    //     ASSERT(new_vertex != core::kInvalidIndex, "Failed to add vertex into mesh\n");
    // }

    // // Constrain all edges
    // for (const auto& side_info : map.side_infos) {
    //     ASSERT(map.mesh.ConstrainEdge(side_info.a_ind, side_info.b_ind),
    //            "Failed to constrain edge %d -> %d\n", (int)(side_info.a_ind),
    //            (int)(side_info.b_ind));
    // }

    // Camera parameters
    common::Vec2f camera_pos = {2.0, 2.0};
    f32 camera_zoom = 50.0;

    // Mouse press params
    bool mouse_is_pressed = false;
    common::Vec2f mouse_pos = {0.0, 0.0};
    common::Vec2f mouse_click_pos = {0.0, 0.0};
    common::Vec2f camera_pos_at_mouse_click = {0.0, 0.0};
    std::optional<usize> selected_vertex_index = std::nullopt;
    std::optional<usize> selected_edge_index = std::nullopt;
    core::QuarterEdge* qe_mouse_face = nullptr;

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
                    selected_edge_index = std::nullopt;
                    if (!selected_vertex_index) {
                        selected_edge_index =
                            map.FindEdgeNearPosition(mouse_click_pos, qe_mouse_face);
                    }
                }
            } else if (event.type == SDL_MOUSEBUTTONUP) {
                if (mouse_is_pressed) {
                    // New release
                    mouse_is_pressed = false;

                    if (selected_vertex_index) {
                        // Check for a selected vertex near the released mouse position
                        std::optional<usize> released_vertex_index =
                            map.FindVertexNearPosition(mouse_pos, qe_mouse_face);
                        if (released_vertex_index.has_value() &&
                            released_vertex_index != selected_vertex_index) {
                            // Join those edges
                            int src = *selected_vertex_index;
                            int dst = *released_vertex_index;
                            if (!map.HasEdge(src, dst)) {
                                map.AddDirectedEdge(src, dst);
                            }

                            // Recompute the mesh if necessary
                            // TODO
                            // if (!map.mesh.HasEdge(*selected_vertex_index,
                            // *released_vertex_index)) {
                            map.InvalidateMesh();
                            // }
                        }
                    }
                }
            } else if (event.type == SDL_MOUSEMOTION) {
                // Move the mouse
                mouse_pos = CameraToGlobal(common::Vec2f(event.motion.x, event.motion.y),
                                           camera_pos, camera_zoom);
                // Update the mouse face
                if (map.HasMesh()) {
                    if (qe_mouse_face != nullptr) {
                        qe_mouse_face =
                            map.GetMesh()->GetEnclosingTriangle(mouse_pos, qe_mouse_face);
                    } else {
                        qe_mouse_face = map.GetMesh()->GetEnclosingTriangle(mouse_pos);
                    }
                } else {
                    qe_mouse_face = nullptr;
                }

                // Pan the camera
                if (mouse_is_pressed && !selected_vertex_index.has_value()) {
                    camera_pos = camera_pos_at_mouse_click + mouse_click_pos - mouse_pos;
                }
            } else if (event.type == SDL_KEYDOWN) {
                if (event.key.keysym.sym == SDLK_e) {
                    ExportGameData(map);
                } else if (event.key.keysym.sym == SDLK_i) {
                    ImportGameData(&map);
                } else if (event.key.keysym.sym == SDLK_DELETE) {
                    // Delete key pressed!
                    if (selected_vertex_index) {
                        std::cout << "Delete vertex! [NOT IMPLEMENTED]" << std::endl;
                        // TODO map.RemoveVertex(*selected_vertex_index);
                        selected_vertex_index = std::nullopt;
                    }
                    if (selected_edge_index) {
                        // Delete that edge.
                        std::cout << "Delete edge!" << std::endl;
                        map.RemoveDirectedEdge(*selected_edge_index);
                        selected_edge_index = std::nullopt;
                    }
                }
            }
        }

        // Clear screen
        SDL_SetRenderDrawColor(renderer, 0x41, 0x41, 0x41, 0xFF);
        SDL_RenderClear(renderer);

        if (map.HasMesh() && qe_mouse_face != nullptr) {
            core::DelaunayMesh* mesh = map.GetMesh();

            // Fill enclosing triangle
            const common::Vec2f& a = mesh->GetTriangleVertex1(qe_mouse_face);
            const common::Vec2f& b = mesh->GetTriangleVertex2(qe_mouse_face);
            const common::Vec2f& c = mesh->GetTriangleVertex3(qe_mouse_face);

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

        if (map.HasMesh()) {
            // Render the mesh
            core::DelaunayMesh* mesh = map.GetMesh();

            SDL_SetRenderDrawColor(renderer, 0xFF, 0x48, 0xCF, 0xFF);

            for (size_t qe_index = 0; qe_index < mesh->NumQuarterEdges(); qe_index++) {
                core::QuarterEdge* qe = mesh->GetQuarterEdge(qe_index);

                if (core::IsPrimalEdge(*qe) && !mesh->IsBoundaryVertex(qe->vertex)) {
                    // Get its opposite side.
                    core::QuarterEdge* qe_sym = mesh->Sym(qe);

                    const core::VertexData* a = qe->vertex;
                    const core::VertexData* b = qe_sym->vertex;
                    if (a > b && !mesh->IsBoundaryVertex(b)) {  // Avoid rendering edges twice
                        auto a_cam = GlobalToCamera(a->vertex, camera_pos, camera_zoom);
                        auto b_cam = GlobalToCamera(b->vertex, camera_pos, camera_zoom);
                        SDL_RenderDrawLine(renderer, (int)(a_cam.x), (int)(a_cam.y), (int)(b_cam.x),
                                           (int)(b_cam.y));
                    }
                }
            }
        }

        {  // Render all side_infos (these are directed)
            SDL_SetRenderDrawColor(renderer, 0x80, 0x80, 0x80, 0xFF);

            for (const auto& side_info : map.GetSideInfos()) {
                common::Vec2f a = map.GetVertices()[side_info.a_ind];
                common::Vec2f b = map.GetVertices()[side_info.b_ind];

                auto a_cam = GlobalToCamera(a, camera_pos, camera_zoom);
                auto b_cam = GlobalToCamera(b, camera_pos, camera_zoom);
                SDL_RenderDrawLine(renderer, (int)(a_cam.x), (int)(a_cam.y), (int)(b_cam.x),
                                   (int)(b_cam.y));

                common::Vec2f c = (a + b) / 2.0;
                common::Vec2f d = c + Rotr(Normalize(b - a)) * 0.2;
                auto c_cam = GlobalToCamera(c, camera_pos, camera_zoom);
                auto d_cam = GlobalToCamera(d, camera_pos, camera_zoom);
                SDL_RenderDrawLine(renderer, (int)(c_cam.x), (int)(c_cam.y), (int)(d_cam.x),
                                   (int)(d_cam.y));
            }
        }

        {  // Render all vertices

            SDL_SetRenderDrawColor(renderer, 0x90, 0x90, 0x90, 0xFF);

            for (const auto& v : map.GetVertices()) {
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

        if (selected_edge_index) {
            // Render our selected edge
            SDL_SetRenderDrawColor(renderer, 0xFF, 0xFF, 0xFF, 0xFF);

            const auto& side_info = map.GetSideInfos()[selected_edge_index.value()];
            common::Vec2f a = map.GetVertices()[side_info.a_ind];
            common::Vec2f b = map.GetVertices()[side_info.b_ind];

            auto a_cam = GlobalToCamera(a, camera_pos, camera_zoom);
            auto b_cam = GlobalToCamera(b, camera_pos, camera_zoom);
            SDL_RenderDrawLine(renderer, (int)(a_cam.x), (int)(a_cam.y), (int)(b_cam.x),
                               (int)(b_cam.y));

            common::Vec2f c = (a + b) / 2.0;
            common::Vec2f d = c + Rotr(Normalize(b - a)) * 0.2;
            auto c_cam = GlobalToCamera(c, camera_pos, camera_zoom);
            auto d_cam = GlobalToCamera(d, camera_pos, camera_zoom);
            SDL_RenderDrawLine(renderer, (int)(c_cam.x), (int)(c_cam.y), (int)(d_cam.x),
                               (int)(d_cam.y));
        }

        if (selected_vertex_index) {
            // Render our selected vertex
            const auto& v = map.GetVertices()[*selected_vertex_index];
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
            if (mouse_is_pressed) {
                SDL_SetRenderDrawColor(renderer, 0xFF, 0xA0, 0xA0, 0xFF);
            } else {
                SDL_SetRenderDrawColor(renderer, 0xFF, 0xFF, 0xFF, 0xFF);
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