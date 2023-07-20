// TODO: Make it easy to add a new edge by:
//   1. Having selected edge panel contain content even if no side info exists
//   2. Giving the option to switch selected edge direction in the selected edge panel
//   3. Giving the option to create a side info in the panel

#include <SDL2/SDL.h>

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "assets_exporter.hpp"
#include "assets_importer.hpp"
#include "delaunay_mesh.hpp"
#include "game_map.hpp"
#include "geometry_utils.hpp"
#include "imgui.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_sdlrenderer2.h"
#include "math_utils.hpp"
#include "palette.hpp"
#include "typedefs.hpp"

#define EDITOR_SCREEN_SIZE_X 1280
#define EDITOR_SCREEN_SIZE_Y 720

#define PLAYER_SCREEN_SIZE_X 640
#define PLAYER_SCREEN_SIZE_Y 360

#define ASSERT(_e, ...)               \
    if (!(_e)) {                      \
        fprintf(stderr, __VA_ARGS__); \
        exit(1);                      \
    }

// ------------------------------------------------------------------------------------------------
common::Vec2f GlobalToCamera(const common::Vec2f& g, const common::Vec2f& camera_pos,
                             f32 camera_zoom) {
    common::Vec2f c = (g - camera_pos) * camera_zoom;
    return common::Vec2f(EDITOR_SCREEN_SIZE_X / 2 + c.x, EDITOR_SCREEN_SIZE_Y / 2 - c.y);
}

// ------------------------------------------------------------------------------------------------
common::Vec2f CameraToGlobal(const common::Vec2f& c, const common::Vec2f& camera_pos,
                             f32 camera_zoom) {
    common::Vec2f g_offset =
        common::Vec2f(c.x - EDITOR_SCREEN_SIZE_X / 2, EDITOR_SCREEN_SIZE_Y / 2 - c.y);
    return g_offset / camera_zoom + camera_pos;
}

// ------------------------------------------------------------------------------------------------
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

// ------------------------------------------------------------------------------------------------
void RenderGrid(SDL_Renderer* renderer, u32 rgba, f32 line_spacing, const common::Vec2f& camera_pos,
                f32 camera_zoom) {
    SetColor(renderer, rgba);

    f32 line_spacing_camera = line_spacing * camera_zoom;
    f32 x_global = camera_pos.x - fmod(camera_pos.x, line_spacing);
    f32 x_camera = EDITOR_SCREEN_SIZE_X / 2 + (x_global - camera_pos.x) * camera_zoom;
    x_camera = fmod(x_camera, line_spacing_camera);
    while (x_camera < EDITOR_SCREEN_SIZE_X) {
        SDL_RenderDrawLine(renderer, (int)(x_camera), (int)(0), (int)(x_camera),
                           (int)(EDITOR_SCREEN_SIZE_Y));
        x_camera += line_spacing_camera;
    }

    f32 y_global = camera_pos.y - fmod(camera_pos.y, line_spacing);
    f32 y_camera = EDITOR_SCREEN_SIZE_Y / 2 - (y_global - camera_pos.y) * camera_zoom;
    y_camera = fmod(y_camera, line_spacing_camera);
    while (y_camera < EDITOR_SCREEN_SIZE_Y) {
        SDL_RenderDrawLine(renderer, (int)(0), (int)(y_camera), (int)(EDITOR_SCREEN_SIZE_X),
                           (int)(y_camera));
        y_camera += line_spacing_camera;
    }
}

// ------------------------------------------------------------------------------------------------
void RenderTextureColumn(u32* pixels, int x, int screen_size_x, int y_lower, int y_upper, int y_lo,
                         int y_hi, f32 x_along_texture, u32 texture_x_offset, u32 texture_y_offset,
                         const core::OldStyleBitmap& bitmap) {
    u32 TEXTURE_SIZE = 64;  // TODO
    f32 TILE_WIDTH = 1.0f;
    f32 PIX_PER_DISTANCE = TEXTURE_SIZE / TILE_WIDTH;

    u32 texture_x = (int)(PIX_PER_DISTANCE * x_along_texture) % TEXTURE_SIZE;
    u32 baseline = bitmap.GetColumnMajorPixelIndex(texture_x + texture_x_offset, texture_y_offset);

    // TODO: Don't squish the texture. Instead
    // assume it has a fixed resolution.
    u32 denom = std::max(1, y_upper - y_lower);
    f32 y_loc = (f32)((y_upper - y_hi) * TEXTURE_SIZE) / denom;
    f32 y_step = (f32)(TEXTURE_SIZE) / denom;
    for (int y = y_hi; y >= y_lo; y--) {
        u32 texture_y = std::min((u32)(y_loc), TEXTURE_SIZE - 1);
        u32 color = bitmap.abgr[texture_y + baseline];
        pixels[(y * screen_size_x) + x] = color;
        y_loc += y_step;
    }
}

// ------------------------------------------------------------------------------------------------
void RenderWallsViaMesh(u32* pixels, f32* wall_raycast_radius, int screen_size_x, int screen_size_y,
                        const core::GameMap& game_map, const core::OldStyleBitmap& bitmap) {
    const core::DelaunayMesh mesh = game_map.GetMesh();

    // Camera data
    common::Vec2f camera_fov = {1.5, 0.84375};
    common::Vec2f camera_pos = {5.0, 5.0};
    common::Vec2f camera_dir = {1.0, 0.0};
    f32 camera_z = 0.4f;
    core::QuarterEdgeIndex qe_camera = mesh.GetEnclosingTriangle(camera_pos);
    if (!core::IsValid(qe_camera)) {
        return;
    }

    f32 half_screen_size = screen_size_y / 2.0f;
    f32 screen_size_y_over_fov_y = screen_size_y / camera_fov.y;

    u32 TEXTURE_SIZE = 64;
    u32 color_ceil = 0xFF222222;
    u32 color_floor = 0xFF444444;

    for (int x = 0; x < screen_size_x; x++) {
        // Camera to pixel column
        const f32 dw =
            camera_fov.x / 2 - (camera_fov.x * x) / screen_size_x;  // TODO: Precompute once.
        const common::Vec2f cp = {camera_dir.x - dw * camera_dir.y,
                                  camera_dir.y + dw * camera_dir.x};

        // Distance from the camera to the column
        const f32 cam_len = common::Norm(cp);

        // Ray direction through this column
        const common::Vec2f dir = cp / cam_len;

        // Start at the camera pos
        core::QuarterEdgeIndex qe_dual = qe_camera;
        common::Vec2f pos = camera_pos;

        // The edge vector of the face that we last crossed
        common::Vec2f v_face = {0.0, 0.0};

        // Step through triangles until we hit a solid triangle
        const core::SideInfo* side_info = nullptr;  // The side info we eventually hit.

        // f32 z_ceil = 999.0;                         // TODO: Use typemax
        // f32 z_floor = -999.0;                       // TODO: Use typemax
        int y_hi = screen_size_y;
        int y_lo = -1;

        int n_steps = 0;
        while (n_steps < 100) {
            n_steps += 1;

            // Grab the enclosing triangle.
            auto [qe_ab, qe_bc, qe_ca] = mesh.GetTriangleQuarterEdges(qe_dual);

            common::Vec2f a = mesh.GetVertex(qe_ab);
            common::Vec2f b = mesh.GetVertex(qe_bc);
            common::Vec2f c = mesh.GetVertex(qe_ca);

            // Project our ray out far enough that it would exit our mesh
            f32 projection_distance = 100.0;  // Ridiculously large
            common::Vec2f pos_next_delta = projection_distance * dir;
            common::Vec2f pos_next = pos + pos_next_delta;

            f32 min_interp = INFINITY;
            core::QuarterEdgeIndex qe_side = {core::kInvalidIndex};

            // See if we cross any of the 3 faces for the triangle we are in,
            // and cross the first segment.
            const f32 eps = 1e-4;
            if (common::GetRightHandedness(a, b, pos_next) < -eps) {
                // We would cross AB
                common::Vec2f v = b - a;
                common::Vec2f w = pos - a;
                float interp_ab = common::Cross(v, w) / common::Cross(pos_next_delta, v);
                if (interp_ab < min_interp) {
                    min_interp = interp_ab;
                    qe_side = qe_ab;
                    v_face = v;
                }
            }
            if (common::GetRightHandedness(b, c, pos_next) < -eps) {
                // We would cross BC
                common::Vec2f v = c - b;
                common::Vec2f w = pos - b;
                float interp_bc = common::Cross(v, w) / common::Cross(pos_next_delta, v);
                if (interp_bc < min_interp) {
                    min_interp = interp_bc;
                    qe_side = qe_bc;
                    v_face = v;
                }
            }
            if (common::GetRightHandedness(c, a, pos_next) < -eps) {
                // We would cross CA
                common::Vec2f v = a - c;
                common::Vec2f w = pos - c;
                float interp_ca = common::Cross(v, w) / common::Cross(pos_next_delta, v);
                if (interp_ca < min_interp) {
                    min_interp = interp_ca;
                    qe_side = qe_ca;
                    v_face = v;
                }
            }

            // Move to the face.
            if (core::IsValid(qe_side)) {
                // Should always be non-null.
                pos += min_interp * pos_next_delta;
                qe_dual = mesh.Rot(qe_side);  // The next face

                auto it = game_map.GetSideInfos().find(qe_side);
                if (it != game_map.GetSideInfos().end()) {
                    side_info = &(it->second);

                    const f32 ray_len = std::max(common::Norm(pos - camera_pos), 0.01f);
                    const f32 gamma = cam_len / ray_len * screen_size_y_over_fov_y;
                    wall_raycast_radius[x] = ray_len;

                    int y_ceil = (int)(half_screen_size + gamma * (side_info->z_ceil - camera_z));
                    int y_upper = (int)(half_screen_size + gamma * (side_info->z_upper - camera_z));
                    int y_lower = (int)(half_screen_size + gamma * (side_info->z_lower - camera_z));
                    int y_floor = (int)(half_screen_size + gamma * (side_info->z_floor - camera_z));

                    // Render the ceiling above the upper texture
                    while (y_hi > y_ceil) {
                        y_hi--;
                        pixels[(y_hi * screen_size_x) + x] = color_ceil;
                    }

                    // Render the upper texture
                    if (y_upper < y_hi) {
                        while (y_hi > y_upper) {
                            y_hi--;
                            pixels[(y_hi * screen_size_x) + x] = 0xFF0000FF;
                        }
                    }

                    // Render the floor below the lower texture
                    while (y_lo < y_floor) {
                        y_lo++;
                        pixels[(y_lo * screen_size_x) + x] = color_floor;
                    }

                    // Render the lower texture
                    if (y_lower > y_lo) {
                        while (y_lo < y_lower) {
                            y_lo++;
                            pixels[(y_lo * screen_size_x) + x] = 0x0000FFFF;
                        }
                    }

                    // Continue on with our projection if the side is passable.
                    if ((side_info->flags & core::kSideInfoFlag_PASSABLE) > 0) {
                        continue;
                    }

                    // The side info has a solid wall.
                    // Calculate where along the segment we intersected.
                    core::QuarterEdgeIndex qe_face_src = mesh.Tor(qe_dual);
                    f32 x_along_texture =
                        common::Norm(v_face) - common::Norm(pos - mesh.GetVertex(qe_face_src));
                    u32 texture_x_offset =
                        (side_info->flags & core::kSideInfoFlag_DARK) > 0 ? TEXTURE_SIZE : 0;
                    u32 texture_y_offset = side_info->texture_info_middle.texture_id * TEXTURE_SIZE;
                    texture_x_offset += side_info->texture_info_middle.x_offset;
                    texture_y_offset += side_info->texture_info_middle.y_offset;
                    RenderTextureColumn(pixels, x, screen_size_x, y_lower, y_upper, y_lo, y_hi,
                                        x_along_texture, texture_x_offset, texture_y_offset,
                                        bitmap);

                    break;
                }

                // Also break if it is the boundary
                if (mesh.IsBoundaryEdge(qe_side)) {
                    break;
                }
            } else {
                break;
            }
        }
    }
}

// ------------------------------------------------------------------------------------------------
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

// ------------------------------------------------------------------------------------------------
void ExportGameData(const core::GameMap& map) {
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Exporting game data" << std::endl;

    core::AssetsExporter exporter;
    std::cout << "Exporting core data" << std::endl;
    exporter.AddEntry(core::ExportColorPalette());

    std::cout << "Exporting map data" << std::endl;
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

// ------------------------------------------------------------------------------------------------
class SDLWindowData {
  public:
    SDLWindowData(const char* title, int size_x, int size_y) :
        screen_size_x(size_x), screen_size_y(size_y) {
        window = SDL_CreateWindow(title, SDL_WINDOWPOS_CENTERED_DISPLAY(1),
                                  SDL_WINDOWPOS_CENTERED_DISPLAY(1), size_x, size_y,
                                  SDL_WINDOW_ALLOW_HIGHDPI);
        ASSERT(window, "Error creating SDL window %s: %s\n", title, SDL_GetError());

        renderer =
            SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC | SDL_RENDERER_ACCELERATED);
        ASSERT(renderer, "Error creating SDL renderer for window %s: %s\n", title, SDL_GetError());

        texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STREAMING,
                                    size_x, size_y);
        ASSERT(texture, "Error creating SDL texture for window %s: %s\n", title, SDL_GetError());
    }

    // Destructor.
    ~SDLWindowData() { SDL_DestroyWindow(window); }

    int screen_size_x;
    int screen_size_y;
    SDL_Window* window;
    SDL_Renderer* renderer;
    SDL_Texture* texture;
};

int main() {
    SDL_version ver;
    SDL_GetVersion(&ver);
    fprintf(stdout, "Running with SDL2 version %d.%d.%d\n", ver.major, ver.minor, ver.patch);

    // Initialize SDL
    ASSERT(SDL_Init(SDL_INIT_VIDEO) == 0, "SDL initialization failed: %s\n", SDL_GetError());

    // From 2.0.18: Enable native IME.
#ifdef SDL_HINT_IME_SHOW_UI
    SDL_SetHint(SDL_HINT_IME_SHOW_UI, "1");
#endif

    // Create our windows
    SDLWindowData editor_window_data("TOOM EDITOR", EDITOR_SCREEN_SIZE_X, EDITOR_SCREEN_SIZE_Y);
    SDLWindowData player_window_data("PLAYER VIEW", PLAYER_SCREEN_SIZE_X, PLAYER_SCREEN_SIZE_Y);

    // Set up Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;   // Enable Gamepad Controls

    // Set up a Dear ImGui style
    ImGui::StyleColorsDark();
    ImGui::GetStyle().FrameRounding = 2.0f;

    // Set up the Dear ImGui Platform/Renderer backends
    ImGui_ImplSDL2_InitForSDLRenderer(editor_window_data.window, editor_window_data.renderer);
    ImGui_ImplSDLRenderer2_Init(editor_window_data.renderer);

    // Import our julia assets
    std::unique_ptr<core::AssetsImporter> julia_assets =
        core::AssetsImporter::LoadFromFile("../toom/assets/assets.bin");
    ASSERT(julia_assets, "Failed to load Julia Assets");

    // Extract our wall textures
    auto bitmap_data_opt = julia_assets->FindEntryData("textures");
    ASSERT(bitmap_data_opt, "Failed to find texture data");
    core::OldStyleBitmap bitmap = core::LoadBitmap(*bitmap_data_opt);

    // Create our map
    core::GameMap map;

    // Camera parameters
    common::Vec2f camera_pos = {2.0, 2.0};
    f32 camera_zoom = 50.0;

    // GUI state
    bool mouse_is_pressed = false;
    common::Vec2f mouse_pos = {0.0, 0.0};
    common::Vec2f mouse_click_pos = {0.0, 0.0};
    common::Vec2f camera_pos_at_mouse_click = {0.0, 0.0};
    core::QuarterEdgeIndex selected_vertex_index = {core::kInvalidIndex};
    core::QuarterEdgeIndex selected_edge_index = {core::kInvalidIndex};
    core::QuarterEdgeIndex qe_mouse_face = map.GetMesh().GetEnclosingTriangle(mouse_pos);

    // Player view data
    u32 player_view_pixels[player_window_data.screen_size_x *
                           player_window_data.screen_size_y];  // row-major
    f32 wall_raycast_radius[player_window_data.screen_size_x];

    bool continue_running = true;
    while (continue_running) {
        // Poll and handle events (inputs, window resize, etc.)
        // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui
        // wants to use your inputs.
        // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main
        //   application, or clear/overwrite your copy of the mouse data.
        // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main
        //   application, or clear/overwrite your copy of the keyboard data.
        // Generally you may always pass all inputs to dear imgui, and hide them from your
        // application based on those two flags.
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            ImGui_ImplSDL2_ProcessEvent(&event);
            if (event.type == SDL_QUIT) {
                continue_running = false;
                break;
            } else if (event.type == SDL_WINDOWEVENT &&
                       event.window.event == SDL_WINDOWEVENT_CLOSE &&
                       (event.window.windowID == SDL_GetWindowID(editor_window_data.window) ||
                        event.window.windowID == SDL_GetWindowID(player_window_data.window))) {
                continue_running = false;
                break;
            } else if (event.type == SDL_MOUSEWHEEL && !io.WantCaptureMouse) {
                if (event.wheel.preciseY > 0) {
                    camera_zoom *= 1.1;
                } else {
                    camera_zoom /= 1.1;
                }
            } else if (event.type == SDL_MOUSEBUTTONDOWN && !io.WantCaptureMouse) {
                if (event.button.button == SDL_BUTTON_LEFT) {
                    if (!mouse_is_pressed) {
                        // New press
                        mouse_is_pressed = true;
                        mouse_click_pos = CameraToGlobal(
                            common::Vec2f(event.button.x, event.button.y), camera_pos, camera_zoom);
                        camera_pos_at_mouse_click = camera_pos;

                        // Check for a selected vertex near the current mouse position
                        selected_edge_index = {core::kInvalidIndex};
                        selected_vertex_index =
                            map.FindVertexNearPosition(mouse_click_pos, qe_mouse_face);

                        // Check for a selected edge near the current mouse position if we did not
                        // select a vertex
                        if (!IsValid(selected_vertex_index)) {
                            selected_edge_index =
                                map.FindEdgeNearPosition(mouse_click_pos, qe_mouse_face);
                        }
                    }
                } else if (event.button.button == SDL_BUTTON_MIDDLE) {
                    // Check for a selected edge.
                    common::Vec2f click_pos = CameraToGlobal(
                        common::Vec2f(event.button.x, event.button.y), camera_pos, camera_zoom);
                    core::QuarterEdgeIndex edge =
                        map.FindEdgeNearPosition(click_pos, qe_mouse_face);
                    if (core::IsValid(edge)) {
                        if (!map.MaybeFlipEdge(edge)) {
                            std::cout << "Failed to flip edge" << std::endl;
                        }
                    }
                }
            } else if (event.type == SDL_MOUSEBUTTONUP && !io.WantCaptureMouse) {
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
            } else if (event.type == SDL_MOUSEMOTION && !io.WantCaptureMouse) {
                // Move the mouse
                mouse_pos = CameraToGlobal(
                    common::Vec2f(event.motion.x, event.motion.y),
                    mouse_is_pressed ? camera_pos_at_mouse_click : camera_pos, camera_zoom);
                // Update the mouse face
                qe_mouse_face = map.GetMesh().GetEnclosingTriangle(mouse_pos, qe_mouse_face);

                if (mouse_is_pressed) {
                    if (core::IsValid(selected_vertex_index)) {
                        // Move the given vertex
                        map.MoveVertexToward(selected_vertex_index, mouse_pos);
                    } else {
                        // Pan the camera
                        camera_pos = camera_pos_at_mouse_click + mouse_click_pos - mouse_pos;
                    }
                }

            } else if (event.type == SDL_KEYDOWN && !io.WantCaptureKeyboard) {
                if (event.key.keysym.sym == SDLK_e) {
                    ExportGameData(map);
                } else if (event.key.keysym.sym == SDLK_i) {
                    ImportGameData(&map);
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

        // ------------------------------------------------------------------------------------------------
        // Clear screen
        SetColor(editor_window_data.renderer, color_background);
        SDL_RenderClear(editor_window_data.renderer);

        if (core::IsValid(qe_mouse_face)) {
            auto renderer = editor_window_data.renderer;
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
            auto renderer = editor_window_data.renderer;
            f32 major_line_spacing = 1.0;
            f32 minor_line_spacing = major_line_spacing / 8.0;

            if (camera_zoom > 75.0) {
                RenderGrid(renderer, 0x404055FF, minor_line_spacing, camera_pos, camera_zoom);
            }

            RenderGrid(renderer, color_light_background, major_line_spacing, camera_pos,
                       camera_zoom);
        }

        {  // Render the mesh
            auto renderer = editor_window_data.renderer;
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

        {  // Render all side_infos (these are directed)
            auto renderer = editor_window_data.renderer;

            const core::DelaunayMesh& mesh = map.GetMesh();
            for (const auto& it : map.GetSideInfos()) {
                SetColor(renderer, color_light_gray);
                if ((it.second.flags & core::kSideInfoFlag_PASSABLE) > 0) {
                    SetColor(renderer, 0x9090C0FF);
                }

                core::QuarterEdgeIndex qe = it.second.qe;
                const common::Vec2f& a = mesh.GetVertex(qe);
                const common::Vec2f& b = mesh.GetVertex(mesh.Sym(qe));

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
            auto renderer = editor_window_data.renderer;
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
            auto renderer = editor_window_data.renderer;
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
            auto renderer = editor_window_data.renderer;
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

        // ------------------------------------------------------------------------------------------------
        // Start the Dear ImGui frame
        ImGui_ImplSDLRenderer2_NewFrame();
        ImGui_ImplSDL2_NewFrame();
        ImGui::NewFrame();

        // Vertex panel
        {
            static bool vertex_panel_drawn_last_frame = false;
            static common::Vec2f vertex_panel_target;
            if (core::IsValid(selected_vertex_index)) {
                // Vertex panel
                vertex_panel_drawn_last_frame = true;
                ImGui::Begin("Vertex");
                ImGui::Text("index:      %lu", selected_vertex_index.i);
                ImGui::Separator();

                const core::DelaunayMesh& mesh = map.GetMesh();
                const auto& v = mesh.GetVertex(selected_vertex_index);
                if (!vertex_panel_drawn_last_frame || !io.WantCaptureMouse) {
                    vertex_panel_target = v;
                }

                ImGui::InputFloat("x", &vertex_panel_target.x, 0.1f);
                ImGui::InputFloat("y", &vertex_panel_target.y, 0.1f);
                map.MoveVertexToward(selected_vertex_index, vertex_panel_target);

                ImGui::End();
            } else {
                vertex_panel_drawn_last_frame = false;
            }
        }

        // Side Info panel
        if (core::IsValid(selected_edge_index)) {
            core::SideInfo* side_info = map.GetEditableSideInfo(selected_edge_index);
            if (side_info != nullptr) {
                ImGui::Begin("SideInfo");
                ImGui::Text("index:      %lu", selected_edge_index.i);
                ImGui::Separator();
                if (ImGui::Button((
                        (side_info->flags & core::kSideInfoFlag_DARK) > 0 ? "Dark" : "Not Dark"))) {
                    side_info->flags ^= core::kSideInfoFlag_DARK;  // toggle
                }
                if (ImGui::Button(((side_info->flags & core::kSideInfoFlag_PASSABLE) > 0
                                       ? "Passable"
                                       : "Not Passable"))) {
                    side_info->flags ^= core::kSideInfoFlag_PASSABLE;  // toggle
                }

                ImGui::Text("flags:      %X", side_info->flags);

                int flags = 0;
                u16 step_u16 = 1;
                i16 step_i16 = 1;
                f32 step_f32 = 0.1f;

                ImGui::Separator();
                if (ImGui::InputScalar("upper texture_id", ImGuiDataType_U16,
                                       (void*)(&side_info->texture_info_upper.texture_id),
                                       (void*)(&step_u16), (void*)(NULL), "%d", flags)) {
                    // Ensure it is in bounds. TODO: Clamp by number of textures we have.
                    if (side_info->texture_info_upper.texture_id > 17) {
                        side_info->texture_info_upper.texture_id = 17;
                    }
                }

                ImGui::InputScalar("upper x_offset", ImGuiDataType_S16,
                                   (void*)(&side_info->texture_info_upper.x_offset),
                                   (void*)(&step_i16), (void*)(NULL), "%d", flags);
                ImGui::InputScalar("upper y_offset", ImGuiDataType_S16,
                                   (void*)(&side_info->texture_info_upper.y_offset),
                                   (void*)(&step_i16), (void*)(NULL), "%d", flags);

                ImGui::Separator();
                if (ImGui::InputScalar("middle texture_id", ImGuiDataType_U16,
                                       (void*)(&side_info->texture_info_middle.texture_id),
                                       (void*)(&step_u16), (void*)(NULL), "%d", flags)) {
                    // Ensure it is in bounds. TODO: Clamp by number of textures we have.
                    if (side_info->texture_info_middle.texture_id > 17) {
                        side_info->texture_info_middle.texture_id = 17;
                    }
                }
                ImGui::InputScalar("middle x_offset", ImGuiDataType_S16,
                                   (void*)(&side_info->texture_info_middle.x_offset),
                                   (void*)(&step_i16), (void*)(NULL), "%d", flags);
                ImGui::InputScalar("middle y_offset", ImGuiDataType_S16,
                                   (void*)(&side_info->texture_info_middle.y_offset),
                                   (void*)(&step_i16), (void*)(NULL), "%d", flags);

                ImGui::Separator();
                if (ImGui::InputScalar("lower texture_id", ImGuiDataType_U16,
                                       (void*)(&side_info->texture_info_lower.texture_id),
                                       (void*)(&step_u16), (void*)(NULL), "%d", flags)) {
                    // Ensure it is in bounds. TODO: Clamp by number of textures we have.
                    if (side_info->texture_info_lower.texture_id > 17) {
                        side_info->texture_info_lower.texture_id = 17;
                    }
                }
                ImGui::InputScalar("lower x_offset", ImGuiDataType_S16,
                                   (void*)(&side_info->texture_info_lower.x_offset),
                                   (void*)(&step_i16), (void*)(NULL), "%d", flags);
                ImGui::InputScalar("lower y_offset", ImGuiDataType_S16,
                                   (void*)(&side_info->texture_info_lower.y_offset),
                                   (void*)(&step_i16), (void*)(NULL), "%d", flags);

                ImGui::Separator();
                if (ImGui::InputScalar("z_ceil", ImGuiDataType_Float, (void*)(&side_info->z_ceil),
                                       (void*)(&step_f32), (void*)(NULL), "%.3f", flags)) {
                    side_info->z_ceil = std::max(side_info->z_ceil, side_info->z_upper);
                }
                if (ImGui::InputScalar("z_upper", ImGuiDataType_Float, (void*)(&side_info->z_upper),
                                       (void*)(&step_f32), (void*)(NULL), "%.3f", flags)) {
                    side_info->z_upper =
                        std::clamp(side_info->z_upper, side_info->z_lower, side_info->z_ceil);
                }
                if (ImGui::InputScalar("z_lower", ImGuiDataType_Float, (void*)(&side_info->z_lower),
                                       (void*)(&step_f32), (void*)(NULL), "%.3f", flags)) {
                    side_info->z_lower =
                        std::clamp(side_info->z_lower, side_info->z_floor, side_info->z_upper);
                }
                if (ImGui::InputScalar("z_floor", ImGuiDataType_Float, &side_info->z_floor,
                                       (void*)(&step_f32), (void*)(NULL), "%.3f", flags)) {
                    side_info->z_floor = std::min(side_info->z_floor, side_info->z_lower);
                }

                ImGui::End();
            }
        }

        // ImGUI Rendering
        ImGui::Render();
        SDL_RenderSetScale(editor_window_data.renderer, io.DisplayFramebufferScale.x,
                           io.DisplayFramebufferScale.y);
        ImGui_ImplSDLRenderer2_RenderDrawData(ImGui::GetDrawData());

        // ------------------------------------------------------------------------------------------------
        {
            // Render the player view.
            RenderWallsViaMesh(player_view_pixels, wall_raycast_radius,
                               player_window_data.screen_size_x, player_window_data.screen_size_y,
                               map, bitmap);

            SDL_UpdateTexture(player_window_data.texture, NULL, player_view_pixels,
                              player_window_data.screen_size_x * 4);
            SDL_RenderCopyEx(player_window_data.renderer, player_window_data.texture, NULL, NULL,
                             0.0, NULL, SDL_FLIP_VERTICAL);
        }

        // SDL_RENDERER_PRESENTVSYNC means this is syncronized with the monitor
        // refresh rate. (30Hz)
        SDL_RenderPresent(editor_window_data.renderer);
        SDL_RenderPresent(player_window_data.renderer);
    }
}