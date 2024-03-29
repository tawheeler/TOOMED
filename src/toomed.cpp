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
#include "input.hpp"
#include "math_utils.hpp"
#include "palette.hpp"
#include "platform_metrics.cpp"
#include "profiler.cpp"
#include "render_assets.hpp"
#include "texture.hpp"
#include "typedefs.hpp"
#include "wad_importer.hpp"

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
bool LoadDoomLevel(core::GameMap* map, const std::string& name,
                   const std::unique_ptr<core::WadImporter>& importer,
                   const core::RenderAssets& render_assets) {
    auto idx_opt = importer->FindEntryDataIndex(name);
    if (!idx_opt) {
        return false;
    }

    // The map data lumps are after the map name lump
    int lump_index = idx_opt.value();

    lump_index += 2;  // skip THINGS
    auto linedefs_opt = importer->GetEntryData(lump_index);
    lump_index += 1;
    auto sidedefs_opt = importer->GetEntryData(lump_index);
    lump_index += 1;
    auto vertexes_opt = importer->GetEntryData(lump_index);
    lump_index += 4;  // skip segs, subsectors, and nodes
    auto sectors_opt = importer->GetEntryData(lump_index);

    if (!linedefs_opt || !sidedefs_opt || !vertexes_opt || !sectors_opt) {
        return false;
    }

    const auto& [linedefs_data, linedefs_data_size] = *linedefs_opt;
    const auto& [sidedefs_data, sidedefs_data_size] = *sidedefs_opt;
    const auto& [vertexes_data, vertexes_data_size] = *vertexes_opt;
    const auto& [sectors_data, sectors_data_size] = *sectors_opt;
    map->LoadFromDoomData(linedefs_data, linedefs_data_size, sidedefs_data, sidedefs_data_size,
                          vertexes_data, vertexes_data_size, sectors_data, sectors_data_size,
                          render_assets);

    return true;
}

// ------------------------------------------------------------------------------------------------
struct CameraState {
    common::Vec2f pos;
    common::Vec2f dir;
    common::Vec2f vel;
    common::Vec2f fov;  // field of vew
    f32 height;         // height above the floor
    f32 z;              // absolute height in global
    f32 omega;
    core::QuarterEdgeIndex qe;
};

void MoveCamera(CameraState* camera_state, const core::KeyBoardState& keyboard_state,
                const core::GameMap& game_map, f32 dt) {
    TimeFunction;

    // ------------------------------------------------------------
    // Process Input
    common::Vec2f input_dir = {0.0,
                               0.0};  // In the body frame, which is right-handed, so y points left.
    if (IsPressed(keyboard_state.w)) {
        input_dir.x += 1.0;
    }
    if (IsPressed(keyboard_state.s)) {
        input_dir.x -= 1.0;
    }
    if (IsPressed(keyboard_state.d)) {
        input_dir.y -= 1.0;
    }
    if (IsPressed(keyboard_state.a)) {
        input_dir.y += 1.0;
    }

    int input_rot_dir = 0;  // Right-hand rotation in plane (CCW)
    if (IsPressed(keyboard_state.q)) {
        input_rot_dir += 1;
    }
    if (IsPressed(keyboard_state.e)) {
        input_rot_dir -= 1;
    }

    // ------------------------------------------------------------
    // Update the velocity and angular speed
    const f32 kPlayerInputAccel = 6.5;
    const f32 kPlayerInputAngularAccel = 8.5;
    const f32 kPlayerMaxSpeed = 5.0;
    const f32 kPlayerMaxOmega = 5.0;
    const f32 kAirFriction = 4.0;
    const f32 kAirFrictionRot = 4.0;
    const f32 kSlidingFriction = 0.95;

    // Note: Speed is in the global frame
    camera_state->vel.x += (camera_state->dir.x * input_dir.x - camera_state->dir.y * input_dir.y) *
                           kPlayerInputAccel * dt;
    camera_state->vel.y += (camera_state->dir.y * input_dir.x + camera_state->dir.x * input_dir.y) *
                           kPlayerInputAccel * dt;
    camera_state->omega += input_rot_dir * kPlayerInputAngularAccel * dt;

    // Clamp the velocity to a maximum magnitude
    f32 speed = common::Norm(camera_state->vel);
    if (speed > kPlayerMaxSpeed) {
        camera_state->vel.x *= kPlayerMaxSpeed / speed;
        camera_state->vel.y *= kPlayerMaxSpeed / speed;
    }
    if (camera_state->omega > kPlayerMaxOmega) {
        camera_state->omega *= kPlayerMaxOmega / camera_state->omega;
    } else if (camera_state->omega < -kPlayerMaxOmega) {
        camera_state->omega *= -kPlayerMaxOmega / camera_state->omega;
    }

    // ------------------------------------------------------------
    // Move the player.
    const core::DelaunayMesh& mesh = game_map.GetMesh();

    f32 dt_remaining = dt;
    while (dt_remaining > 0.0f) {
        common::Vec2f pos_delta = dt_remaining * camera_state->vel;
        common::Vec2f pos_next = camera_state->pos + pos_delta;

        // Exit if the step is too small, to avoid math problems.
        if (SqNorm(pos_delta) < 1e-6) {
            break;
        }

        // Check to see if we cross over any mesh edges.
        auto [qe_ab, qe_bc, qe_ca] = mesh.GetTriangleQuarterEdges(camera_state->qe);

        common::Vec2f a = mesh.GetVertex(qe_ab);
        common::Vec2f b = mesh.GetVertex(qe_bc);
        common::Vec2f c = mesh.GetVertex(qe_ca);

        // The fraction of the distance we will move
        f32 interp = 1.0f;
        core::QuarterEdgeIndex qe_side = {core::kInvalidIndex};
        common::Vec2f v_face = {0.0, 0.0};

        const f32 eps = 1e-4;
        if (common::GetRightHandedness(a, b, pos_next) < -eps) {
            // We would cross AB
            common::Vec2f v = b - a;
            common::Vec2f w = camera_state->pos - a;
            float interp_ab = common::Cross(v, w) / common::Cross(pos_delta, v);
            if (interp_ab < interp) {
                interp = interp_ab;
                qe_side = qe_ab;
                v_face = v;
            }
        }
        if (common::GetRightHandedness(b, c, pos_next) < -eps) {
            // We would cross BC
            common::Vec2f v = c - b;
            common::Vec2f w = camera_state->pos - b;
            float interp_bc = common::Cross(v, w) / common::Cross(pos_delta, v);
            if (interp_bc < interp) {
                interp = interp_bc;
                qe_side = qe_bc;
                v_face = v;
            }
        }
        if (common::GetRightHandedness(c, a, pos_next) < -eps) {
            // We would cross CA
            common::Vec2f v = a - c;
            common::Vec2f w = camera_state->pos - c;
            float interp_ca = common::Cross(v, w) / common::Cross(pos_delta, v);
            if (interp_ca < interp) {
                interp = interp_ca;
                qe_side = qe_ca;
                v_face = v;
            }
        }

        // Move the player
        // If we would have crossed any edge, this merely moves us up to the boundary instead
        camera_state->pos += interp * pos_delta;
        dt_remaining *= (1.0f - interp);
        dt_remaining -= 1e-5;  // eps factor for safety

        if (core::IsValid(qe_side)) {
            // We would have crossed into another triangle.

            // Determine whether to continue on or stop at the edge.
            bool stop_at_edge = false;
            bool is_portal = false;
            const core::SideInfo* side_info = game_map.GetSideInfo(qe_side);
            if (side_info != nullptr) {
                // TODO: Prevent proceeding if the z-difference is sufficiently large
                const bool is_passable = ((side_info->flags & core::kSideInfoFlag_PASSABLE) > 0);
                stop_at_edge = !is_passable;
                is_portal = is_passable && core::IsValid(side_info->qe_portal);
            }

            // prevent stopping for now
            stop_at_edge = false;

            if (is_portal) {
                // Calculate where along the segment we intersected.
                f32 v_face_len = common::Norm(v_face);
                f32 x_along_texture = v_face_len - common::Norm(pos_next - mesh.GetVertex(qe_side));

                // Update our quarter edge
                camera_state->qe = mesh.Tor(side_info->qe_portal);

                // The vector along the new face is
                common::Vec2f c = mesh.GetVertex(side_info->qe_portal);
                common::Vec2f d = mesh.GetVertex(mesh.Sym(side_info->qe_portal));
                common::Vec2f cd = d - c;
                cd = Normalize(cd);

                // Our new position is:
                camera_state->pos = c + x_along_texture * cd;

                f32 angle_ab = std::atan2(v_face.y, v_face.x);
                f32 angle_cd = std::atan2(cd.y, cd.x);
                f32 angle_rot = common::CounterClockwiseAngleDistance(angle_ab, angle_cd) + M_PI;
                f32 cos_rot = std::cos(angle_rot);
                f32 sin_rot = std::sin(angle_rot);

                // Now to get the new direction, we need to rotate out of the new face
                camera_state->dir = {cos_rot * camera_state->dir.x - sin_rot * camera_state->dir.y,
                                     sin_rot * camera_state->dir.x + cos_rot * camera_state->dir.y};

                // We also need to rotate our speed
                camera_state->vel = {cos_rot * camera_state->vel.x - sin_rot * camera_state->vel.y,
                                     sin_rot * camera_state->vel.x + cos_rot * camera_state->vel.y};

            } else if (stop_at_edge) {
                // The new triangle is solid, so do not change triangles.
                // Lose all velocity into the boundary surface.
                // This results in the vector along the face.
                camera_state->vel = VectorProjection(camera_state->vel, v_face);
                camera_state->vel.x *= kSlidingFriction;
                camera_state->vel.y *= kSlidingFriction;
            } else {
                // Accept the new triangle
                camera_state->qe = mesh.Rot(qe_side);
            }
        }
    }

    // camera_state->qe = mesh.GetEnclosingTriangle(camera_state->pos, camera_state->qe);

    // Update height
    const core::SideInfo* side_info = game_map.GetSideInfo(mesh.Rot(camera_state->qe));
    if (side_info != nullptr) {
        const core::Sector* sector = game_map.GetSector(side_info->sector_id);
        ASSERT(sector, "sector is null!");
        camera_state->z = camera_state->height + sector->z_floor;
    } else {
        // TODO: Remove.
        // camera_state->z = camera_state->height;
    }

    // Update the camera's rotational heading
    f32 theta = atan2(camera_state->dir.y, camera_state->dir.x);
    theta += camera_state->omega * dt;
    camera_state->dir = common::Vec2f(cos(theta), sin(theta));

    // Apply air friction
    f32 air_friction_decay = exp(-kAirFriction * dt);
    camera_state->vel.x *= air_friction_decay;
    camera_state->vel.y *= air_friction_decay;
    camera_state->omega *= exp(-kAirFrictionRot * dt);
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
struct RenderData {
    u32 n_pixels_per_world_unit;  // Number of pixels across a texture on a 1 unit span should be.
    u32 screen_size_x;            // The number of pixels across the canvas is
    u32 screen_size_y;            // The number of pixels high the canvas is
    u32 half_screen_size_y;
    f32 darkness_per_world_dist;  // for the 32 colormaps, by increasing darkness
    u32 max_render_steps;         // Maximum number of triangles we will traverse
    u32 color_ceil;
    u32 color_floor;
};

struct ActiveSpanData {
    common::Vec2f hit;   // Location of ray_dir_lo's intersection with the ground.
    common::Vec2f step;  // Amount to shift 'hit' per pix column moved toward hi's intersection.
    u16 sector_id;       // The sector id for this span.
    u16 flat_id;         // The flat id for this span (floor or ceiling)
    u8 colormap_index;
    bool is_active;  // Whether there is an active span.
    int x_start;     // The camera x pixel index where this span starts.
    int x_end;       // The most recent camera x pixel where this span currently ends.
};

void RenderPatchColumn(u32* pixels, int x_screen, int y_lower, int y_upper, int y_lo, int y_hi,
                       f32 x_along_texture, u32 x_base_offset, u32 y_base_offset, f32 column_height,
                       const doom::Patch& patch, const core::Palette& palette,
                       const core::Colormap& colormap, const RenderData& render_data) {
    TimeFunction;
    // TODO: for now, ignore y_base_offset.

    // The index into the patch of the column that we want to draw.
    u32 x_patch = ((int)(render_data.n_pixels_per_world_unit * x_along_texture) + x_base_offset) %
                  patch.size_x;

    // y_lower = screen y coordinate of bottom of column (can exceed screen bounds)
    // y_upper = screen y coordinate of top of column (can exceed screen bounds)
    // y_lo    = screen y coordinate where we end drawing (does not exceed screen bounds)
    // y_hi    = screen y coordinate where we start drawing (does not exceed screen bounds)
    // column_height = real-world height of the painting surface

    // How many patch pixels high the column is.
    f32 y_patch_height_pix = render_data.n_pixels_per_world_unit * column_height;

    // y_patch = m * y_screen + b for converting screen to patch coordinate
    //                         0 = m * y_upper + b   -> b = -m * y_upper
    //    y_patch_height_pix - 1 = m * y_lower + b
    //                           = m * y_lower - m * y_upper
    //                           = m * (y_lower - y_upper)
    // m = (y_patch_height_pix - 1) / (y_lower - y_upper)

    f32 m = (f32)(y_patch_height_pix - 1.0f) / (y_lower - y_upper);
    if (m >= 0.0) {
        return;  // This can happen if the y patch height is zero.
    }

    // The number of (continuous) patch pixels y changes per screen pixel
    f32 y_patch_step_per_screen_pixel = m;  // If y_screen goes up by 1, y_patch goes up this much
    f32 y_screen_step_per_patch_pixel = 1.0f / m;

    f32 y_patch = 0.0f;
    f32 y_screen = y_upper;
    while (y_screen > y_lo) {
        u32 column_offset = patch.column_offsets[x_patch];
        if (patch.post_data[column_offset] == 0xFF) {
            break;
        }
        while (patch.post_data[column_offset] != 0xFF) {
            u8 y_patch_delta = patch.post_data[column_offset];
            column_offset++;
            int post_length = patch.post_data[column_offset];
            column_offset++;

            // skip transparent pixels
            y_patch += y_patch_delta;
            y_screen += y_patch_delta * y_screen_step_per_patch_pixel;

            // process the post. We have `post_length` pixels to draw
            while (post_length > 0 && y_screen > y_lo) {
                // Keep decreasing y_screen (and increasing y_patch) as long as we are within the
                // post data.

                // Render pixels
                u8 palette_index = patch.post_data[column_offset];
                palette_index = colormap.map[palette_index];
                u32 abgr = *(u32*)(palette.rgbs + 3 * palette_index) | 0xFF000000;

                // Render this color for all screen pixels that map to y_patch
                u16 y_patch_discrete = (u16)y_patch;
                int y_screen_discrete = (int)y_screen;
                while ((u16)y_patch == y_patch_discrete) {
                    if (y_screen_discrete < y_hi && y_screen_discrete > y_lo) {
                        pixels[(y_screen_discrete * render_data.screen_size_x) + x_screen] = abgr;
                    }
                    y_screen_discrete -= 1;
                    y_screen -= 1.0f;
                    y_patch -= y_patch_step_per_screen_pixel;
                }

                u16 patch_delta = (u16)y_patch - y_patch_discrete;
                while (patch_delta > 0) {
                    column_offset += 1;
                    post_length -= 1;
                    patch_delta -= 1;
                }
            }
        }
    }
}

// ------------------------------------------------------------------------------------------------
void RenderSpan(u32* pixels, ActiveSpanData* span, int y_screen,
                const std::vector<doom::Flat>& flats, const core::Palette& palette,
                const std::vector<core::Colormap>& colormaps, const RenderData& render_data) {
    TimeFunction;

    // Advance `hit` to the current location
    span->hit += span->step * span->x_start;

    int pixels_index = y_screen * render_data.screen_size_x + span->x_start;
    for (int x_span = span->x_start; x_span < span->x_end; x_span++) {
        // Tile width is 1.0f
        int x_ind_hit = (int)(floorf(span->hit.x / 1.0f));
        int y_ind_hit = (int)(floorf(span->hit.y / 1.0f));
        f32 x_rem_hit = span->hit.x - 1.0f * x_ind_hit;
        f32 y_rem_hit = span->hit.y - 1.0f * y_ind_hit;
        u32 texture_x = (int)(x_rem_hit / 1.0f * 64.0f);  // Floor textures are always 64x64
        u32 texture_y = (int)(y_rem_hit / 1.0f * 64.0f);

        const auto& flat = flats[span->flat_id];
        u8 palette_index = flat.data[texture_x + 64 * texture_y];
        palette_index = colormaps[span->colormap_index].map[palette_index];
        u32 abgr = *(u32*)(palette.rgbs + 3 * palette_index) | 0xFF000000;
        pixels[pixels_index] = abgr;

        // Step
        span->hit += span->step;
        pixels_index += 1;
    }
}

// ------------------------------------------------------------------------------------------------
void Raycast(u32* pixels, f32* raycast_distances, ActiveSpanData* active_spans,
             const core::GameMap& game_map, const CameraState& camera,
             const std::vector<doom::Patch>& patches, const std::vector<doom::Flat>& flats,
             const core::Palette& palette, const std::vector<core::Colormap>& colormaps,
             const RenderData& render_data, common::Vec2f pos, common::Vec2f dir,
             core::QuarterEdgeIndex qe_dual, const int x,
             const f32 cam_len_times_screen_size_y_over_fov_y, common::Vec2f ray_dir_lo,
             common::Vec2f ray_dir_hi) {
    TimeFunction;

    const core::DelaunayMesh mesh = game_map.GetMesh();
    const int half_screen_size_y = render_data.half_screen_size_y;

    constexpr f32 kProjectionDistance = 100.0;  // Ridiculously large
    constexpr f32 kGetRightHandedNessEps = 1e-4f;

    common::Vec2f camera_dir = camera.dir;
    common::Vec2f camera_pos = camera.pos;
    f32 camera_z = camera.z;

    int y_lo = -1;
    int y_hi = render_data.screen_size_y;

    // TODO: Figure out graphics artifact when looking up stairs, likely due to patch rendering or
    // y_lo.

    bool recurse = true;
    u32 n_steps = 0;
    while (recurse && n_steps <= render_data.max_render_steps && y_lo < y_hi) {
        recurse = false;
        n_steps += 1;

        // Grab the enclosing triangle.
        auto [qe_ab, qe_bc, qe_ca] = mesh.GetTriangleQuarterEdges(qe_dual);

        common::Vec2f a = mesh.GetVertex(qe_ab);
        common::Vec2f b = mesh.GetVertex(qe_bc);
        common::Vec2f c = mesh.GetVertex(qe_ca);

        // Project our ray out far enough that it would exit our mesh
        common::Vec2f pos_next_delta = kProjectionDistance * dir;
        common::Vec2f pos_next = pos + pos_next_delta;

        f32 min_interp = INFINITY;
        core::QuarterEdgeIndex qe_side = {core::kInvalidIndex};

        // The edge vector of the face that we last crossed
        common::Vec2f v_face;

        // See if we cross any of the 3 faces for the triangle we are in,
        // and cross the first segment.
        {
            ProfileBlock sdl_poll_events("Raycast_EdgeChecks", __COUNTER__ + 1);

            if (common::GetRightHandedness(a, b, pos_next) < -kGetRightHandedNessEps) {
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
            if (common::GetRightHandedness(b, c, pos_next) < -kGetRightHandedNessEps) {
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
            if (common::GetRightHandedness(c, a, pos_next) < -kGetRightHandedNessEps) {
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
        }

        if (core::IsValid(qe_side)) {
            // Move to the face
            pos_next = pos + min_interp * pos_next_delta;
            qe_dual = mesh.Rot(qe_side);

            // Accumulate distance traveled
            raycast_distances[x] += min_interp * kProjectionDistance;

            const core::SideInfo* side_info = game_map.GetSideInfo(qe_side);
            if (side_info != nullptr) {
                ProfileBlock sdl_poll_events("Raycast_WithValidSideInfo", __COUNTER__ + 1);

                const f32 ray_len = std::max(raycast_distances[x], 0.01f);
                const f32 gamma = cam_len_times_screen_size_y_over_fov_y / ray_len;

                const core::Sector* sector = game_map.GetSector(side_info->sector_id);
                f32 z_ceil = sector->z_ceil;
                f32 z_upper = sector->z_ceil;
                f32 z_lower = sector->z_floor;
                f32 z_floor = sector->z_floor;

                // Get the height on the other side, if it is passable.
                const bool is_passable = ((side_info->flags & core::kSideInfoFlag_PASSABLE) > 0);
                bool is_portal = core::IsValid(side_info->qe_portal);
                if (is_passable && !is_portal) {
                    core::QuarterEdgeIndex qe_sym = mesh.Sym(qe_side);
                    const core::SideInfo* side_info_sym = game_map.GetSideInfo(qe_sym);

                    if (side_info_sym != nullptr) {
                        const core::Sector* sector_sym =
                            game_map.GetSector(side_info_sym->sector_id);
                        z_lower = sector_sym->z_floor;
                        z_upper = sector_sym->z_ceil;
                    } else {
                        std::cout << "Unexpected nullptr qe_sym!" << std::endl;
                    }
                }

                // TODO: Check how it is possible that y_floor ever be above the halfway point - we
                //       should never see such floors.
                int y_ceil = (int)(half_screen_size_y + gamma * (z_ceil - camera_z));
                int y_upper = (int)(half_screen_size_y + gamma * (z_upper - camera_z));
                int y_lower = (int)(half_screen_size_y + gamma * (z_lower - camera_z));
                int y_floor = (int)(half_screen_size_y + gamma * (z_floor - camera_z));

                // Calculate where along the segment we intersected
                f32 v_face_len = common::Norm(v_face);
                f32 x_along_texture = v_face_len - common::Norm(pos_next - mesh.GetVertex(qe_side));

                // Determine the light level
                u8 colormap_index =
                    sector->light_level + (u8)(render_data.darkness_per_world_dist * ray_len);

                // Make faces that run closer to north-south brighter, and faces running closer
                // to east-west darker. (cos > 0.7071). We're using the law of cosines.
                bool face_is_closer_to_east_west = std::abs(v_face.x / v_face_len) > 0.7071;
                if (face_is_closer_to_east_west) {
                    colormap_index += 1;  // darker
                } else if (colormap_index > 0) {
                    colormap_index -= 1;  // lighter
                }

                u8 max_colormap_index = 32 - 1;
                colormap_index = std::min(colormap_index, max_colormap_index);
                const core::Colormap& colormap = colormaps[colormap_index];

                // Render the ceiling above the upper texture
                y_ceil = std::max(y_ceil, y_lo);
                while (y_hi > y_ceil) {
                    y_hi--;
                    if (y_hi < half_screen_size_y) {
                        y_hi = y_ceil;
                        break;
                    }

                    ActiveSpanData& active_span = active_spans[y_hi];

                    // Render the span if we close it out - that is, the sector is different than
                    // what it was before
                    bool start_new_span = !active_span.is_active ||
                                          active_span.sector_id != side_info->sector_id ||
                                          active_span.x_end < x - 1;

                    if (start_new_span && active_span.is_active) {
                        // Render the span
                        RenderSpan(pixels, &active_span, y_hi, flats, palette, colormaps,
                                   render_data);
                    }

                    if (start_new_span) {
                        // Start a new span
                        active_span.is_active = true;
                        active_span.x_start = x;
                        active_span.x_end = x;
                        active_span.sector_id = side_info->sector_id;
                        active_span.flat_id = sector->flat_id_ceil;

                        // @efficiency: precompute
                        f32 zpp = (y_hi - render_data.screen_size_y / 2.0f) *
                                  (camera.fov.y / render_data.screen_size_y);

                        // distance of the ray from the player through y_hi at the ceiling.
                        f32 radius = (z_ceil - camera_z) / zpp;
                        radius = std::max(radius, 0.01f);

                        active_span.colormap_index =
                            sector->light_level +
                            (u8)(render_data.darkness_per_world_dist * radius);
                        active_span.colormap_index =
                            std::min(active_span.colormap_index, max_colormap_index);

                        // Location of the 1st ray's intersection
                        active_span.hit.x = camera_pos.x + radius * ray_dir_lo.x;
                        active_span.hit.y = camera_pos.y + radius * ray_dir_lo.y;

                        // Each step is
                        // (hit_x2 - hit_x) / SCREEN_SIZE_X;
                        // = ((camera->pos.x + radius * ray_dir_lo_x) - (camera->pos.x + radius *
                        // ray_dir_lo_x)) / SCREEN_SIZE_X = (radius * ray_dir_lo_x - (radius *
                        // ray_dir_lo_x)) / SCREEN_SIZE_X = radius * (ray_dir_hi_x - ray_dir_lo_x) /
                        // SCREEN_SIZE_X
                        // @efficiency - could precompute the delta divided by sceen size.
                        active_span.step.x =
                            radius * (ray_dir_hi.x - ray_dir_lo.x) / render_data.screen_size_x;
                        active_span.step.y =
                            radius * (ray_dir_hi.y - ray_dir_lo.y) / render_data.screen_size_x;
                    } else {
                        active_span.x_end = x;  // @efficiency - need to do this either way
                    }
                }

                // Render the upper texture
                if (y_hi > y_upper) {
                    f32 column_height = z_ceil - z_upper;
                    int y_upper_clamped = std::max(y_upper, y_lo);
                    RenderPatchColumn(pixels, x, y_upper, y_ceil, y_upper_clamped, y_hi,
                                      x_along_texture, side_info->texture_info_upper.x_offset,
                                      side_info->texture_info_upper.y_offset, column_height,
                                      patches[side_info->texture_info_upper.texture_id], palette,
                                      colormap, render_data);
                    y_hi = y_upper_clamped;
                }

                // Render the floor below the lower texture
                y_floor = std::min(y_floor, y_hi);
                while (y_floor > y_lo) {
                    y_lo++;
                    if (y_lo >= half_screen_size_y) {
                        y_lo = y_floor;
                        break;
                    }

                    ActiveSpanData& active_span = active_spans[y_lo];

                    // Render the span if we close it out - that is, the sector is different than
                    // what it was before
                    bool start_new_span = !active_span.is_active ||
                                          active_span.sector_id != side_info->sector_id ||
                                          active_span.x_end < x - 1;

                    if (start_new_span && active_span.is_active) {
                        // Render the span
                        RenderSpan(pixels, &active_span, y_lo, flats, palette, colormaps,
                                   render_data);
                    }

                    if (start_new_span) {
                        // Start a new span
                        active_span.is_active = true;
                        active_span.x_start = x;
                        active_span.x_end = x;
                        active_span.sector_id = side_info->sector_id;
                        active_span.flat_id = sector->flat_id_floor;

                        // @efficiency: precompute
                        f32 zpp = (render_data.screen_size_y / 2.0f - y_lo) *
                                  (camera.fov.y / render_data.screen_size_y);

                        // distance of the ray from the player through y_lo at the floor.
                        f32 radius = (camera_z - z_floor) / zpp;
                        radius = std::max(radius, 0.01f);

                        active_span.colormap_index =
                            sector->light_level +
                            (u8)(render_data.darkness_per_world_dist * radius);
                        active_span.colormap_index =
                            std::min(active_span.colormap_index, max_colormap_index);

                        // Location of the 1st ray's intersection
                        active_span.hit.x = camera_pos.x + radius * ray_dir_lo.x;
                        active_span.hit.y = camera_pos.y + radius * ray_dir_lo.y;

                        // Each step is(hit_x2 - hit_x) / SCREEN_SIZE_X;
                        // = ((camera->pos.x + radius * ray_dir_lo_x) - (camera->pos.x + radius *
                        // ray_dir_lo_x)) / SCREEN_SIZE_X = (radius * ray_dir_lo_x - (radius *
                        // ray_dir_lo_x)) / SCREEN_SIZE_X = radius * (ray_dir_hi_x - ray_dir_lo_x) /
                        // SCREEN_SIZE_X
                        // @efficiency - could precompute the delta divided by screen size.
                        active_span.step.x =
                            radius * (ray_dir_hi.x - ray_dir_lo.x) / render_data.screen_size_x;
                        active_span.step.y =
                            radius * (ray_dir_hi.y - ray_dir_lo.y) / render_data.screen_size_x;
                    } else {
                        active_span.x_end = x;  // @efficiency - need to do this either way
                    }
                }

                // Render the lower texture
                if (y_lower > y_lo) {
                    f32 column_height = z_lower - z_floor;
                    int y_lower_clamped = std::min(y_lower, y_hi);
                    RenderPatchColumn(pixels, x, y_floor, y_lower, y_lo, y_lower_clamped,
                                      x_along_texture, side_info->texture_info_lower.x_offset,
                                      side_info->texture_info_lower.y_offset, column_height,
                                      patches[side_info->texture_info_lower.texture_id], palette,
                                      colormap, render_data);
                    y_lo = y_lower_clamped;
                }

                // Recurse if the side is passable.
                if (is_passable) {
                    common::Vec2f dir_next = dir;

                    if (is_portal) {
                        qe_dual = mesh.Tor(side_info->qe_portal);

                        // Do the transform on pos_next and dir

                        // The vector along the new face is
                        common::Vec2f c = mesh.GetVertex(side_info->qe_portal);
                        common::Vec2f d = mesh.GetVertex(mesh.Sym(side_info->qe_portal));
                        common::Vec2f dc = d - c;
                        f32 dc_len = common::Norm(dc);

                        // Our new position is:
                        pos_next = c + (x_along_texture / dc_len) * dc;

                        // A -> B is v_face (and v_face_len)
                        // P -> Q (our ray) is dir
                        f32 cos_theta = common::Dot(v_face, dir) / v_face_len;
                        f32 sin_theta = abs(common::Cross(v_face, dir) / v_face_len);

                        // Now to get the new direction, we need to rotate out of the new face.
                        dir_next = {(-cos_theta * dc.x - sin_theta * dc.y) / dc_len,
                                    (sin_theta * dc.x - cos_theta * dc.y) / dc_len};

                        // Update the camera position to be where it would be, as if we were on the
                        // other side of the destination portal looking through it.
                        camera_pos = pos_next - dir_next * ray_len;

                        // Need to also update our ray directions
                        // Doing this correctly requires rotating the original camera.
                        // We don't want dir_next, but rather the overall camera's dir.
                        // @efficiency - precompute all viewing column angles.
                        f32 column_angle_run = ((int)(render_data.screen_size_x / 2) - x) *
                                               camera.fov.x / render_data.screen_size_x;
                        f32 column_angle_hypot = std::hypot(column_angle_run, 1.0f);
                        f32 sin_pix_column_angle = column_angle_run / column_angle_hypot;
                        f32 cos_pix_column_angle = 1.0f / column_angle_hypot;

                        // Back out our camera direction by rotating back by the viewing column
                        // angle.
                        camera_dir = Rot(dir_next, -sin_pix_column_angle, cos_pix_column_angle);

                        // raycast direction for leftmost pixel column, x = 0
                        f32 half_camera_width = camera.fov.x / 2.0f;
                        ray_dir_lo = {camera_dir.x - half_camera_width * camera_dir.y,
                                      camera_dir.y + half_camera_width * camera_dir.x};
                        // raycast direction for rightmost pixel column, x = screen_size_x
                        ray_dir_hi = {camera_dir.x + half_camera_width * camera_dir.y,
                                      camera_dir.y - half_camera_width * camera_dir.x};

                        // The camera height on the other side.
                        const core::SideInfo* side_info_portal =
                            game_map.GetSideInfo(side_info->qe_portal);
                        const core::Sector* sector_side = game_map.GetSector(side_info->sector_id);
                        const core::Sector* sector_portal =
                            game_map.GetSector(side_info_portal->sector_id);
                        camera_z += (sector_portal->z_floor - sector_side->z_floor);
                    }

                    // Recurse.
                    recurse = true;
                    dir = dir_next;
                    pos = pos_next;
                } else {
                    // The side info has a solid wall.
                    f32 column_height = z_upper - z_lower;
                    RenderPatchColumn(pixels, x, y_lower, y_upper, y_lo, y_hi, x_along_texture,
                                      side_info->texture_info_middle.x_offset,
                                      side_info->texture_info_middle.y_offset, column_height,
                                      patches[side_info->texture_info_middle.texture_id], palette,
                                      colormap, render_data);
                }
            } else {
                // Failed to get side_info.
                // If this is a normal edge, recurse as long as it is not the boundary edge
                recurse = !mesh.IsBoundaryEdge(qe_side);
                pos = pos_next;
            }
        }
    }
}

// ------------------------------------------------------------------------------------------------
void RenderScene(u32* pixels, f32* raycast_distances, f32* dw, ActiveSpanData* active_spans,
                 const core::GameMap& game_map, const CameraState& camera,
                 const core::RenderAssets& render_assets, const RenderData& render_data) {
    TimeFunction;

    // Unpack
    const std::vector<doom::Patch>& patches = render_assets.patches;
    const std::vector<doom::Flat>& flats = render_assets.flats;
    const core::Palette& palette = render_assets.palettes[0];
    const std::vector<core::Colormap>& colormaps = render_assets.colormaps;

    // Zero out the raycast distances
    memset(raycast_distances, 0, sizeof(f32) * render_data.screen_size_x);
    memset(active_spans, 0, sizeof(ActiveSpanData) * render_data.screen_size_y);

    // Camera data
    if (!core::IsValid(camera.qe)) {
        return;
    }

    f32 screen_size_y_over_fov_y = render_data.screen_size_y / camera.fov.y;

    f32 half_camera_width = camera.fov.x / 2.0f;

    // raycast direction for leftmost pixel column, x = 0
    common::Vec2f ray_dir_lo = {camera.dir.x - half_camera_width * camera.dir.y,
                                camera.dir.y + half_camera_width * camera.dir.x};
    // raycast direction for rightmost pixel column, x = screen_size_x
    common::Vec2f ray_dir_hi = {camera.dir.x + half_camera_width * camera.dir.y,
                                camera.dir.y - half_camera_width * camera.dir.x};

    for (u32 x = 0; x < render_data.screen_size_x; x++) {
        // Camera to pixel column
        const f32 delta_w = dw[x];
        const common::Vec2f cp = {camera.dir.x - delta_w * camera.dir.y,
                                  camera.dir.y + delta_w * camera.dir.x};

        // Distance from the camera to the column
        const f32 cam_len = common::Norm(cp);
        const f32 cam_len_times_screen_size_y_over_fov_y = cam_len * screen_size_y_over_fov_y;

        // Ray direction through this column
        const common::Vec2f dir = cp / cam_len;

        // Step through triangles until we hit a solid triangle
        Raycast(pixels, raycast_distances, active_spans, game_map, camera, patches, flats, palette,
                colormaps, render_data, camera.pos, dir, camera.qe, x,
                cam_len_times_screen_size_y_over_fov_y, ray_dir_lo, ray_dir_hi);
    }

    // Render all remaining floor and ceiling spans
    for (int y = 0; y < (int)(render_data.screen_size_y); y++) {
        auto& active_span = active_spans[y];
        if (!active_span.is_active) {
            continue;
        }
        RenderSpan(pixels, &active_span, y, flats, palette, colormaps, render_data);
    }
}

// ------------------------------------------------------------------------------------------------
void ImportGameData(core::GameMap* map, core::RenderAssets* render_assets) {
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Importing game data" << std::endl;

    core::AssetsExporter exporter;
    exporter.LoadAssetsFile("../toom/assets/toomed.bin");
    std::cout << "Num entries:" << exporter.NumEntries() << std::endl;

    bool succeeded = core::ImportPalettes(&(render_assets->palettes), exporter);
    if (!succeeded) {
        std::cout << "Failed to load palettes!" << std::endl;
    }

    succeeded = core::ImportColormaps(&(render_assets->colormaps), exporter);
    if (!succeeded) {
        std::cout << "Failed to load colormaps!" << std::endl;
    }

    succeeded = doom::ImportPatches(&(render_assets->patches), exporter);
    if (!succeeded) {
        std::cout << "Failed to load patches!" << std::endl;
    }

    succeeded = doom::ImportFlats(&(render_assets->flats), exporter);
    if (!succeeded) {
        std::cout << "Failed to load flats!" << std::endl;
    }

    succeeded = map->Import(exporter);
    if (!succeeded) {
        std::cout << "Failed to load game map! Clearing possibly-corrupted map." << std::endl;
        map->Clear();
    }

    std::cout << "DONE" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
}

// ------------------------------------------------------------------------------------------------
void ExportGameData(const core::GameMap& map, const core::RenderAssets& render_assets) {
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Exporting game data" << std::endl;

    core::AssetsExporter exporter;
    std::cout << "Exporting core data" << std::endl;
    exporter.AddEntry(core::ExportPalettes(render_assets.palettes));
    exporter.AddEntry(core::ExportColormaps(render_assets.colormaps));
    exporter.AddEntry(doom::ExportPatches(render_assets.patches));
    exporter.AddEntry(doom::ExportFlats(render_assets.flats));

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

    // Load our DOOM assets
    std::unique_ptr<core::WadImporter> doom_assets =
        core::WadImporter::LoadFromFile("../toom/assets/DOOM.WAD");
    ASSERT(doom_assets, "Failed to load DOOM Assets");

    // Import our julia assets
    // std::unique_ptr<core::AssetsImporter> julia_assets =
    //     core::AssetsImporter::LoadFromFile("../toom/assets/assets.bin");
    // ASSERT(julia_assets, "Failed to load Julia Assets");

    // Extract our wall textures
    // auto bitmap_data_opt = julia_assets->FindEntryData("textures");
    // ASSERT(bitmap_data_opt, "Failed to find texture data");
    // core::OldStyleBitmap bitmap = core::LoadBitmap(*bitmap_data_opt);
    // for (const doom::Patch& patch :
    //      core::ExtractPatches(bitmap, "julia", render_assets.palettes[0])) {
    //     render_assets.patches.push_back(patch);
    // }

    // Create our render assets
    core::RenderAssets render_assets;
    // render_assets.palettes = doom::ParseDoomPalettes(doom_assets);
    // ASSERT(render_assets.palettes.size() == doom::kNumDoomPalettes, "Failed to load DOOM
    // palettes"); render_assets.colormaps = doom::ParseDoomColormaps(doom_assets);
    // ASSERT(render_assets.colormaps.size() == doom::kNumDoomColormaps,
    //        "Failed to load DOOM colormaps");
    // render_assets.patches = doom::ParseDoomTextures(doom_assets);
    // ASSERT(render_assets.patches.size() > 0, "Failed to load DOOM patches");
    // render_assets.flats = doom::ParseDoomFlats(doom_assets);
    // ASSERT(render_assets.flats.size() > 0, "Failed to load DOOM flats");

    // Create our map
    core::GameMap map;

    // Import our game data
    ImportGameData(&map, &render_assets);

    // Generate some basic flats for me
    // {
    //     doom::Flat flat;
    //     flat.name = "mono_92";
    //     memset(flat.data, 92, sizeof(flat.data));
    //     render_assets.flats.push_back(flat);
    // }
    // {
    //     doom::Flat flat;
    //     flat.name = "mono_100";
    //     memset(flat.data, 100, sizeof(flat.data));
    //     render_assets.flats.push_back(flat);
    // }
    // {
    //     doom::Flat flat;
    //     flat.name = "mono_110";
    //     memset(flat.data, 110, sizeof(flat.data));
    //     render_assets.flats.push_back(flat);
    // }

    // Load a doom level
    // LoadDoomLevel(&map, "E1M1", doom_assets, render_assets);

    // Export mesh for Julia
    // {
    //     // First, print out the mesh vertices
    //     const core::DelaunayMesh& mesh = map.GetMesh();

    //     printf("vertices = Vector{Float64}[\n");
    //     std::map<core::VertexIndex, int> vertex_to_index;  // 0 - indexed
    //     core::VertexIndex i_vertex = mesh.GetFirstVertexIndex();
    //     while (core::IsValid(i_vertex)) {
    //         auto v = mesh.GetVertex(i_vertex);
    //         vertex_to_index[i_vertex] = vertex_to_index.size();
    //         printf("\t[%.3f, %.3f],\n", v.x, v.y);
    //         i_vertex = mesh.GetNext(i_vertex);
    //     }
    //     printf("];\n");

    //     printf("edges = Tuple{Int, Int}[\n");
    //     int n_printed = 0;
    //     core::QuarterEdgeIndex qe = mesh.GetFirstQuarterEdgeIndex();
    //     while (core::IsValid(qe)) {
    //         if (mesh.IsPrimal(qe)) {
    //             // Get its opposite side.
    //             core::QuarterEdgeIndex qe_sym = mesh.Sym(qe);

    //             const core::VertexData& a = mesh.GetVertexData(qe);
    //             const core::VertexData& b = mesh.GetVertexData(qe_sym);
    //             if (a.i_self > b.i_self) {  // Avoid rendering edges twice
    //                 printf("\t(%d, %d),", vertex_to_index[a.i_self], vertex_to_index[b.i_self]);
    //                 n_printed += 1;
    //                 if (n_printed % 5 == 0) {
    //                     printf("\n");
    //                 }
    //             }
    //         }
    //         // Get the next one
    //         qe = mesh.GetNext(qe);
    //     }
    //     printf("];\n");

    //     // Print out whether each edge is traversible or not.
    //     // An edge is traversible if it is not constrained or if it is a sideinfo that is
    //     // traversible where the sector sides don't change in z more than a given threshold.
    //     constexpr f32 kMaxTraversibleDeltaZ = 0.5;
    //     printf("edge_traversibility = BitVector([\n");
    //     qe = mesh.GetFirstQuarterEdgeIndex();
    //     n_printed = 0;
    //     while (core::IsValid(qe)) {
    //         if (mesh.IsPrimal(qe)) {
    //             // Get its opposite side.
    //             core::QuarterEdgeIndex qe_sym = mesh.Sym(qe);

    //             const core::VertexData& a = mesh.GetVertexData(qe);
    //             const core::VertexData& b = mesh.GetVertexData(qe_sym);
    //             if (a.i_self > b.i_self) {
    //                 // This is a new edge. Determine traversibility.
    //                 bool is_traversible = true;

    //                 if (mesh.IsConstrained(qe)) {
    //                     is_traversible = false;

    //                     // Check for a SideInfo.
    //                     const core::SideInfo* sideinfo_qe = map.GetSideInfo(qe);
    //                     const core::SideInfo* sideinfo_qe_sym = map.GetSideInfo(qe_sym);

    //                     if (sideinfo_qe != nullptr && sideinfo_qe_sym != nullptr) {
    //                         // Potentially traversible. Check the passable flag.
    //                         if ((sideinfo_qe->flags & core::kSideInfoFlag_PASSABLE) > 0) {
    //                             // Finally, verify the sector height change.
    //                             const core::Sector* sector =
    //                             map.GetSector(sideinfo_qe->sector_id); f32 z1 = sector->z_floor;
    //                             const core::Sector* sector_sym =
    //                                 map.GetSector(sideinfo_qe_sym->sector_id);
    //                             f32 z2 = sector_sym->z_floor;
    //                             if (std::abs(z2 - z1) < kMaxTraversibleDeltaZ) {
    //                                 is_traversible = true;
    //                             }
    //                         }
    //                     }
    //                 }

    //                 printf("\t%s,", is_traversible ? "true" : "false");
    //                 n_printed += 1;
    //                 if (n_printed % 5 == 0) {
    //                     printf("\n");
    //                 }
    //             }
    //         }
    //         // Get the next one
    //         qe = mesh.GetNext(qe);
    //     }
    //     printf("]);\n");

    //     // Print out the centers of every triangle.
    //     printf("dual_vertices = Vector{Float64}[\n");
    //     std::map<core::QuarterEdgeIndex, int> dual_qe_to_index;  // 1 - indexed
    //     qe = mesh.GetFirstQuarterEdgeIndex();
    //     while (core::IsValid(qe)) {
    //         if (mesh.IsDual(qe)) {
    //             // Check to make sure this quarter edge is the lowest among the ones that share
    //             its
    //             // origin.
    //             if (qe < mesh.Next(qe) && qe < mesh.Prev(qe)) {
    //                 const auto [qe_ab, qe_bc, qe_ca] = mesh.GetTriangleQuarterEdges(qe);
    //                 const common::Vec2f& a = mesh.GetVertex(qe_ab);
    //                 const common::Vec2f& b = mesh.GetVertex(qe_bc);
    //                 const common::Vec2f& c = mesh.GetVertex(qe_ca);
    //                 const common::Vec2f v = (a + b + c) / 3.0f;

    //                 printf("\t[%.3f, %.3f],\n", v.x, v.y);
    //                 dual_qe_to_index[qe] = dual_qe_to_index.size();
    //             }
    //         }
    //         // Get the next one
    //         qe = mesh.GetNext(qe);
    //     }
    //     printf("];\n");

    //     // Print out the traversible dual edges
    //     printf("dual_edges = Tuple{Int, Int}[\n");
    //     qe = mesh.GetFirstQuarterEdgeIndex();
    //     while (core::IsValid(qe)) {
    //         if (mesh.IsDual(qe)) {
    //             // Get its opposite side.
    //             core::QuarterEdgeIndex qe_sym = mesh.Sym(qe);

    //             core::QuarterEdgeIndex qe_a = std::min(std::min(qe, mesh.Next(qe)),
    //             mesh.Prev(qe)); core::QuarterEdgeIndex qe_b =
    //                 std::min(std::min(qe_sym, mesh.Next(qe_sym)), mesh.Prev(qe_sym));

    //             // Check to see whether the edge is traversible.
    //             const core::VertexData& a = mesh.GetVertexData(mesh.Rot(qe));
    //             const core::VertexData& b = mesh.GetVertexData(mesh.Rot(qe_sym));
    //             if (a.i_self > b.i_self) {
    //                 // This is a new edge. Determine traversibility.
    //                 bool is_traversible = true;

    //                 if (mesh.IsConstrained(mesh.Rot(qe))) {
    //                     is_traversible = false;
    //                     // Check for a SideInfo.
    //                     const core::SideInfo* sideinfo_qe = map.GetSideInfo(mesh.Rot(qe));
    //                     const core::SideInfo* sideinfo_qe_sym =
    //                     map.GetSideInfo(mesh.Rot(qe_sym));

    //                     if (sideinfo_qe != nullptr && sideinfo_qe_sym != nullptr) {
    //                         // Potentially traversible. Check the passable flag.
    //                         if ((sideinfo_qe->flags & core::kSideInfoFlag_PASSABLE) > 0) {
    //                             // Finally, verify the sector height change.
    //                             const core::Sector* sector =
    //                             map.GetSector(sideinfo_qe->sector_id); f32 z1 = sector->z_floor;
    //                             const core::Sector* sector_sym =
    //                                 map.GetSector(sideinfo_qe_sym->sector_id);
    //                             f32 z2 = sector_sym->z_floor;
    //                             if (std::abs(z2 - z1) < kMaxTraversibleDeltaZ) {
    //                                 is_traversible = true;
    //                             }
    //                         }
    //                     }
    //                 }

    //                 if (is_traversible) {
    //                     printf("\t(%d, %d),\n", dual_qe_to_index[qe_a], dual_qe_to_index[qe_b]);
    //                 }
    //             }
    //         }
    //         // Get the next one
    //         qe = mesh.GetNext(qe);
    //     }
    //     printf("];\n");

    //     // Print out which primal edge was traversed for a dual edge
    //     printf("primal_edges = Tuple{Int, Int}[\n");
    //     qe = mesh.GetFirstQuarterEdgeIndex();
    //     while (core::IsValid(qe)) {
    //         if (mesh.IsDual(qe)) {
    //             // Get its opposite side.
    //             core::QuarterEdgeIndex qe_sym = mesh.Sym(qe);

    //             core::QuarterEdgeIndex qe_a = std::min(std::min(qe, mesh.Next(qe)),
    //             mesh.Prev(qe)); core::QuarterEdgeIndex qe_b =
    //                 std::min(std::min(qe_sym, mesh.Next(qe_sym)), mesh.Prev(qe_sym));

    //             // Check to see whether the edge is traversible.
    //             const core::VertexData& a = mesh.GetVertexData(mesh.Rot(qe));
    //             const core::VertexData& b = mesh.GetVertexData(mesh.Rot(qe_sym));
    //             if (a.i_self > b.i_self) {
    //                 // This is a new edge. Determine traversibility.
    //                 bool is_traversible = true;

    //                 if (mesh.IsConstrained(mesh.Rot(qe))) {
    //                     is_traversible = false;
    //                     // Check for a SideInfo.
    //                     const core::SideInfo* sideinfo_qe = map.GetSideInfo(mesh.Rot(qe));
    //                     const core::SideInfo* sideinfo_qe_sym =
    //                     map.GetSideInfo(mesh.Rot(qe_sym));

    //                     if (sideinfo_qe != nullptr && sideinfo_qe_sym != nullptr) {
    //                         // Potentially traversible. Check the passable flag.
    //                         if ((sideinfo_qe->flags & core::kSideInfoFlag_PASSABLE) > 0) {
    //                             // Finally, verify the sector height change.
    //                             const core::Sector* sector =
    //                             map.GetSector(sideinfo_qe->sector_id); f32 z1 = sector->z_floor;
    //                             const core::Sector* sector_sym =
    //                                 map.GetSector(sideinfo_qe_sym->sector_id);
    //                             f32 z2 = sector_sym->z_floor;
    //                             if (std::abs(z2 - z1) < kMaxTraversibleDeltaZ) {
    //                                 is_traversible = true;
    //                             }
    //                         }
    //                     }
    //                 }

    //                 if (is_traversible) {
    //                     printf("\t(%d, %d),\n", vertex_to_index[a.i_self],
    //                            vertex_to_index[b.i_self]);
    //                 }
    //             }
    //         }
    //         // Get the next one
    //         qe = mesh.GetNext(qe);
    //     }
    //     printf("];\n");
    // }

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
    core::KeyBoardState keyboard_state;
    core::ClearKeyboardState(&keyboard_state);

    CameraState player_cam = {};
    player_cam.pos = {5.0, 5.0};
    // player_cam.pos = {-10.25, -50.5};
    // player_cam.pos = {-4.9, -50.5};
    player_cam.dir = {1.0, 0.0};
    player_cam.fov = {1.5, 0.84375};
    player_cam.height = 49.5f / 64.0f;  // 0.4f;
    player_cam.z = player_cam.height;
    player_cam.qe = map.GetMesh().GetEnclosingTriangle(player_cam.pos);

    RenderData render_data;
    render_data.n_pixels_per_world_unit = 64;
    render_data.screen_size_x = player_window_data.screen_size_x;
    render_data.screen_size_y = player_window_data.screen_size_y;
    render_data.half_screen_size_y = render_data.screen_size_y >> 1;
    render_data.darkness_per_world_dist = 1.5f;
    render_data.max_render_steps = 64;
    render_data.color_ceil = 0xFF222222;
    render_data.color_floor = 0xFF444444;

    // Player view data
    u32 player_view_pixels[player_window_data.screen_size_x *
                           player_window_data.screen_size_y];  // row-major
    f32 raycast_distances[player_window_data.screen_size_x];
    ActiveSpanData active_spans[render_data.screen_size_y];

    // This precomputation assumes that the fov does not change
    f32 dw[player_window_data.screen_size_x];
    for (int x = 0; x < player_window_data.screen_size_x; x++) {
        dw[x] = player_cam.fov.x / 2.0f - (player_cam.fov.x * x) / render_data.screen_size_x;
    }

    // Estimate our CPU frequency
    u64 cpu_freq = EstimateCPUTimerFreq(100);

    bool continue_running = true;
    while (continue_running) {
        BeginProfile();

        {
            // Poll and handle events (inputs, window resize, etc.)
            // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear
            // imgui wants to use your inputs.
            // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main
            //   application, or clear/overwrite your copy of the mouse data.
            // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your
            // main
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
                            mouse_click_pos =
                                CameraToGlobal(common::Vec2f(event.button.x, event.button.y),
                                               camera_pos, camera_zoom);
                            camera_pos_at_mouse_click = camera_pos;

                            // Check for a selected vertex near the current mouse position
                            selected_edge_index = {core::kInvalidIndex};
                            selected_vertex_index = map.FindVertexNearPosition(
                                mouse_click_pos, qe_mouse_face, 0.3f / camera_zoom * 50.0f);

                            // Check for a selected edge near the current mouse position if we did
                            // not select a vertex
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
                            const bool did_not_drag =
                                common::Norm(mouse_pos - mouse_click_pos) < 0.2;
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
                    if (event.key.keysym.sym == SDLK_o) {
                        ExportGameData(map, render_assets);
                    } else if (event.key.keysym.sym == SDLK_i) {
                        ImportGameData(&map, &render_assets);
                        player_cam.qe = map.GetMesh().GetEnclosingTriangle(
                            player_cam.pos);  // TODO: Move somewhere better
                    } else if (event.key.keysym.sym == SDLK_w) {
                        keyboard_state.w = core::KeyboardKeyState_Pressed;
                    } else if (event.key.keysym.sym == SDLK_s) {
                        keyboard_state.s = core::KeyboardKeyState_Pressed;
                    } else if (event.key.keysym.sym == SDLK_d) {
                        keyboard_state.d = core::KeyboardKeyState_Pressed;
                    } else if (event.key.keysym.sym == SDLK_a) {
                        keyboard_state.a = core::KeyboardKeyState_Pressed;
                    } else if (event.key.keysym.sym == SDLK_q) {
                        keyboard_state.q = core::KeyboardKeyState_Pressed;
                    } else if (event.key.keysym.sym == SDLK_e) {
                        keyboard_state.e = core::KeyboardKeyState_Pressed;
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
                } else if (event.type == SDL_KEYUP && !io.WantCaptureKeyboard) {
                    if (event.key.keysym.sym == SDLK_w) {
                        keyboard_state.w = core::KeyboardKeyState_Released;
                    } else if (event.key.keysym.sym == SDLK_s) {
                        keyboard_state.s = core::KeyboardKeyState_Released;
                    } else if (event.key.keysym.sym == SDLK_d) {
                        keyboard_state.d = core::KeyboardKeyState_Released;
                    } else if (event.key.keysym.sym == SDLK_a) {
                        keyboard_state.a = core::KeyboardKeyState_Released;
                    } else if (event.key.keysym.sym == SDLK_q) {
                        keyboard_state.q = core::KeyboardKeyState_Released;
                    } else if (event.key.keysym.sym == SDLK_e) {
                        keyboard_state.e = core::KeyboardKeyState_Released;
                    }
                }
            }
        }

        u32 color_white = 0xFFFFFFFF;
        u32 color_background = 0x414141FF;
        u32 color_light_background = 0x454545FF;
        u32 color_light_gray = 0x909090FF;
        u32 color_qe_constrained = 0x52BCC4FF;
        u32 color_qe_normal = 0x587F82FF;
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
                common::Vec2f d = c + Rotr(Normalize(b - a)) * 0.1;
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
            common::Vec2f c = (a + b) / 2.0;
            common::Vec2f d = c + Rotr(Normalize(b - a)) * 0.1;
            auto c_cam = GlobalToCamera(c, camera_pos, camera_zoom);
            auto d_cam = GlobalToCamera(d, camera_pos, camera_zoom);
            SDL_RenderDrawLine(renderer, (int)(c_cam.x), (int)(c_cam.y), (int)(d_cam.x),
                               (int)(d_cam.y));
        }

        if (IsValid(selected_vertex_index)) {
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

        {  // Render the player camera position
            auto renderer = editor_window_data.renderer;

            auto dir_rotr = common::Rotr(player_cam.dir);
            auto a_cam = GlobalToCamera(player_cam.pos, camera_pos, camera_zoom);
            auto b_cam = GlobalToCamera(player_cam.pos + 0.5 * player_cam.dir - 0.4 * dir_rotr,
                                        camera_pos, camera_zoom);
            auto c_cam = GlobalToCamera(player_cam.pos + 0.5 * player_cam.dir + 0.4 * dir_rotr,
                                        camera_pos, camera_zoom);

            SDL_Vertex triangle[3];
            triangle[0] = {
                SDL_FPoint{a_cam.x, a_cam.y},
                SDL_Color{155, 55, 55, 155},
                SDL_FPoint{0},
            };
            triangle[1] = {
                SDL_FPoint{b_cam.x, b_cam.y},
                SDL_Color{155, 55, 55, 155},
                SDL_FPoint{0},
            };
            triangle[2] = {
                SDL_FPoint{c_cam.x, c_cam.y},
                SDL_Color{155, 55, 55, 155},
                SDL_FPoint{0},
            };
            SDL_RenderGeometry(renderer, nullptr, triangle, 3, nullptr, 0);
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
                const core::DelaunayMesh& mesh = map.GetMesh();

                // Vertex panel
                vertex_panel_drawn_last_frame = true;
                ImGui::Begin("Vertex");
                ImGui::Text("qe index:      %lu", selected_vertex_index.i);
                ImGui::Text("vx index:      %lu",
                            mesh.GetVertexData(selected_vertex_index).i_self.i);
                ImGui::Separator();

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

        {
            ImGui::Begin("EnforceEdge");
            const core::DelaunayMesh& mesh = map.GetMesh();

            static int int_src = 0;
            static int int_dst = 0;
            ImGui::InputInt("qe src: ", &int_src);
            ImGui::InputInt("qe dst: ", &int_dst);

            static bool success = true;
            if (ImGui::Button("enforce edge")) {
                core::QuarterEdgeIndex qe_src = {int_src};
                core::QuarterEdgeIndex qe_dst = {int_dst};

                if (mesh.IsPrimal(qe_src)) {
                    success = map.MaybeEnforceEdge(qe_src, qe_dst);
                } else {
                    success = false;
                }
            }
            ImGui::Text("success: %s", success ? "SUCCESS" : "FAILURE");

            ImGui::End();
        }

        // Side Info panel
        if (core::IsValid(selected_edge_index)) {
            const core::DelaunayMesh& mesh = map.GetMesh();

            ImGui::Begin("SideInfo");

            // Present the option to switch side infos
            ImGui::Text("index:      %lu", selected_edge_index.i);
            if (ImGui::Button("prev")) {
                selected_edge_index = mesh.Prev(selected_edge_index);
            }
            ImGui::SameLine();
            if (ImGui::Button("next")) {
                selected_edge_index = mesh.Next(selected_edge_index);
            }
            ImGui::SameLine();
            if (ImGui::Button("sym")) {
                selected_edge_index = mesh.Sym(selected_edge_index);
            }

            ImGui::Separator();
            core::SideInfo* side_info = map.GetEditableSideInfo(selected_edge_index);
            if (side_info != nullptr) {
                if (ImGui::Button(((side_info->flags & core::kSideInfoFlag_PASSABLE) > 0
                                       ? "Passable"
                                       : "Not Passable"))) {
                    side_info->flags ^= core::kSideInfoFlag_PASSABLE;  // toggle

                    // If we set it to passable, ensure that the opposite side is also a passable
                    // side info.
                    if (side_info->flags & core::kSideInfoFlag_PASSABLE) {
                        core::QuarterEdgeIndex qe_sym = mesh.Sym(selected_edge_index);
                        map.AddSideInfo(qe_sym);
                        map.GetEditableSideInfo(qe_sym)->flags |= core::kSideInfoFlag_PASSABLE;
                    }
                }

                ImGui::Text("flags:      %X", side_info->flags);

                int flags = 0;
                u16 step_u16 = 1;
                i16 step_i16 = 1;
                f32 step_f32 = 0.1f;

                ImGui::Separator();
                ImGui::Text(
                    "upper texture name: %s",
                    render_assets.patches[side_info->texture_info_upper.texture_id].name.c_str());
                if (ImGui::InputScalar("upper texture_id", ImGuiDataType_U16,
                                       (void*)(&side_info->texture_info_upper.texture_id),
                                       (void*)(&step_u16), (void*)(NULL), "%d", flags)) {
                    // Ensure it is in bounds.
                    usize n_patches = render_assets.patches.size();
                    if (side_info->texture_info_upper.texture_id >= n_patches) {
                        side_info->texture_info_upper.texture_id = n_patches - 1;
                    }
                }

                ImGui::InputScalar("upper x_offset", ImGuiDataType_S16,
                                   (void*)(&side_info->texture_info_upper.x_offset),
                                   (void*)(&step_i16), (void*)(NULL), "%d", flags);
                ImGui::InputScalar("upper y_offset", ImGuiDataType_S16,
                                   (void*)(&side_info->texture_info_upper.y_offset),
                                   (void*)(&step_i16), (void*)(NULL), "%d", flags);

                ImGui::Separator();
                const doom::Patch& middle_patch =
                    render_assets.patches[side_info->texture_info_middle.texture_id];
                ImGui::Text("middle patch name: %s", middle_patch.name.c_str());
                ImGui::Text("middle patch size x: %d", middle_patch.size_x);
                ImGui::Text("middle patch size y: %d", middle_patch.size_y);
                if (ImGui::InputScalar("middle texture_id", ImGuiDataType_U16,
                                       (void*)(&side_info->texture_info_middle.texture_id),
                                       (void*)(&step_u16), (void*)(NULL), "%d", flags)) {
                    // Ensure it is in bounds.
                    usize n_patches = render_assets.patches.size();
                    if (side_info->texture_info_middle.texture_id >= n_patches) {
                        side_info->texture_info_middle.texture_id = n_patches - 1;
                    }
                }
                ImGui::InputScalar("middle x_offset", ImGuiDataType_S16,
                                   (void*)(&side_info->texture_info_middle.x_offset),
                                   (void*)(&step_i16), (void*)(NULL), "%d", flags);
                ImGui::InputScalar("middle y_offset", ImGuiDataType_S16,
                                   (void*)(&side_info->texture_info_middle.y_offset),
                                   (void*)(&step_i16), (void*)(NULL), "%d", flags);

                ImGui::Separator();
                ImGui::Text(
                    "lower texture name: %s",
                    render_assets.patches[side_info->texture_info_lower.texture_id].name.c_str());
                if (ImGui::InputScalar("lower texture_id", ImGuiDataType_U16,
                                       (void*)(&side_info->texture_info_lower.texture_id),
                                       (void*)(&step_u16), (void*)(NULL), "%d", flags)) {
                    // Ensure it is in bounds.
                    usize n_patches = render_assets.patches.size();
                    if (side_info->texture_info_lower.texture_id >= n_patches) {
                        side_info->texture_info_lower.texture_id = n_patches - 1;
                    }
                }
                ImGui::InputScalar("lower x_offset", ImGuiDataType_S16,
                                   (void*)(&side_info->texture_info_lower.x_offset),
                                   (void*)(&step_i16), (void*)(NULL), "%d", flags);
                ImGui::InputScalar("lower y_offset", ImGuiDataType_S16,
                                   (void*)(&side_info->texture_info_lower.y_offset),
                                   (void*)(&step_i16), (void*)(NULL), "%d", flags);

                ImGui::Separator();

                if (ImGui::InputScalar("portal qe index", ImGuiDataType_U64,
                                       (void*)(&side_info->qe_portal.i), (void*)(NULL),
                                       (void*)(NULL), "%lu",
                                       ImGuiInputTextFlags_EnterReturnsTrue)) {
                    // Ensure the other quarter edge exists, that both are passable, and that both
                    // are the same length.
                    core::SideInfo* side_info_portal =
                        map.GetEditableSideInfo(side_info->qe_portal);
                    if (side_info->qe == side_info->qe_portal || side_info_portal == nullptr) {
                        // Invalid!
                        side_info->qe_portal = {core::kInvalidIndex};
                        std::cout << "Both are equal (" << (side_info->qe == side_info->qe_portal)
                                  << ") or side_info_portal is null ("
                                  << (side_info_portal == nullptr) << ")" << std::endl;
                    } else {
                        // Check for same length
                        auto sqdist_ab = common::SqNorm(mesh.GetVertex(side_info->qe) -
                                                        mesh.GetVertex(mesh.Sym(side_info->qe)));
                        auto sqdist_cd =
                            common::SqNorm(mesh.GetVertex(side_info->qe_portal) -
                                           mesh.GetVertex(mesh.Sym(side_info->qe_portal)));
                        if (abs(sqdist_ab - sqdist_cd) > 0.1) {
                            std::cout << "Not of the same length! " << sqdist_ab << " vs. "
                                      << sqdist_cd << std::endl;
                            side_info->qe_portal = {core::kInvalidIndex};
                        } else {
                            // Ensure both are passable. Note that the far side may or may not be
                            // passable.
                            side_info->flags |= core::kSideInfoFlag_PASSABLE;
                            side_info_portal->flags |= core::kSideInfoFlag_PASSABLE;

                            // Ensure the other one is set to the same value.
                            side_info_portal->qe_portal = side_info->qe;

                            std::cout << "Connected portals!" << std::endl;
                        }
                    }
                }

                ImGui::Separator();

                bool updated_sector_id = false;
                if (ImGui::InputScalar("sector index", ImGuiDataType_U16,
                                       (void*)(&side_info->sector_id), (void*)(&step_u16),
                                       (void*)(NULL), "%d", flags)) {
                    updated_sector_id = true;
                    // Ensure that it lies in bounds
                    if (side_info->sector_id > map.GetMaxSectorIndex()) {
                        side_info->sector_id = map.GetMaxSectorIndex();
                    }
                }
                if (ImGui::Button("New sector index")) {
                    u16 sector_index = map.AddSector();
                    side_info->sector_id = sector_index;
                    updated_sector_id = true;
                }
                if (updated_sector_id) {
                    // Walk around the triangle and ensure all edges have the same sector id.
                    core::QuarterEdgeIndex qe = mesh.Prev(mesh.Sym(selected_edge_index));
                    while (qe != selected_edge_index) {
                        core::SideInfo* side_info2 = map.GetEditableSideInfo(qe);
                        if (side_info2 != nullptr) {
                            side_info2->sector_id = side_info->sector_id;
                        }

                        qe = mesh.Prev(mesh.Sym(qe));
                    }
                }

                core::Sector* sector = map.GetEditableSector(side_info->sector_id);
                if (sector == nullptr) {
                    ImGui::Text("SECTOR INDEX IS INVALID");
                } else {
                    ImGui::Text("sector flags:      %X", sector->flags);
                    if (ImGui::InputScalar("z_ceil", ImGuiDataType_Float, (void*)(&sector->z_ceil),
                                           (void*)(&step_f32), (void*)(NULL), "%.3f", flags)) {
                        sector->z_ceil = std::max(sector->z_ceil, sector->z_floor);
                    }
                    if (ImGui::InputScalar("z_floor", ImGuiDataType_Float, &sector->z_floor,
                                           (void*)(&step_f32), (void*)(NULL), "%.3f", flags)) {
                        sector->z_floor = std::min(sector->z_floor, sector->z_ceil);
                    }
                    if (ImGui::InputScalar("ceiling flat_id", ImGuiDataType_U16,
                                           (void*)(&sector->flat_id_ceil), (void*)(&step_u16),
                                           (void*)(NULL), "%d", flags)) {
                        // Ensure it is in bounds.
                        usize n_flats = render_assets.flats.size();
                        if (sector->flat_id_ceil >= n_flats) {
                            sector->flat_id_ceil = n_flats - 1;
                        }
                    }
                    if (ImGui::InputScalar("floor flat_id", ImGuiDataType_U16,
                                           (void*)(&sector->flat_id_floor), (void*)(&step_u16),
                                           (void*)(NULL), "%d", flags)) {
                        // Ensure it is in bounds.
                        usize n_flats = render_assets.flats.size();
                        if (sector->flat_id_floor >= n_flats) {
                            sector->flat_id_floor = n_flats - 1;
                        }
                    }
                }

            } else {
                ImGui::Text("there is no sideinfo for this edge");
                if (ImGui::Button("Create SideInfo")) {
                    bool success = map.AddSideInfo(selected_edge_index);
                    if (!success) {
                        std::cout << "Failed to add a side info" << std::endl;
                    }
                }
            }
            ImGui::End();
        }

        // ImGUI Rendering
        ImGui::Render();
        SDL_RenderSetScale(editor_window_data.renderer, io.DisplayFramebufferScale.x,
                           io.DisplayFramebufferScale.y);
        ImGui_ImplSDLRenderer2_RenderDrawData(ImGui::GetDrawData());

        // ------------------------------------------------------------------------------------------------
        {
            // Update the camera
            f32 dt = 1.0f / 30.0f;  // TODO
            MoveCamera(&player_cam, keyboard_state, map, dt);

            // Render the player view.
            RenderScene(player_view_pixels, raycast_distances, dw, active_spans, map, player_cam,
                        render_assets, render_data);

            SDL_UpdateTexture(player_window_data.texture, NULL, player_view_pixels,
                              player_window_data.screen_size_x * 4);
            SDL_RenderCopyEx(player_window_data.renderer, player_window_data.texture, NULL, NULL,
                             0.0, NULL, SDL_FLIP_VERTICAL);
        }

        // Decay the keyboard state
        DecayKeyboardState(&keyboard_state);

        // SDL_RENDERER_PRESENTVSYNC means this is syncronized with the monitor
        // refresh rate. (30Hz)
        SDL_RenderPresent(editor_window_data.renderer);
        SDL_RenderPresent(player_window_data.renderer);

        EndAndPrintProfile(cpu_freq);

        // Maintain a running average of the total tick time.
        static f64 total_cpu_elapsed_ms_smoothed = 60.000;  // Init close to current reading
        constexpr f64 kTotalCpuElapsedSmoothingConst = 0.95;
        u64 total_cpu_elapsed = kGlobalProfiler.tsc_end - kGlobalProfiler.tsc_start;
        f64 total_cpu_elapsed_ms = 1000.0 * (f64)total_cpu_elapsed / (f64)cpu_freq;
        total_cpu_elapsed_ms_smoothed =
            total_cpu_elapsed_ms_smoothed * kTotalCpuElapsedSmoothingConst +
            total_cpu_elapsed_ms * (1.0 - kTotalCpuElapsedSmoothingConst);
        printf("Total time smoothed: %0.4fms\n", total_cpu_elapsed_ms_smoothed);

        ResetProfiler();
    }
}