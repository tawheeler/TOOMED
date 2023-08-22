#pragma once

#include <string.h>

#include <cstdio>

#include "platform_metrics.cpp"
#include "typedefs.hpp"

#define ArrayCount(Array) (sizeof(Array) / sizeof((Array)[0]))

#ifndef PROFILER
#define PROFILER 1
#endif

#if PROFILER

struct ProfileAnchor {
    u64 tsc_elapsed_exclusive;  // Does _not_ include children
    u64 tsc_elapsed_inclusive;  // Does include children
    u64 hit_count;
    char const* label;
};

static ProfileAnchor kGlobalProfilerAnchors[4096];
static u32 kGlobalProfilerParent;

struct ProfileBlock {
    ProfileBlock(char const* label_, u32 anchor_index_) {
        parent_index = kGlobalProfilerParent;
        anchor_index = anchor_index_;
        label = label_;

        ProfileAnchor* anchor = kGlobalProfilerAnchors + anchor_index_;
        old_tsc_elapsed_inclusive = anchor->tsc_elapsed_inclusive;

        kGlobalProfilerParent = anchor_index;
        tsc_start = ReadCPUTimer();
    }

    ~ProfileBlock(void) {
        u64 elapsed = ReadCPUTimer() - tsc_start;
        kGlobalProfilerParent = parent_index;

        ProfileAnchor* parent = kGlobalProfilerAnchors + parent_index;
        ProfileAnchor* anchor = kGlobalProfilerAnchors + anchor_index;

        parent->tsc_elapsed_exclusive -= elapsed;
        anchor->tsc_elapsed_exclusive += elapsed;
        anchor->tsc_elapsed_inclusive = old_tsc_elapsed_inclusive + elapsed;
        anchor->hit_count += 1;

        /* NOTE(casey): This write happens every time solely because there is no
            straightforward way in C++ to have the same ease-of-use. In a better programming
            language, it would be simple to have the anchor points gathered and labeled at compile
            time, and this repetative write would be eliminated. */
        anchor->label = label;
    }

    char const* label;
    u64 tsc_start;
    u64 old_tsc_elapsed_inclusive;
    u32 parent_index;
    u32 anchor_index;
};

#define NameConcat2(A, B) A##B
#define NameConcat(A, B) NameConcat2(A, B)
#define TimeBlock(Name) ProfileBlock NameConcat(Block, __LINE__)(Name, __COUNTER__ + 1);
#define ProfilerEndOfCompilationUnit                                \
    static_assert(__COUNTER__ < ArrayCount(kGlobalProfilerAnchors), \
                  "Number of profile points exceeds size of ProfileAnchors array")

static void PrintTimeElapsed(u64 total_tsc_elapsed, ProfileAnchor* anchor) {
    f64 percent = 100.0 * ((f64)anchor->tsc_elapsed_exclusive / (f64)total_tsc_elapsed);
    printf("  %-20s[%6lu]: %10lu (%.2f%%", anchor->label, anchor->hit_count,
           anchor->tsc_elapsed_exclusive, percent);
    if (anchor->tsc_elapsed_inclusive != anchor->tsc_elapsed_exclusive) {
        f64 percent_with_children =
            100.0 * ((f64)anchor->tsc_elapsed_inclusive / (f64)total_tsc_elapsed);
        printf(", %.2f%% w/children", percent_with_children);
    }
    printf(")\n");
}

static void PrintAnchorData(u64 total_cpu_elapsed) {
    for (u32 i_anchor = 0; i_anchor < ArrayCount(kGlobalProfilerAnchors); ++i_anchor) {
        ProfileAnchor* anchor = kGlobalProfilerAnchors + i_anchor;
        if (anchor->tsc_elapsed_inclusive) {
            PrintTimeElapsed(total_cpu_elapsed, anchor);
        }
    }
}

#else

#define TimeBlock(...)
#define PrintAnchorData(...)
#define ProfilerEndOfCompilationUnit

#endif

struct Profiler {
    u64 tsc_start;
    u64 tsc_end;
};
static Profiler kGlobalProfiler;

#define TimeFunction TimeBlock(__func__)

static void BeginProfile(void) { kGlobalProfiler.tsc_start = ReadCPUTimer(); }

// u64 cpu_freq = EstimateCPUTimerFreq(100);
static void EndAndPrintProfile(u64 cpu_freq) {
    kGlobalProfiler.tsc_end = ReadCPUTimer();

    u64 total_cpu_elapsed = kGlobalProfiler.tsc_end - kGlobalProfiler.tsc_start;

    if (cpu_freq) {
        printf("\nTotal time: %0.4fms (CPU freq %lu)\n",
               1000.0 * (f64)total_cpu_elapsed / (f64)cpu_freq, cpu_freq);
    }

    PrintAnchorData(total_cpu_elapsed);
}

static void ResetProfiler(void) {
    memset(kGlobalProfilerAnchors, 0, sizeof(kGlobalProfilerAnchors));
}

// Starting time:
// Interestingly, RenderPatchColumn takes almost no time. Its all in RenderWallsInner.

// Total time: 261.2309ms (CPU freq 1992020900)
//   MoveCamera          [     1]:       6751 (0.00%)
//   RenderGrid          [     1]:       9434 (0.00%)
//   RenderPatchColumn   [  1141]:    5681896 (1.09%)
//   RenderWallsInner    [ 12433]:  498245124 (95.75%, 96.84% w/children)
//   RenderWalls         [     1]:     218863 (0.04%, 96.88% w/children)
//   SDLPollEvents       [     1]:      65963 (0.01%)

// Total time: 282.9411ms (CPU freq 1991977930)
//   GetRightHandedness  [ 48790]:    5595253 (0.99%)
//   RenderPatchColumn   [  1082]:    5637169 (1.00%)
//   RenderWallsInner    [ 12213]:  533552651 (94.67%, 97.74% w/children)
//   RenderWallsInner_WithValidSideInfo[  3344]:    6089453 (1.08%, 83.27% w/children)
//   RenderWalls         [     1]:     241069 (0.04%, 97.78% w/children)

// After making raycasting a while loop rather than recursive
// Total time: 33.0318ms (CPU freq 1991988590)
//   GetRightHandedness  [ 32890]:    3095990 (4.71%)
//   RenderPatchColumn   [  1293]:    5331569 (8.10%)
//   RenderWallsInner    [   640]:   34651635 (52.66%, 74.61% w/children)
//   RenderWallsInner_EdgeChecks[  8222]:    2856858 (4.34%, 6.16% w/children)
//   RenderWallsInner_WithValidSideInfo[  2257]:    3172907 (4.82%, 12.92% w/children)
//   RenderWalls         [     1]:     129132 (0.20%, 74.81% w/children)

// After adding floor and ceiling rendering
// Total time: 68.4885ms (CPU freq 1991996640)
//   GetRightHandedness  [ 38878]:    2110694 (1.55%)
//   RenderPatchColumn   [  1623]:   14123995 (10.35%)
//   RenderSpan          [   502]:    6697976 (4.91%)
//   Raycast             [   640]:   70439031 (51.63%, 76.19% w/children)
//   Raycast_EdgeChecks  [ 12959]:    6179480 (4.53%, 6.06% w/children)
//   Raycast_WithValidSideInfo[  3390]:    7078638 (5.19%, 18.50% w/children)
//   RenderScene         [     1]:     242694 (0.18%, 78.32% w/children)
// Total time smoothed: 69.2931ms

// After skipping one triangle check
// Total time: 66.7742ms (CPU freq 1992002070)
//   GetRightHandedness  [ 25919]:    1360790 (1.02%)
//   RenderPatchColumn   [  1623]:   13847854 (10.41%)
//   RenderSpan          [   502]:    7105657 (5.34%)
//   Raycast             [   640]:   72094385 (54.20%, 78.40% w/children)
//   Raycast_EdgeChecks  [ 12959]:    5390282 (4.05%, 5.07% w/children)
//   Raycast_WithValidSideInfo[  3390]:    7211444 (5.42%, 19.13% w/children)
//   RenderScene         [     1]:     319789 (0.24%, 80.68% w/children)
// Total time smoothed: 67.0397ms