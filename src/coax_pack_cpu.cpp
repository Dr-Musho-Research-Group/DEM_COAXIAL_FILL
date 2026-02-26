/*
 * coax_pack_cpu.cpp — 3D hard-sphere packing in a coaxial (annular) column
 *
 * Similar style/structure to md_cpu.cpp but:
 *  - 3D domain (x,y,z) with inner/outer cylindrical walls and finite length L
 *  - gravity-driven fill from the top; diameter distribution sampling
 *  - hard-sphere collisions with simple tangential damping
 *  - optional vertical shaking (container acceleration) during fill
 *  - optional top ram/piston compression at end
 *
 * CLI:
 *   coax_pack_cpu
 *     natoms_max         (int)   max spheres we will inject total
 *     dt                  time step [s]
 *     niter              (int)   number of steps
 *     dump_interval      (int)   dump every N steps (XYZ/VTK)
 *     debug              (int)   0/1
 *     seed               (int)   RNG seed (<=0 -> time-based)
 *     Rin                 inner radius  [m]
 *     Rout                outer radius  [m]
 *     L                   length (z)    [m]
 *     flux               (int)   spheres per second injected while filling
 *     g                   gravity magnitude [m/s^2] (points -z)
 *     shake_hz            shake frequency [Hz]; 0 => off
 *     shake_amp           shake amplitude [m] (vertical)
 *     fill_time           seconds to inject (afterwards stop injecting)
 *     ram_start_t         when to start ram (s from t=0)
 *     ram_duration        duration of ram [s]
 *     ram_speed           downward piston speed [m/s] (top plate)
 *
 * Output:
 *   atoms_<iter>.xyz, atom.<iter>.vtk (positions/velocities)
 *
 * Notes:
 *   - Neighbor search via uniform 3D grid (cell-linked list),
 *     cell size ~ max_diameter → ~27-cell checks.
 *   - Contact: impulse-like normal restitution + tangential damping.
 *   - Walls: inner cylinder (r >= Rin+ri), outer cylinder (r <= Rout-ri),
 *            bottom plane (z >= ri), top piston plane (z <= z_top(t)-ri).
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>
#include <random>
#include <string>
#include <array>
#include <limits>
#include <cstdint>

#include <unordered_map>
#include <cctype>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// --------------------- Threading (OpenMP preferred) ---------------------
// This code targets MSVC/Windows and GCC/Clang/Linux. Prefer OpenMP for
// per-particle parallel loops. If OpenMP is not enabled, execution falls back
// to serial loops (still functional, just slower).

#if defined(_OPENMP)
  #include <omp.h>
#endif

// 0 => OpenMP default, 1 => serial, >1 => request that many threads (OpenMP)
static int g_threads = 0;

static inline int effective_threads(){
#if defined(_OPENMP)
    if (g_threads > 0) return g_threads;
    return omp_get_max_threads();
#else
    return 1;
#endif
}

static inline void apply_thread_request(){
#if defined(_OPENMP)
    if (g_threads > 0) omp_set_num_threads(g_threads);
#endif
}


// ------------------------------ Types ------------------------------
using real = double;

struct Vec3 {
    real x,y,z;
};
struct Particle {
    Vec3 p;    // position
    Vec3 v;    // velocity
    Vec3 a;    // acceleration (reset each step)
    real r;   // radius
    real m;   // mass
    uint8_t type_id = 0;
    // --- new: sleep state ---
    bool asleep = false;      // if true: skip pair & wall collisions
    int  calm_steps = 0;      // consecutive steps below thresholds
    Vec3 prev_p = {0,0,0};    // last position (for drift check)
};

struct Monitor {
    real kinE;
    int   n;
    real avgZ;
};

static inline Vec3 operator+(const Vec3&a,const Vec3&b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
static inline Vec3 operator-(const Vec3&a,const Vec3&b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
static inline Vec3 operator*(const Vec3&a,real s){ return {a.x*s,a.y*s,a.z*s}; }
static inline real dot(const Vec3&a,const Vec3&b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline real norm(const Vec3&a){ return std::sqrt(dot(a,a)); }

// ------------------------- User/Global Params ----------------------
int   g_debug = 0;
int   g_seed  = 0;
// packing cutoff by volume fraction (disabled if <= 0)
real g_phi_target = 0.0;
// live gravity used in integration (we can flip to 0 when phi_target is reached)
real g_g_run = 0.0;
bool g_phi_reached = false;

// ---------------- Early-stop after phi target ----------------
// If phi_target is reached we already disable injection and switch to zero-g.
// These knobs optionally stop the run once the ensemble is "settled".
// Defaults keep old behavior (run until niter).
real g_stop_vrms = 0.0;            // [m/s] stop if RMS speed below this (0 disables)
real g_stop_vmax = 0.0;            // [m/s] stop if max speed below this (0 disables)
real g_stop_sleep_frac = 0.0;      // [0..1] stop if >= this fraction asleep (0 disables)
int  g_stop_check_interval = 200;  // [steps] how often to evaluate stop criteria
int  g_stop_checks_required = 10;  // number of consecutive successful checks

int   g_natoms_max;
real g_dt;
int   g_niter;
int   g_dump_interval;
int   g_xyz_interval = -1;       // -1 => follow dump_interval
int   g_vtk_interval = -1;       // -1 => follow dump_interval
int   g_vtk_domain_interval = 0; // 0 => off
int   g_vtk_domain_segments = 96; // cylinder tessellation

real g_Rin, g_Rout, g_L;
int   g_flux;                  // spheres/sec
real g_g = 9.81f;             // m/s^2 downward

// shaking (container vertical acceleration via frame motion)
real g_shake_hz = 1000;
real g_shake_amp = 9.81;      // meters

// shaking (container acceleration). Existing flags:
//   --shake_hz, --shake_amp  (legacy: used as vertical amplitude)
// New (optional) flags allow independent horizontal/vertical shaking:
//   --shake_amp_x, --shake_amp_y, --shake_amp_z
// Compatibility:
//   --shake_xy_legacy 1 keeps the historical behavior where the vertical
//   shake acceleration was also applied to x and y. Set to 0 to use the
//   independent axis amplitudes.
real g_shake_amp_x = 0.0;  // meters
real g_shake_amp_y = 0.0;  // meters
real g_shake_amp_z = 0.0;  // meters (alias of g_shake_amp)
int  g_shake_xy_legacy = 1;

// contact geometry cushion: expands the effective collision radius by this amount
// (useful for adding a "soft" clearance without changing the physical particle size)
real g_cushion = 0.0; // meters

// optional wall spring (penalty) model: adds a spring-damper normal response
// on top of the existing impulse + positional correction wall handling.
// Set g_wall_k <= 0 to disable.
real g_wall_k = 0.0;        // N/m
real g_wall_zeta = 0.20;    // damping ratio (0..1)
real g_wall_dvmax = 5.0;    // cap on per-step dv from spring-damper [m/s]

// ---------------- Repulsive (non-contact) forces ----------------
// Optional short-range repulsion to help particles "spread" during filling.
// This is additive and does not replace the existing collision handling.
//
// Enable by setting k > 0 and a positive range. Defaults keep old behavior.
real g_repulse_range = 0.0;   // [m] active range in *gap* space (0 disables)
real g_repulse_k_pp  = 0.0;   // [N/m] particle-particle repulsion strength
real g_repulse_k_pw  = 0.0;   // [N/m] particle-wall repulsion strength

int  g_repulse_use_mass = 0; // 0: k is accel scale [m/s^2]; 1: k behaves like spring [N/m] divided by mass
real g_repulse_dvmax = 0.0;  // [m/s] per-step velocity clamp from repulsion (0 disables)
// injection initial velocity (set via CLI)
real g_inject_vx = 0.0;
real g_inject_vy = 0.0;
real g_inject_vz = -0.5;

// piston (top wall) normal velocity component along the contact normal
// (updated once per step in the main loop; used for wall spring damping)
real g_piston_wall_vn = 0.0;

// drag
real g_air_visc = 1.8e-5;      // Pa·s (air)
real g_drag_mult = 1.0;        // 0 to disable, 1 for physical-ish

// injection
real g_fill_time = 0.0;      // seconds to inject

// ram (top piston)
real g_ram_t0 = 0.0;
real g_ram_dt = 0.0;
real g_ram_speed = 0.0;      // downward (positive)

// materials / contacts
real g_density = 2000.0;     // kg/m^3 (default)
real g_e_pp = 0.45f;          // restitution particle-particle
real g_e_pw = 0.45f;          // restitution particle-wall
real g_tangent_damp = 0.85f;  // simple tangential vel damping

// global linear velocity damping (optional)
// Applies v <- v * exp(-lin_damp * dt) each step (0 disables).
real g_lin_damp = 0.0; // [1/s]

// neighbor grid
real g_cell_h;                // cell size ~ max_diam
int   gx, gy, gz;
Vec3  g_minB, g_maxB;

// Make sleeping less eager and waking easier
real g_sleep_v   = 5e-6;      // was 2e-3; require smaller speeds before sleeping
real g_sleep_dz  = 2e-8;      // was 5e-7; allow a bit more micro drift per step
int  g_sleep_N   = 5000;        // was 25; sustain calm longer before sleeping
real g_wake_vn   = 5e-4;      // was 2e-4; even gentle taps wake sleepers
real g_wake_band = 2.5;       // was 6.0; only wake near the actual piston zone

// successive damping (multi-sweep velocity decay)
int   g_sd_sweeps     = 2;       // was 4
real  g_sd_alpha0     = 0.03;    // was 0.05
real  g_sd_alpha1     = 0.10;    // was 0.15
real  g_sd_vcut       = 1e-6;    // was 5e-8; avoid freezing “whisper” motion
real  g_sd_top_guardD = 0.25;     // damp less near the lid

// -- wall contact knobs --
real g_e_bottom       = 0.20;    // a bit bouncier so energy isn’t all in-plane
real g_floor_mu_t     = 0.55;    // was 0.75; less skate, but not glue
real g_floor_roll     = 0.20;    // was 0.30
real g_floor_beta     = 0.25;    // keep mild Baumgarte damping on penetration
real g_vstick_floor   = 0.010;   // was 0.005; low-speed impacts stick sooner
real g_vstatic_floor  = 0.020;   // was 5 (!) m/s; snap-to-rest only at tiny slip speeds
real g_sd_tan_boost   = 0.10;    // smaller extra XY damping near floor
real g_wall_mu_t   = 0.40;  // tangential friction for cylindrical walls
real g_wall_roll   = 0.10;  // rolling loss on walls
real g_wall_beta   = 0.50;  // Baumgarte fraction on walls (0.0..0.5 reasonable)

// Optional: geometric wall roughness on the cylindrical walls.
// Implemented as a deterministic sinusoidal variation of the effective inner/outer
// wall radii vs (theta,z). This provides wall "texture" that can reduce strong
// near-wall layering/segregation for polydisperse packs.
// Set g_wall_rough_amp <= 0 to disable.
real g_wall_rough_amp   = 0.0;   // [m] peak roughness amplitude
int  g_wall_rough_mth   = 8;     // azimuthal modes (>=1)
int  g_wall_rough_mz    = 3;     // axial modes (>=0)
real g_wall_rough_phth  = 0.0;   // [rad] azimuthal phase (set from seed)
real g_wall_rough_phz   = 0.0;   // [rad] axial phase (set from seed)
// -- particle-particle contact knobs --
real g_pair_mu_s     = 0.60;   // static friction (stick if within cap)
real g_pair_mu_k     = 0.40;   // kinetic friction  (<= mu_s)
real g_pair_vstatic  = 0.02;   // m/s; snap tangential slip below this
real g_pair_beta     = 0.20;   // Baumgarte fraction for positional stabilization
real g_pair_sep_cap  = 0.25;   // cap separation to this frac of (A.r + B.r)


// Each particle type gets its own (diam_um, cumu)
struct SizeDist {
    std::vector<real> diam_um;
    std::vector<real> cumu;     // monotonically nondecreasing, last = 1
};

// Four types: 0 (baseline) + three additional
static std::array<SizeDist,4> g_types = {};

// Composition (fractions); will be normalized at runtime
//                             La       Sr      Fe      Co
static real g_type_frac[4] = {0.4667, 0.2020, 0.2945, 0.0368};

// Optional: per-type density (defaults to global density)
static real g_density_type[4] = {0,0,0,0}; // 0 => use g_density

// Initialize defaults for type 0 to match legacy arrays
static void init_default_distributions_once() {
    static bool inited = false;
    if (inited) return;
    inited = true;
    //La2O3
    g_types[0].diam_um = {
        1.279151515,1.67730303,2.075454545,2.473606061,
        2.871757576,3.269909091,3.668060606,4.066212121,
        4.464363636,4.862515152,5.260666667,5.658818182,
        6.056969697,6.455121212,6.853272727,7.251424242,
        7.649575758,8.047727273,8.445878788,8.844030303,
        9.242181818,9.640333333,10.03848485,10.43663636,
        10.83478788,11.23293939,11.63109091,12.02924242,
        12.42739394,12.82554545,13.22369697,13.62184848,14.02
    };
    g_types[0].cumu = {
        0.081371888,0.115647757,0.159102372,0.212064199,
        0.27411818,0.344014775,0.419701917,0.49849134,
        0.577339844,0.653197383,0.723356443,0.785736889,
        0.839057226,0.882871634,0.917483148,0.943767896,
        0.962957522,0.976425694,0.985512895,0.991407176,
        0.995082628,0.997285907,0.998555626,0.999259062,
        0.999633706,0.999825526,0.999919942,0.999964619,
        0.999984942,0.999993829,0.999997565,0.999999075,0.999999662
    };
    //SrO
    g_types[1].diam_um = {
        3.876095745,7.752191489,11.62828723,15.50438298,
        19.38047872,23.25657447,27.13267021,31.00876596,
        34.8848617,38.76095745,42.63705319,46.51314894,
        50.38924468,54.26534043,58.14143617,62.01753191,
        65.89362766,69.7697234,73.64581915,77.52191489,
        81.39801064,85.27410638,89.15020213,93.02629787,
        96.90239362,100.7784894,104.6545851,108.5306809,
        112.4067766,116.2828723,120.1589681,124.0350638,
        127.9111596,131.7872553,135.6633511,139.5394468,
        143.4155426,147.2916383,151.167734,155.0438298,
        158.9199255,162.7960213,166.672117
    };
    g_types[1].cumu = {
        0.162790218,0.199961103,0.241818325,0.288022579,
        0.338018849,0.391050765,0.446192587,0.502396813,
        0.558553427,0.613555231,0.666362905,0.716063457,
        0.761916637,0.803385467,0.840149079,0.872098179,
        0.89931536,0.922043909,0.940649538,0.955579579,
        0.967323733,0.976379526,0.983224554,0.988296423,
        0.991980298,0.994603221,0.996433893,0.997686397,
        0.998526421,0.999078687,0.999434605,0.999659455,
        0.999798701,0.999883232,0.999933535,0.999962878,
        0.999979657,0.999989063,0.999994231,0.999997015,
        0.999998485,0.999999245,0.999999631
    };
    //Fe2O3
    g_types[2].diam_um = {
        0.552982759,1.105965517,1.658948276,2.211931034,
        2.764913793,3.317896552,3.87087931,4.423862069,
        4.976844828,5.529827586,6.082810345,6.635793103,
        7.188775862,7.741758621,8.294741379,8.847724138,
        9.400706897,9.953689655,10.50667241,11.05965517,
        11.61263793,12.16562069,12.71860345,13.27158621,
        13.82456897,14.37755172,14.93053448,15.48351724,
        16.0365,16.58948276,17.14246552,17.69544828,
        18.24843103,18.80141379,19.35439655,19.90737931,
        20.46036207,21.01334483,21.56632759,22.11931034,
        22.6722931,23.22527586,23.77825862,24.33124138,
        24.88422414,25.4372069,25.99018966,26.54317241,
        27.09615517,27.64913793,28.20212069,28.75510345,
        29.30808621,29.86106897,30.41405172,30.96703448,
        31.52001724,32.073
    };
    g_types[2].cumu = {
        0.086571103,0.107758337,0.132446772,0.160779115,
        0.192800567,0.228443149,0.267515262,0.309697952,
        0.35454886,0.401514239,0.449948651,0.49914126,
        0.548346946,0.596819995,0.643847832,0.688782267,
        0.731065986,0.770252505,0.806018459,0.838167834,
        0.866628488,0.89144192,0.912747747,0.930764625,
        0.945769429,0.958076386,0.968017627,0.975926238,
        0.982122494,0.986903597,0.990536856,0.993256013,
        0.995260216,0.996715068,0.997755145,0.998487433,
        0.998995203,0.999341958,0.999575169,0.999729638,
        0.999830402,0.999895137,0.999936095,0.999961617,
        0.99997728,0.999986746,0.99999238,0.999995683,
        0.99999759,0.999998674,0.999999281,0.999999616,
        0.999999798,0.999999895,0.999999946,0.999999973,
        0.999999987,0.999999993
    };
    //Co3O4
    g_types[3].diam_um = {
        1.7381,2.1292,2.5203,2.9114,3.3025,3.6936,
        4.0847,4.4758,4.8669,5.258,5.6491,6.0402,
        6.4313,6.8224,7.2135,7.6046,7.9957,8.3868,
        8.7779,9.169,9.5601,9.9512,10.3423,10.7334,
        11.1245,11.5156,11.9067,12.2978,12.6889,13.08
    };
    g_types[3].cumu = {
        0.055120161,0.082628135,0.119377722,0.166372915,
        0.223898575,0.291301232,0.36689739,0.448055096,
        0.531455307,0.613492709,0.690736575,0.76035488,
        0.820415552,0.870013556,0.909218912,0.938883188,
        0.960367879,0.975262569,0.985146749,0.991425257,
        0.995242757,0.997464571,0.99870235,0.999362411,
        0.999699335,0.999863955,0.999940947,0.999975414,
        0.999990184,0.999996242
    };
    g_density_type[0] = 2650.0;  // kg/m^3
    g_density_type[1] = 3700.0;
    g_density_type[2] = 5250.0;
    g_density_type[3] = 4130.0;
}

// --------------------------- Utilities -----------------------------
static inline int clampi(int v,int lo,int hi){ return v<lo?lo:(v>hi?hi:v); }

// NEW: portable clamp for arithmetic types (no <algorithm> / C++17 needed)
template <typename T>
static inline T clamp_val(T v, T lo, T hi) { return v < lo ? lo : (v > hi ? hi : v); }


// helper map for element symbols by type_id
static const char* element_symbol(uint8_t tid) {
    switch (tid) {
        case 0: return "La";
        case 1: return "Sr";
        case 2: return "Fe";
        case 3: return "Co";
        default: return "X";
    }
}

void dump_xyz(const std::vector<Particle>& P, int iter){
    char fn[64]; std::snprintf(fn,sizeof(fn),"atoms_%d.xyz",iter);
    FILE* f = std::fopen(fn,"w");
    if(!f){ std::perror("xyz"); return; }
    std::fprintf(f,"%zu\n", P.size());
    std::fprintf(f,"iter %d\n", iter);
    for(const auto& p: P)
        std::fprintf(f,"%s %.9g %.9g %.9g %.9g\n",
                     element_symbol(p.type_id),
                     p.p.x, p.p.y, p.p.z, 2*p.r);
    std::fclose(f);
}


// VTK domain geometry (annulus walls + bottom + current top plane).
// Outputs a simple POLYDATA surface mesh (triangulated quads) for visualization.
// NOTE: This is visualization-only. It does not represent a CFD mesh.
void dump_domain_vtk(real top_z, int iter){
    if (g_vtk_domain_interval <= 0) return;
    char fn[96]; std::snprintf(fn,sizeof(fn),"domain.%d.vtk",iter);
    FILE* f = std::fopen(fn,"w");
    if(!f){ std::perror("vtk-domain"); return; }

    const int N = std::max(12, g_vtk_domain_segments);
    const real z0 = 0.0;
    const real z1 = top_z; // current piston plane
    const real rin = g_Rin;
    const real rout = g_Rout;

    // points: inner ring (z0,z1) + outer ring (z0,z1)
    const int nPts = 4*N;
    std::fprintf(f,"# vtk DataFile Version 2.0\nCoax domain\nASCII\nDATASET POLYDATA\n");
    std::fprintf(f,"POINTS %d float\n", nPts);

    auto emit_ring = [&](real r, real z){
        for(int i=0;i<N;++i){
            real th = (2.0*M_PI*real(i))/real(N);
            real x = r*std::cos(th);
            real y = r*std::sin(th);
            std::fprintf(f,"%g %g %g\n", x, y, z);
        }
    };

    // order: inner(z0), inner(z1), outer(z0), outer(z1)
    emit_ring(rin, z0);
    emit_ring(rin, z1);
    emit_ring(rout, z0);
    emit_ring(rout, z1);

    auto idxPt = [&](int ring, int i)->int{
        i = (i%N + N)%N;
        return ring*N + i;
    };

    const int nTris = 8*N; // 2 tris per quad, 4 quad sets
    std::fprintf(f,"POLYGONS %d %d\n", nTris, nTris*4);

    auto emit_quad_as_tris = [&](int a,int b,int c,int d){
        std::fprintf(f,"3 %d %d %d\n", a,b,c);
        std::fprintf(f,"3 %d %d %d\n", a,c,d);
    };

    for(int i=0;i<N;++i){
        int i1 = (i+1)%N;

        // inner wall
        emit_quad_as_tris(idxPt(0,i), idxPt(0,i1), idxPt(1,i1), idxPt(1,i));
        // outer wall
        emit_quad_as_tris(idxPt(2,i), idxPt(3,i), idxPt(3,i1), idxPt(2,i1));
        // bottom cap (annulus)
        emit_quad_as_tris(idxPt(0,i), idxPt(2,i), idxPt(2,i1), idxPt(0,i1));
        // top cap (annulus)
        emit_quad_as_tris(idxPt(1,i), idxPt(1,i1), idxPt(3,i1), idxPt(3,i));
    }

    std::fclose(f);
}

void dump_vtk(const std::vector<Particle>& P, int iter){
    char fn[64]; std::snprintf(fn,sizeof(fn),"atom.%d.vtk",iter);
    FILE* f = std::fopen(fn,"w");
    if(!f){ std::perror("vtk"); return; }
    std::fprintf(f,"# vtk DataFile Version 2.0\nCoax packing\nASCII\nDATASET POLYDATA\n");
    std::fprintf(f,"POINTS %zu float\n", P.size());
    for(const auto& a: P)
        std::fprintf(f,"%g %g %g\n", a.p.x, a.p.y, a.p.z);

    // velocity vectors
    std::fprintf(f,"POINT_DATA %zu\nVECTORS velocity float\n", P.size());
    for(const auto& a: P)
        std::fprintf(f,"%g %g %g\n", a.v.x, a.v.y, a.v.z);

    // radius scalar field
    std::fprintf(f,"SCALARS radius float 1\nLOOKUP_TABLE default\n");
    for(const auto& a: P)
        std::fprintf(f,"%g\n", a.r);

    // type_id as categorical integer field
    std::fprintf(f,"SCALARS type_id int 1\nLOOKUP_TABLE default\n");
    for(const auto& a: P)
        std::fprintf(f,"%d\n", (int)a.type_id);

    // optional label string (for visualization in ParaView)
    std::fprintf(f,"FIELD FieldData 1\n");
    std::fprintf(f,"element_symbol 1 %zu string\n", P.size());
    for(const auto& a: P)
        std::fprintf(f,"%s\n", element_symbol(a.type_id));

    std::fclose(f);
}


// Pick a particle type by g_type_frac[] (fractions; normalized)
inline int sample_type(std::mt19937& rng){
    std::uniform_real_distribution<real> U(0.0,1.0);
    real fsum = g_type_frac[0] + g_type_frac[1] + g_type_frac[2] + g_type_frac[3];
    if (fsum <= 0) return 0;
    real r = U(rng) * fsum;
    real acc = g_type_frac[0];
    if (r < acc) return 0;
    acc += g_type_frac[1];
    if (r < acc) return 1;
    acc += g_type_frac[2];
    if (r < acc) return 2;
    return 3;
}

// sample diameter from a given type's cumulative list (meters)
real sample_diameter_for_type(std::mt19937& rng, int type_id){
    const auto &dist = g_types[ (type_id>=0 && type_id<4) ? type_id : 0 ];
    const auto &D = dist.diam_um;
    const auto &C = dist.cumu;
    if (D.empty() || C.empty()) {
        // fallback to type 0 if misconfigured
        const auto &D0 = g_types[0].diam_um;
        const auto &C0 = g_types[0].cumu;
        std::uniform_real_distribution<real> U(0.0,1.0);
        real r = U(rng);
        size_t i=0; for(; i<C0.size(); ++i) if (r < C0[i]) break;
        if(i>=D0.size()) i = D0.size()-1;
        return 1e-6f * D0[i];
    }
    std::uniform_real_distribution<real> U(0.0,1.0);
    real r = U(rng);
    size_t i=0; for(; i<C.size(); ++i) if (r < C[i]) break;
    if(i>=D.size()) i = D.size()-1;
    return 1e-6f * D[i];
}

// (legacy wrapper – kept for compatibility; uses type 0)
real sample_diameter(std::mt19937& rng){
    return sample_diameter_for_type(rng, 0);
}

static real max_diam_from_types_um(){
    real mx = 0.0;
    for(int t=0;t<4;++t){
        for (real d : g_types[t].diam_um) if (d > mx) mx = d;
    }
    return mx > 0 ? mx : 5.0; // fallback 5 µm if misconfigured
}

static inline void successive_damping(std::vector<Particle>& P, real top_z)
{
    if (g_sd_sweeps <= 0) return;

    // do not overdamp fresh injectees just above the top plane
    real max_diam = 1e-6 * max_diam_from_types_um();
    max_diam += 2.0 * g_cushion;
    // account for effective collision diameter inflation from cushion
    const real z_guard  = top_z + g_sd_top_guardD * max_diam;
    const real z_guard_lo = top_z - 2.0 * g_sd_top_guardD * max_diam; // band *below* the lid

    for (int s = 0; s < g_sd_sweeps; ++s) {
        const real w     = (g_sd_sweeps > 1) ? (real(s) / real(g_sd_sweeps - 1)) : 0.0;
        const real alpha = g_sd_alpha0 + (g_sd_alpha1 - g_sd_alpha0) * w; // ramp up
        const real keep  = clamp_val<real>(1.0 - alpha, (real)0.0, (real)1.0);
#if defined(_OPENMP)
        #pragma omp parallel for
        for (int i = 0; i < (int)P.size(); ++i) {
            auto &a = P[i];
            if (a.asleep) continue;       // don’t touch sleepers
            if (a.p.z > z_guard_lo) continue; // spare freshly injected/free-fallers

            a.v.x *= keep; a.v.y *= keep; a.v.z *= keep;

            // extra tangential damping when touching / almost touching the floor
            const bool near_floor = (a.p.z <= ((a.r + g_cushion) + 1e-12));
            const bool near_top = (a.p.z >= z_guard);
            const real keep_xy = near_top ? keep : (near_floor ? keep * (1.0 - g_sd_tan_boost) : keep);
            const real keep_z  = keep;
            a.v.x *= clamp_val<real>(keep_xy, 0.0, 1.0);
            a.v.y *= clamp_val<real>(keep_xy, 0.0, 1.0);
            a.v.z *= clamp_val<real>(keep_z,  0.0, 1.0);

            // micro cutoff to kill jitter
            if (std::abs(a.v.x) < g_sd_vcut) a.v.x = 0.0;
            if (std::abs(a.v.y) < g_sd_vcut) a.v.y = 0.0;
            if (std::abs(a.v.z) < g_sd_vcut) a.v.z = 0.0;
        }
#else
        struct SD_Ctx { std::vector<Particle>* P; real keep, z_guard, z_guard_lo; };
        auto sd_worker = [](int i, void* vctx){
            SD_Ctx* c = (SD_Ctx*)vctx;
            auto &a = (*(c->P))[i];
            if (a.asleep) return;
            if (a.p.z > c->z_guard_lo) return;

            a.v.x *= c->keep; a.v.y *= c->keep; a.v.z *= c->keep;

            const bool near_floor = (a.p.z <= ((a.r + g_cushion) + 1e-12));
            const bool near_top = (a.p.z >= c->z_guard);
            const real keep_xy = near_top ? c->keep : (near_floor ? c->keep * (1.0 - g_sd_tan_boost) : c->keep);
            const real keep_z  = c->keep;
            a.v.x *= clamp_val<real>(keep_xy, 0.0, 1.0);
            a.v.y *= clamp_val<real>(keep_xy, 0.0, 1.0);
            a.v.z *= clamp_val<real>(keep_z,  0.0, 1.0);

            if (std::abs(a.v.x) < g_sd_vcut) a.v.x = 0.0;
            if (std::abs(a.v.y) < g_sd_vcut) a.v.y = 0.0;
            if (std::abs(a.v.z) < g_sd_vcut) a.v.z = 0.0;
        };
        SD_Ctx ctx{&P, keep, z_guard, z_guard_lo};
        for (int i = 0; i < (int)P.size(); ++i) sd_worker(i, &ctx);
#endif
    }
}

static inline real vertical_shake_accel(real t)
{
    if (g_shake_hz <= 0.0 || g_shake_amp == 0.0)
        return 0.0;

    // Base frequency
    const real w = 2.0 * M_PI * g_shake_hz;

    // Static RNG (thread-safe via thread_local)
    thread_local std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<real> noise_phase(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<real> noise_amp(0.8, 1.2);  // ±20% amplitude variation

    // Introduce slow phase drift and random amplitude per call
    const real random_phase = noise_phase(rng);
    const real amp = g_shake_amp * noise_amp(rng);

    // Optional: small frequency jitter (±5%)
    std::uniform_real_distribution<real> noise_freq(0.95, 1.05);
    const real wf = w * noise_freq(rng);

    return -(wf * wf) * amp * std::sin(wf * t + random_phase);
}


static inline void shake_accels(real t, real& ax, real& ay, real& az)
{
    ax = ay = az = 0.0;
    if (g_shake_hz <= 0.0) return;

    const real amp_z = (g_shake_amp_z != 0.0 ? g_shake_amp_z : g_shake_amp);
    const real w = 2.0 * M_PI * g_shake_hz;

    thread_local std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<real> noise_phase(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<real> noise_amp(0.8, 1.2);
    std::uniform_real_distribution<real> noise_freq(0.95, 1.05);

    const real random_phase = noise_phase(rng);
    const real wf = w * noise_freq(rng);

    auto accel_from_amp = [&](real amp)->real{
        if (amp == 0.0) return 0.0;
        const real a = amp * noise_amp(rng);
        return -(wf * wf) * a * std::sin(wf * t + random_phase);
    };

    az = accel_from_amp(amp_z);

    if (g_shake_xy_legacy) {
        ax = az;
        ay = az;
    } else {
        ax = accel_from_amp(g_shake_amp_x);
        ay = accel_from_amp(g_shake_amp_y);
    }
}


// uniform point in annulus cross-section for injection (x,y), at z ~ L + margin
inline void sample_injection_xy(std::mt19937& rng, real rp, real& x, real& y){
    std::uniform_real_distribution<real> U01(0.0, 1.0);

    // effective inner/outer radii available to the particle center
    const real rp_eff = rp + g_cushion;
    const real Rin_eff  = g_Rin  + rp_eff;
    const real Rout_eff = g_Rout - rp_eff;

    // if invalid (no room), clamp or throw
    if (Rout_eff <= Rin_eff){
        x = y = std::numeric_limits<real>::quiet_NaN();
        return;
    }

    // area-uniform in the annulus: sample r^2 uniformly
    const real Rin2  = Rin_eff  * Rin_eff;
    const real Rout2 = Rout_eff * Rout_eff;
    const real r     = std::sqrt( Rin2 + (Rout2 - Rin2) * U01(rng) );
    const real th    = (2.0*M_PI) * U01(rng);

    x = r * std::cos(th);
    y = r * std::sin(th);
}

// build grid indexing helpers
inline Vec3 worldMin(){ return g_minB; }
inline int  idx3(int ix,int iy,int iz){ return (iz*gy + iy)*gx + ix; }

void rebuild_grid(const std::vector<Particle>& P,
                  std::vector<int>& head,
                  std::vector<int>& next)
{
    std::fill(head.begin(), head.end(), -1);
    for (int i=0;i<(int)P.size();++i){
        const auto& a = P[i];
        int ix = (int)std::floor((a.p.x - g_minB.x)/g_cell_h);
        int iy = (int)std::floor((a.p.y - g_minB.y)/g_cell_h);
        int iz = (int)std::floor((a.p.z - g_minB.z)/g_cell_h);
        ix = clampi(ix,0,gx-1); iy = clampi(iy,0,gy-1); iz = clampi(iz,0,gz-1);
        int h = idx3(ix,iy,iz);
        next[i] = head[h];
        head[h] = i;
    }
}

// --------------------- Repulsive forces (optional) ---------------------
// Adds short-range repulsive accelerations when a particle is within
// g_repulse_range (in *gap* space) of another particle or a wall.
//
// Smooth ramp (C1):
//   s = clamp(1 - gap/range, 0..1)
//   if repulse_use_mass=1: a = (k/m) * (s*s) * n
//   else               : a = k * (s*s) * n   (k in [m/s^2])
// where n points away from the neighbor/wall.
static inline void apply_repulsion_accel(std::vector<Particle>& P,
                                        const std::vector<int>& head,
                                        const std::vector<int>& next,
                                        real top_z)
{
    if (g_repulse_range <= 0.0) return;
    const bool use_pp = (g_repulse_k_pp > 0.0);
    const bool use_pw = (g_repulse_k_pw > 0.0);
    if (!use_pp && !use_pw) return;

    const real eps = (real)1e-15;
    const real range = g_repulse_range;

    // Accumulate per-particle accelerations to avoid races.
    std::vector<Vec3> acc(P.size(), Vec3{0,0,0});

    if (use_pp) {
        const real kpp = g_repulse_k_pp;

#if defined(_OPENMP)
        #pragma omp parallel for schedule(static)
#endif
        for (int i = 0; i < (int)P.size(); ++i) {
            // compute i's cell
            int ix = clampi(int(std::floor((P[i].p.x - g_minB.x)/g_cell_h)), 0, gx-1);
            int iy = clampi(int(std::floor((P[i].p.y - g_minB.y)/g_cell_h)), 0, gy-1);
            int iz = clampi(int(std::floor((P[i].p.z - g_minB.z)/g_cell_h)), 0, gz-1);

            for (int dz=-1; dz<=1; ++dz){
                int jz = clampi(iz+dz,0,gz-1);
                for (int dy=-1; dy<=1; ++dy){
                    int jy = clampi(iy+dy,0,gy-1);
                    for (int dx=-1; dx<=1; ++dx){
                        int jx = clampi(ix+dx,0,gx-1);
                        int h = idx3(jx,jy,jz);
                        for (int j=head[h]; j!=-1; j=next[j]){
                            if (j<=i) continue;

                            // Optional: still apply repulsion to sleepers (it is gentle),
                            // but skip if both are asleep to avoid useless work.
                            if (P[i].asleep && P[j].asleep) continue;

                            Vec3 d = P[j].p - P[i].p;
                            real dist2 = dot(d,d);
                            if (dist2 <= eps) continue;
                            real dist = std::sqrt(dist2);

                            const real R = (P[i].r + g_cushion) + (P[j].r + g_cushion);
                            const real gap = dist - R;
                            if (gap <= 0.0) continue;          // overlap handled by collisions
                            if (gap >= range) continue;        // out of range

                            const real invd = (real)1.0 / dist;
                            const Vec3 n = d * invd;           // from i -> j
                            const real s = std::max((real)0.0, std::min((real)1.0, (real)1.0 - gap / range));
                            real a_scale_i = kpp * (s*s);
                            real a_scale_j = kpp * (s*s);
                            if (g_repulse_use_mass){
                                a_scale_i /= P[i].m;
                                a_scale_j /= P[j].m;
                            }
                            const Vec3 ai = n * (-a_scale_i);  // push i away from j
                            const Vec3 aj = n * ( +a_scale_j); // push j away from i

#if defined(_OPENMP)
                            #pragma omp atomic
#endif
                            acc[i].x += ai.x;
#if defined(_OPENMP)
                            #pragma omp atomic
#endif
                            acc[i].y += ai.y;
#if defined(_OPENMP)
                            #pragma omp atomic
#endif
                            acc[i].z += ai.z;

#if defined(_OPENMP)
                            #pragma omp atomic
#endif
                            acc[j].x += aj.x;
#if defined(_OPENMP)
                            #pragma omp atomic
#endif
                            acc[j].y += aj.y;
#if defined(_OPENMP)
                            #pragma omp atomic
#endif
                            acc[j].z += aj.z;
                        }
                    }
                }
            }
        }
    }

    if (use_pw) {
        const real kpw = g_repulse_k_pw;
#if defined(_OPENMP)
        #pragma omp parallel for
#endif
        for (int i = 0; i < (int)P.size(); ++i) {
            const Particle& a = P[i];
            const real r_eff = a.r + g_cushion;

            // Cylindrical walls
            const real rxy = std::hypot(a.p.x, a.p.y);
            if (rxy > eps) {
                const real invr = (real)1.0 / rxy;
                const Vec3 rhat = { a.p.x * invr, a.p.y * invr, 0.0 };

                // inner wall: r >= Rin + r_eff (push outward)
                {
                    const real rmin = g_Rin + r_eff;
                    const real gap = rxy - rmin;
                    if (gap > 0.0 && gap < range) {
                        const real s = std::max((real)0.0, std::min((real)1.0, (real)1.0 - gap / range));
                        real a_scale = kpw * (s*s);
                        if (g_repulse_use_mass) a_scale /= a.m;
                        acc[i].x += rhat.x * (+a_scale);
                        acc[i].y += rhat.y * (+a_scale);
                    }
                }
                // outer wall: r <= Rout - r_eff (push inward)
                {
                    const real rmax = g_Rout - r_eff;
                    const real gap = rmax - rxy;
                    if (gap > 0.0 && gap < range) {
                        const real s = std::max((real)0.0, std::min((real)1.0, (real)1.0 - gap / range));
                        real a_scale = kpw * (s*s);
                        if (g_repulse_use_mass) a_scale /= a.m;
                        acc[i].x += rhat.x * (-a_scale);
                        acc[i].y += rhat.y * (-a_scale);
                    }
                }
            }

            // bottom plane: z >= r_eff (push up)
            {
                const real gap = a.p.z - r_eff;
                if (gap > 0.0 && gap < range) {
                    const real s = std::max((real)0.0, std::min((real)1.0, (real)1.0 - gap / range));
                    real a_scale = kpw * (s*s);
                        if (g_repulse_use_mass) a_scale /= a.m;
                    acc[i].z += (+a_scale);
                }
            }
            // top plane: z <= top_z - r_eff (push down)
            {
                const real zmax = top_z - r_eff;
                const real gap = zmax - a.p.z;
                if (gap > 0.0 && gap < range) {
                    const real s = std::max((real)0.0, std::min((real)1.0, (real)1.0 - gap / range));
                    real a_scale = kpw * (s*s);
                        if (g_repulse_use_mass) a_scale /= a.m;
                    acc[i].z += (-a_scale);
                }

            // Corner anti-trap: outer-wall with bottom/top plane
            // Helps prevent particles from pinning at the OD edge (r ~ Rout-r, z ~ r or z ~ top_z-r).
            if (rxy > eps) {
                const real invr = (real)1.0 / rxy;
                const Vec3 rhat = { a.p.x * invr, a.p.y * invr, 0.0 };

                const real gap_out = (g_Rout - r_eff) - rxy;     // >0 means inside domain
                const real gap_bot = a.p.z - r_eff;              // >0 above bottom plane
                const real gap_top = (top_z - r_eff) - a.p.z;    // >0 below top plane

                auto add_corner = [&](real gap_r, real gap_z, real nz){
                    if (gap_r > 0.0 && gap_z > 0.0 && gap_r < range && gap_z < range) {
                        // blend strength based on proximity to both surfaces
                        const real sr = std::max((real)0.0, std::min((real)1.0, (real)1.0 - gap_r / range));
                        const real sz = std::max((real)0.0, std::min((real)1.0, (real)1.0 - gap_z / range));
                        const real s  = sr * sz;
                        // direction away from corner: inward + up (bottom) or inward + down (top)
                        Vec3 n = { -rhat.x, -rhat.y, nz };
                        const real nn = std::sqrt(n.x*n.x + n.y*n.y + n.z*n.z);
                        if (nn > eps) { n.x /= nn; n.y /= nn; n.z /= nn; }
                        real a_scale = kpw * (s*s);
                        if (g_repulse_use_mass) a_scale /= a.m;
                        acc[i].x += n.x * a_scale;
                        acc[i].y += n.y * a_scale;
                        acc[i].z += n.z * a_scale;
                    }
                };

                // bottom outer corner: inward + up
                add_corner(gap_out, gap_bot, +1.0);
                // top outer corner: inward + down
                add_corner(gap_out, gap_top, -1.0);
            }
            }
        }
    }

    // Apply accumulated accelerations
#if defined(_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i < (int)P.size(); ++i) {
        if (!P[i].asleep) {
            {
            Vec3 aadd = acc[i];
            if (g_repulse_dvmax > 0.0){
                // Limit per-step delta-v from repulsion to keep things stable with tiny particle masses
                const real dvmax = g_repulse_dvmax;
                const real amag = std::sqrt(aadd.x*aadd.x + aadd.y*aadd.y + aadd.z*aadd.z);
                if (amag > 0.0){
                    const real dv = amag * g_dt;
                    if (dv > dvmax){
                        const real s = dvmax / dv;
                        aadd.x *= s; aadd.y *= s; aadd.z *= s;
                    }
                }
            }
            P[i].a.x += aadd.x;
            P[i].a.y += aadd.y;
            P[i].a.z += aadd.z;
        }
        }
    }
}

// -------------------------- Collisions -----------------------------
inline void collide_pair(Particle& A, Particle& B, real e_common, real tdamp)
{
    // If both asleep, nothing to do (but allow “wake” below if hit hard enough)
    if (A.asleep && B.asleep) return;

    const real eps = (real)1e-12;

    // Geometry
    Vec3 d = B.p - A.p;
    real dist2 = dot(d, d);
    const real R = (A.r + g_cushion) + (B.r + g_cushion);
    if (dist2 <= 0.0) dist2 = 0.0;
    if (dist2 == 0.0) {
        // Perfect overlap: pick a stable direction from relative velocity or +x
        Vec3 vrel0 = B.v - A.v;
        real vr = norm(vrel0);
        d = (vr > eps ? vrel0 * ( (real)1.0 / vr ) : Vec3{1,0,0});
        dist2 = dot(d, d);
    }
    real dist = std::sqrt(dist2);
    if (dist >= R) return; // no contact

    // Contact normal (from A to B), unit length
    const real invdist = (dist > eps ? (real)1.0 / dist : (real)0.0);
    const Vec3 n = d * invdist;

    // Penetration depth and mass-weighted positional correction (stabilized)
    const real pen = (R - dist);
    if (pen > 0.0f)
    {
        const real sep_cap = g_pair_sep_cap * R; // cap to avoid teleports
        const real sep = (pen < sep_cap ? pen : sep_cap);

        // Mass weights (heavier body moves less)
        const real invmA = (real)1.0 / A.m;
        const real invmB = (real)1.0 / B.m;
        const real wA = invmA / (invmA + invmB);
        const real wB = invmB / (invmA + invmB);

        // Add mild Baumgarte bias: move a bit more than purely geometric sep
        const real sep_bias = sep * ( (real)1.0 + g_pair_beta );

        A.p = A.p - n * ( wA * sep_bias );
        B.p = B.p + n * ( wB * sep_bias );
    }

    // Relative velocity (post position projection)
    Vec3 vrel = B.v - A.v;
    real vn = dot(vrel, n);

    // Wake up on sufficiently strong approach
    if (vn < 0.0) {
        if (A.asleep && (-vn) > g_wake_vn) { A.asleep = false; A.calm_steps = 0; }
        if (B.asleep && (-vn) > g_wake_vn) { B.asleep = false; B.calm_steps = 0; }
    }

    // If separating, no impulses needed
    if (vn >= 0.0) return;

    // Effective mass for impulses
    const real invmA = (real)1.0 / A.m;
    const real invmB = (real)1.0 / B.m;
    const real m_eff = (real)1.0 / (invmA + invmB);

    // ---------- Normal impulse (restitution) ----------
    // Use caller-provided restitution e_common (0..1). Sticky threshold handled by vstatic below.
    const real Jn_mag = -( (real)1.0 + e_common ) * vn * m_eff; // positive magnitude
    const Vec3 Jn = n * Jn_mag;

    A.v = A.v - Jn * invmA;
    B.v = B.v + Jn * invmB;

    // Recompute relative velocity after normal impulse; split into components
    vrel = B.v - A.v;
    vn   = dot(vrel, n);
    Vec3 vt = vrel - n * vn;
    real vt_mag = norm(vt);

    // ---------- Tangential friction (Coulomb, static then kinetic) ----------
    if (vt_mag > (real)0.0)
    {
        // Desired tangential impulse to eliminate slip completely (static target):
        // Jt_des = - m_eff * vt
        Vec3 Jt_des = vt * ( -m_eff );

        // Static friction cap: up to mu_s * |Jn|
        const real Jt_static_cap = g_pair_mu_s * std::abs(Jn_mag);
        const real Jt_des_mag    = norm(Jt_des);

        if (Jt_des_mag <= Jt_static_cap)
        {
            // Static: cancel slip fully
            A.v = A.v - Jt_des * invmA;
            B.v = B.v + Jt_des * invmB;
        }
        else
        {
            // Kinetic: apply capped impulse opposing slip
            const real Jt_kin_cap = g_pair_mu_k * std::abs(Jn_mag);
            Vec3 t_hat = Jt_des * ( (real)1.0 / (Jt_des_mag + eps) ); // = -vt/|vt|
            Vec3 Jt = t_hat * Jt_kin_cap; // already points to reduce slip
            A.v = A.v - Jt * invmA;
            B.v = B.v + Jt * invmB;
        }
    }

    // ---------- Static “snap” for tiny residual slip ----------
    // (Keeps packs from slowly shearing at mm/s forever)
    vrel = B.v - A.v;
    vn   = dot(vrel, n);
    vt   = vrel - n * vn;
    vt_mag = norm(vt);

    if (vt_mag < g_pair_vstatic)
    {
        // Impulse to zero tangential slip exactly
        Vec3 Jsnap = vt * ( -m_eff );
        A.v = A.v - Jsnap * invmA;
        B.v = B.v + Jsnap * invmB;
    }

    // ---------- Successive-damping hook on tangential slip (if tdamp>0) ----------
    // Reduce remaining slip by a fractional amount tdamp (0..1) via impulse form.
    if (tdamp > (real)0.0)
    {
        vrel = B.v - A.v;
        vn   = dot(vrel, n);
        vt   = vrel - n * vn;
        vt_mag = norm(vt);

        if (vt_mag > (real)0.0)
        {
            const real damp = std::max((real)0.0, std::min((real)1.0, tdamp));
            Vec3 Jd = vt * ( -m_eff * damp );  // impulse that shrinks slip by ~damp
            A.v = A.v - Jd * invmA;
            B.v = B.v + Jd * invmB;
        }
    }
}


inline void collide_walls(Particle& a, real top_z, real e_common)
{
    const real eps = (real)1e-12;
    const real TWO_PI = (real)6.2831853071795864769;

auto apply_wall_spring = [&](const Vec3& n, real pen, real wall_vn)
{
    if (g_wall_k <= 0.0 || pen <= 0.0) return;

    const real vn_rel = dot(a.v, n) - wall_vn;
    const real k = g_wall_k;
    const real c = 2.0 * g_wall_zeta * std::sqrt(std::max((real)0.0, k * a.m));
    const real an = (k * pen - c * vn_rel) / a.m;

    real dvn = an * g_dt;
    if (dvn >  g_wall_dvmax) dvn =  g_wall_dvmax;
    if (dvn < -g_wall_dvmax) dvn = -g_wall_dvmax;

    a.v = a.v + n * dvn;
};


    auto unit_rhat = [&](real x, real y) -> std::pair<real,real> {
        real r = std::hypot(x, y);
        if (r > eps) return {x / r, y / r};
        // fallback: use velocity dir if available, else +x
        real vxy = std::hypot(a.v.x, a.v.y);
        if (vxy > eps) return {a.v.x / vxy, a.v.y / vxy};
        return { (real)1.0, (real)0.0 };
    };

    // Deterministic geometric roughness for cylindrical walls.
    // Returns an additive radial offset (meters) applied to the *inner* wall radius
    // and subtracted from the *outer* wall radius (so positive values intrude into
    // the domain from both sides).
    auto wall_rough = [&](real x, real y, real z) -> real {
        if (g_wall_rough_amp <= (real)0.0) return (real)0.0;
        const int mth = std::max(1, g_wall_rough_mth);
        const int mz  = std::max(0, g_wall_rough_mz);
        real th = std::atan2(y, x);
        // keep theta continuous in [0,2pi)
        if (th < (real)0.0) th += TWO_PI;

        real f = std::sin((real)mth * th + g_wall_rough_phth);
        if (mz > 0 && g_L > eps) {
            f *= std::sin((real)mz * (TWO_PI * (z / g_L)) + g_wall_rough_phz);
        }
        return g_wall_rough_amp * f;
    };

    auto resolve_contact = [&](const Vec3& n, real pen,
                               real e_n,     // normal restitution for this contact
                               real beta,    // Baumgarte fraction
                               real mu_t,    // tangential friction coeff
                               real roll,    // tangential rolling loss
                               real v_static // static snap threshold (m/s)
                              )
    {
        // 1) project position out of penetration
        a.p.x += n.x * pen;
        a.p.y += n.y * pen;
        a.p.z += n.z * pen;

        // 2) split velocity into normal + tangential
        real vn = a.v.x*n.x + a.v.y*n.y + a.v.z*n.z;

        // 3) normal response (only if approaching the surface)
        if (vn < 0.0)
        {
            // Baumgarte-like extra normal damping proportional to penetration depth
            const real cn = beta / g_dt;                 // [1/s]
            real v_rebound = -e_n * vn - cn * pen;       // in normal direction

            // stick tiny rebounds at the floor/top as needed (optional)
            // (we keep this general here; your floor-specific vstick still applies below)
            // a.v += (v_rebound - vn) * n;
            real dvn = (v_rebound - vn);
            a.v.x += dvn * n.x;
            a.v.y += dvn * n.y;
            a.v.z += dvn * n.z;
        }
        else
        {
            // If numerics gave a tiny outward vn, no normal action needed
        }

        // 4) tangential friction (Coulomb-like) + mild rolling loss
        //   vt = v - (v·n) n  (recompute vn after possible normal change)
        real vn2 = a.v.x*n.x + a.v.y*n.y + a.v.z*n.z;
        Vec3 vt  = { a.v.x - vn2*n.x, a.v.y - vn2*n.y, a.v.z - vn2*n.z };
        real vt_mag = std::sqrt(vt.x*vt.x + vt.y*vt.y + vt.z*vt.z);

        if (vt_mag > (real)0.0)
        {
            // Scale slip by a normal-impulse proxy; use |vn_before| if you want impact-based,
            // here we use |vn2| to remain stable even on gentle sustained contacts.
            real slip_scale = (real)1.0 - mu_t * (1.0 + e_n) * std::max((real)0.0, -vn2) / (vt_mag + eps);
            if (slip_scale < (real)0.0) slip_scale = (real)0.0; // stick limit
            // apply slip scaling
            a.v.x = vn2*n.x + vt.x * slip_scale;
            a.v.y = vn2*n.y + vt.y * slip_scale;
            a.v.z = vn2*n.z + vt.z * slip_scale;

            // mild rolling resistance
            a.v.x *= (1.0 - roll);
            a.v.y *= (1.0 - roll);
            a.v.z *= (1.0 - roll);
        }

        // 5) static snap (if tangential nearly stopped, kill it exactly)
        vn2 = a.v.x*n.x + a.v.y*n.y + a.v.z*n.z;
        vt  = { a.v.x - vn2*n.x, a.v.y - vn2*n.y, a.v.z - vn2*n.z };
        vt_mag = std::sqrt(vt.x*vt.x + vt.y*vt.y + vt.z*vt.z);
        if (vt_mag < v_static)
        {
            // keep normal component, zero tangential exactly
            a.v.x = vn2*n.x;
            a.v.y = vn2*n.y;
            a.v.z = vn2*n.z;
        }
    };

    const real r_eff = a.r + g_cushion;

    // Clamp roughness amplitude to avoid pathological inversion of the annulus.
    // This is conservative: allow up to 25% of the local free gap.
    const real free_gap = (g_Rout - g_Rin) - (real)2.0 * r_eff;
    const real rough_cap = std::max((real)0.0, (real)0.25 * free_gap);

    // ---------------------- INNER CYLINDER (r = g_Rin + a.r) ----------------------
    {
        real rough = wall_rough(a.p.x, a.p.y, a.p.z);
        if (rough >  rough_cap) rough =  rough_cap;
        if (rough < -rough_cap) rough = -rough_cap;
        const real rmin = g_Rin + r_eff + rough;
        real rxy = std::hypot(a.p.x, a.p.y);
        if (rxy < rmin)
        {
            const real pen = (rmin - rxy);            // >0 into wall
            auto [rx, ry] = unit_rhat(a.p.x, a.p.y);  // r̂
            Vec3 n = { +rx, +ry, 0.0 };               // normal points outward (toward allowed region)
            apply_wall_spring(n, pen, (real)0.0);
            // Snap radius exactly to rmin by using position projection in resolve_contact
            resolve_contact(n, pen, /*e_n=*/e_common,
                            /*beta=*/g_wall_beta, /*mu_t=*/g_wall_mu_t,
                            /*roll=*/g_wall_roll, /*v_static=*/(real)0.02);
            // ensure exactly on the surface to avoid drift
            const real s = (rmin + eps) / std::hypot(a.p.x, a.p.y);
            a.p.x *= s; a.p.y *= s;
        }
    }

    // ---------------------- OUTER CYLINDER (r = g_Rout - a.r) --------------------
    {
        real rough = wall_rough(a.p.x, a.p.y, a.p.z);
        if (rough >  rough_cap) rough =  rough_cap;
        if (rough < -rough_cap) rough = -rough_cap;
        const real rmax = g_Rout - r_eff - rough;
        real rxy = std::hypot(a.p.x, a.p.y);
        if (rxy > rmax)
        {
            const real pen = (rxy - rmax);            // >0 into wall
            auto [rx, ry] = unit_rhat(a.p.x, a.p.y);  // r̂
            Vec3 n = { -rx, -ry, 0.0 };               // normal points inward (toward allowed region)
            apply_wall_spring(n, pen, (real)0.0);
            resolve_contact(n, pen, /*e_n=*/e_common,
                            /*beta=*/g_wall_beta, /*mu_t=*/g_wall_mu_t,
                            /*roll=*/g_wall_roll, /*v_static=*/(real)0.02);
            // snap radius exactly to rmax
            const real s = (rmax - eps) / std::hypot(a.p.x, a.p.y);
            a.p.x *= s; a.p.y *= s;
        }
    }

    // ---------------------- BOTTOM PLANE (z = a.r) --------------------------------
    {
        const real zbot = r_eff;
        if (a.p.z < zbot)
        {
            const real pen = (zbot - a.p.z);
            Vec3 n = { 0.0, 0.0, +1.0 };              // up, toward allowed region
            apply_wall_spring(n, pen, (real)0.0);
            resolve_contact(n, pen, /*e_n=*/g_e_bottom,
                            /*beta=*/g_floor_beta, /*mu_t=*/g_floor_mu_t,
                            /*roll=*/g_floor_roll, /*v_static=*/g_vstatic_floor);

            // belt-and-suspenders floor clamp
            if (a.p.z < zbot) a.p.z = zbot;
            if (a.v.z < 0)    a.v.z = 0;
        }
    }

    // ---------------------- TOP PISTON (z = top_z - a.r) --------------------------
    {
        const real ztop = top_z - r_eff;
        if (a.p.z > ztop)
        {
            const real pen = (a.p.z - ztop);
            Vec3 n = { 0.0, 0.0, -1.0 };              // down, toward allowed region
            apply_wall_spring(n, pen, g_piston_wall_vn);
            resolve_contact(n, pen, /*e_n=*/e_common,
                            /*beta=*/g_wall_beta, /*mu_t=*/g_wall_mu_t,
                            /*roll=*/g_wall_roll, /*v_static=*/(real)0.02);

            // exact snap to the plane
            if (a.p.z > ztop) a.p.z = ztop;
            if (a.v.z > 0)    a.v.z = 0;
        }
    }
}


// ----------------------- Integration Loop -------------------------
int main(int argc, char** argv){
    auto usage = [&](int rc){
        std::fprintf(stderr,
            "Usage (legacy positional):\n"
            "  %s natoms_max dt niter dump_interval debug seed Rin Rout L flux g shake_hz shake_amp fill_time ram_start ram_duration ram_speed phi_target\n"
            "\n"
            "Usage (flags, order-independent):\n"
            "  %s [--natoms_max N] [--dt DT] [--niter N] [--dump_interval N]\n"
            "     [--debug 0|1] [--seed S] [--rin R] [--rout R] [--length L]\n"
            "     [--flux N] [--gravity G] [--shake_hz HZ] [--shake_amp A]\n"
            "     [--fill_time T] [--ram_start T0] [--ram_duration DT] [--ram_speed V]\n"
            "     [--phi_target PHI]\n"
            "     [--repulse_range R] [--repulse_k_pp K] [--repulse_k_pw K] [--repulse_use_mass 0|1] [--repulse_dvmax V]\n"
            "     [--stop_vrms V] [--stop_vmax V] [--stop_sleep_frac F]\n"
            "     [--stop_check_interval N] [--stop_checks_required N]\n"
            "     [--lin_damp D] [--e_pp E] [--e_pw E] [--tangent_damp D]\n"

            "     [--wall_rough_amp A] [--wall_rough_mth M] [--wall_rough_mz M]\n"
            "     [--threads N]\n"
            "     [--xyz_interval N] [--vtk_interval N] [--vtk_domain_interval N] [--vtk_domain_segments N]\n",
            argv[0], argv[0]);
        return rc;
    };

    // -------------------- defaults (match batch-script style) --------------------
    g_debug = 0;
    g_seed  = 42;
    g_natoms_max = 25000;
    g_dt = 2e-6;
    g_niter = 100000;
    g_dump_interval = 500;

    g_Rin = 22e-6;
    g_Rout = 40e-6;
    g_L = 380e-6;

    g_flux = 60000;
    g_g = 9.81;
    g_shake_hz = 1000.0;
    g_shake_amp = 5e-6;

    g_shake_amp_x = 0.0;
    g_shake_amp_y = 0.0;
    g_shake_amp_z = 0.0;
    g_shake_xy_legacy = 1;

    g_cushion = 0.0;
    g_wall_k = 0.0;
    g_wall_zeta = 0.20;
    g_wall_dvmax = 5.0;

    g_lin_damp = 0.0; // [1/s]

    g_repulse_range = 0.0;
    g_repulse_k_pp  = 0.0;
    g_repulse_k_pw  = 0.0;
    g_repulse_use_mass = 0;
    g_repulse_dvmax = 0.0;

    g_wall_rough_amp = 0.0;
    g_wall_rough_mth = 8;
    g_wall_rough_mz  = 3;
    g_wall_rough_phth = 0.0;
    g_wall_rough_phz  = 0.0;

    g_inject_vx = 0.0;
    g_inject_vy = 0.0;
    g_inject_vz = -0.5;

    g_fill_time = 8.0;
    g_ram_t0 = 0.0;
    g_ram_dt = 0.0;
    g_ram_speed = 0.0;

    g_phi_target = 0.0;

    g_stop_vrms = 0.0;
    g_stop_vmax = 0.0;
    g_stop_sleep_frac = 0.0;
    g_stop_check_interval = 200;
    g_stop_checks_required = 10;

    g_threads = 0;
    g_xyz_interval = -1;
    g_vtk_interval = -1;
    g_vtk_domain_interval = 0;
    g_vtk_domain_segments = 96;

    if (argc <= 1) return usage(1);

    // -------------------- parse --------------------
    std::vector<std::string> pos;
    std::unordered_map<std::string, std::string> kv;

    auto is_flag = [&](const char* s)->bool{
        return s && s[0]=='-' && s[1]=='-';
    };
    auto to_lower = [&](std::string s)->std::string{
        for(char& c: s) c = (char)std::tolower((unsigned char)c);
        return s;
    };

    for (int i=1; i<argc; ++i){
        const char* a = argv[i];
        if (!a) continue;
        if (std::strcmp(a, "--help")==0 || std::strcmp(a, "-h")==0){
            return usage(0);
        }
        if (is_flag(a)){
            std::string key = a+2;
            std::string val = "1";
            auto eq = key.find('=');
            if (eq != std::string::npos){
                val = key.substr(eq+1);
                key = key.substr(0, eq);
            } else {
                // if next token exists and is not another flag, treat it as the value
                if (i+1 < argc && !(argv[i+1][0]=='-' && argv[i+1][1]=='-')){
                    val = argv[++i];
                }
            }
            kv[to_lower(key)] = val;
        } else {
            pos.emplace_back(a);
        }
    }

    auto get_i = [&](const char* k, int& dst){
        auto it = kv.find(k);
        if (it != kv.end()) dst = std::atoi(it->second.c_str());
    };
    auto get_r = [&](const char* k, real& dst){
        auto it = kv.find(k);
        if (it != kv.end()) dst = (real)std::atof(it->second.c_str());
    };

    // Legacy positional mode if any positional args were provided.
    // This preserves old scripts like: coax_pack_cpu.exe 20000 2e-6 ...
    if (!pos.empty()){
        if ((int)pos.size() < 18){
            std::fprintf(stderr,"Error: legacy positional CLI expects 18 values (got %zu)\n", pos.size());
            return usage(1);
        }
        int iarg=0;
        g_natoms_max    = std::atoi(pos[iarg++].c_str());
        g_dt            = std::atof(pos[iarg++].c_str());
        g_niter         = std::atoi(pos[iarg++].c_str());
        g_dump_interval = std::atoi(pos[iarg++].c_str());
        g_debug         = std::atoi(pos[iarg++].c_str());
        g_seed          = std::atoi(pos[iarg++].c_str());
        g_Rin           = std::atof(pos[iarg++].c_str());
        g_Rout          = std::atof(pos[iarg++].c_str());
        g_L             = std::atof(pos[iarg++].c_str());
        g_flux          = std::atoi(pos[iarg++].c_str());
        g_g             = std::abs(std::atof(pos[iarg++].c_str()));
        g_shake_hz      = std::atof(pos[iarg++].c_str());
        g_shake_amp     = std::atof(pos[iarg++].c_str());
        g_fill_time     = std::atof(pos[iarg++].c_str());
        g_ram_t0        = std::atof(pos[iarg++].c_str());
        g_ram_dt        = std::atof(pos[iarg++].c_str());
        g_ram_speed     = std::atof(pos[iarg++].c_str());
        g_phi_target    = std::atof(pos[iarg++].c_str()); // <=0 disables cutoff
    } else {
        // Flag mode (order-independent)
        get_i("natoms_max", g_natoms_max);
        get_r("dt", g_dt);
        get_i("niter", g_niter);
        get_i("dump_interval", g_dump_interval);
        get_i("debug", g_debug);
        get_i("seed", g_seed);
        get_r("rin", g_Rin);
        get_r("rout", g_Rout);
        get_r("length", g_L);
        get_i("flux", g_flux);
        get_r("gravity", g_g);
        get_r("shake_hz", g_shake_hz);
        get_r("shake_amp", g_shake_amp);
        get_r("shake_amp_x", g_shake_amp_x);
        get_r("shake_amp_y", g_shake_amp_y);
        get_r("shake_amp_z", g_shake_amp_z);
        get_i("shake_xy_legacy", g_shake_xy_legacy);

        get_r("cushion", g_cushion);
        get_r("wall_k", g_wall_k);
        get_r("wall_zeta", g_wall_zeta);
        get_r("wall_dvmax", g_wall_dvmax);

        get_r("lin_damp", g_lin_damp);
        get_r("e_pp", g_e_pp);
        get_r("e_pw", g_e_pw);
        get_r("tangent_damp", g_tangent_damp);

        get_r("wall_rough_amp", g_wall_rough_amp);
        get_i("wall_rough_mth", g_wall_rough_mth);
        get_i("wall_rough_mz",  g_wall_rough_mz);

        get_r("inject_vx", g_inject_vx);
        get_r("inject_vy", g_inject_vy);
        get_r("inject_vz", g_inject_vz);
        get_r("fill_time", g_fill_time);
        get_r("ram_start", g_ram_t0);
        get_r("ram_duration", g_ram_dt);
        get_r("ram_speed", g_ram_speed);
        get_r("phi_target", g_phi_target);

        get_r("repulse_range", g_repulse_range);
        get_r("repulse_k_pp",  g_repulse_k_pp);
        get_r("repulse_k_pw",  g_repulse_k_pw);
        get_i("repulse_use_mass", g_repulse_use_mass);
        get_r("repulse_dvmax", g_repulse_dvmax);

        get_r("stop_vrms", g_stop_vrms);
        get_r("stop_vmax", g_stop_vmax);
        get_r("stop_sleep_frac", g_stop_sleep_frac);
        get_i("stop_check_interval", g_stop_check_interval);
        get_i("stop_checks_required", g_stop_checks_required);

        get_i("threads", g_threads);
        get_i("nthreads", g_threads);

        get_i("xyz_interval", g_xyz_interval);
        get_i("vtk_interval", g_vtk_interval);
        get_i("vtk_domain_interval", g_vtk_domain_interval);
        get_i("vtk_domain_segments", g_vtk_domain_segments);

        // accept some aliases from older batch scripts
        get_r("r_in", g_Rin);
        get_r("r_out", g_Rout);
        get_r("shake_ampx", g_shake_amp_x);
        get_r("shake_ampy", g_shake_amp_y);
        get_r("shake_ampz", g_shake_amp_z);
    }

    // If roughness is enabled, derive phases from the RNG seed so runs are repeatable.
    // (User does not need to supply phases; changing seed changes the texture.)
    if (g_wall_rough_amp > (real)0.0) {
        const double TWO_PI = 6.2831853071795864769;
        uint32_t s = (uint32_t)g_seed;
        auto next_u01 = [&]()->double{
            // LCG (Numerical Recipes)
            s = 1664525u * s + 1013904223u;
            return (double)(s) / (double)UINT32_MAX;
        };
        g_wall_rough_phth = (real)(TWO_PI * next_u01());
        g_wall_rough_phz  = (real)(TWO_PI * next_u01());
        g_wall_rough_mth = std::max(1, g_wall_rough_mth);
        g_wall_rough_mz  = std::max(0, g_wall_rough_mz);
    }

#if defined(_OPENMP)
    if (g_threads > 0) omp_set_num_threads(g_threads);
#endif

    // resolve output intervals (default: follow dump_interval)
    if (g_xyz_interval < 0) g_xyz_interval = g_dump_interval;
    if (g_vtk_interval < 0) g_vtk_interval = g_dump_interval;

    real phi = 0.0; // volume fraction estimate

    // init multi-type defaults (copies type 0 into 1..3 unless you override)
init_default_distributions_once();
    // normalize fractions once; allow user to set g_type_frac[] later
    real fsum = g_type_frac[0] + g_type_frac[1] + g_type_frac[2] + g_type_frac[3];
    if (fsum <= 0) { g_type_frac[0]=1.0; g_type_frac[1]=g_type_frac[2]=g_type_frac[3]=0.0; fsum=1.0; }

    if (g_seed<=0) g_seed = int(std::time(nullptr));
    std::mt19937 rng(g_seed);

    // domain bbox for neighbor grid
    real max_diam = 1e-6 * max_diam_from_types_um();
    max_diam += 2.0 * g_cushion;
    // account for effective collision diameter inflation from cushion
    g_cell_h = std::max(max_diam, 0.5*max_diam) * 1.05; // start near 1.05*D

    g_minB = {-g_Rout - max_diam, -g_Rout - max_diam, -max_diam};
    g_maxB = { g_Rout + max_diam,  g_Rout + max_diam,  g_L + max_diam};

    // Cap total number of cells ~ O(natoms_max)
    auto spanX = (g_maxB.x - g_minB.x);
    auto spanY = (g_maxB.y - g_minB.y);
    auto spanZ = (g_maxB.z - g_minB.z);

    // Budget ~ 8 * natoms_max cells (tweak as desired)
    double cell_budget = std::max( (double) (8LL * std::max(1,g_natoms_max)), 1e5 ); // also avoid too tiny
    double cell_vol = (spanX*spanY*spanZ) / cell_budget;
    double h_target = std::cbrt(std::max(cell_vol, 1e-18));

    // final cell size: at least max_diam, and at least h_target
    g_cell_h = std::max( g_cell_h, h_target );

    gx = std::max(4, int(std::ceil(spanX/g_cell_h)));
    gy = std::max(4, int(std::ceil(spanY/g_cell_h)));
    gz = std::max(8, int(std::ceil(spanZ/g_cell_h)));

    std::vector<Particle> P;
    P.reserve(std::min(g_natoms_max, 200000));

    // neighbor grid buffers
    std::vector<int> head(gx*gy*gz,-1);
    std::vector<int> next(g_natoms_max, -1);

    // domain volume (constant) for phi
    const real domVol = M_PI * (g_Rout*g_Rout - g_Rin*g_Rin) * g_L;
    // running particle-volume accumulator (O(1) update on injection)
    real totalVol_running = 0.0;

    // piston kinematics
    auto topPlane = [&](real t)->real{
        if (t < g_ram_t0) return g_L;
        real tau = std::min(std::max(t - g_ram_t0, real(0)), g_ram_dt);
        real top = g_L - g_ram_speed * tau;

        const real eps = 1e-12;
        if (top < (0.0 + eps)) top = (0.0 + eps);
        return top;
    };

    // injection book-keeping
    int injected_total = 0;
    real t = 0.0;

    std::printf("Coax pack: Rin=%.4g m Rout=%.4g m L=%.4g m, fill_time=%.2fs, flux=%d 1/s\n",
                g_Rin, g_Rout, g_L, g_fill_time, g_flux);

    for (int it=0; it<=g_niter; ++it, t += g_dt){
        // 1) Inject new spheres while t < fill_time
        // 0) If phi target reached, disable injection and switch to zero-g "settle"
        if (!g_phi_reached && g_phi_target > 0.0) {
            // phi is maintained via totalVol_running; cheap to check every step
           phi = (domVol > 0.0 ? totalVol_running / domVol : 0.0);
            if (phi >= g_phi_target) {
                g_phi_reached = true;
                g_g_run = 0.0;         // zero gravity
                //g_shake_hz = 0.0;      // stop shaking
                g_ram_speed = 0.0;     // freeze top lid motion
                g_ram_dt = 0.0;
                // wake everyone once so they can relax to rest under zero-g + damping
                #pragma omp parallel for
                for (int i = 0; i < (int)P.size(); ++i) {
                    P[i].asleep = false;
                    P[i].calm_steps = 0;
                }
                if (g_debug) std::printf("[phi-cutoff] target=%.5f reached at t=%.6g (N=%zu)\n",
                                          g_phi_target, t, P.size());
            }
        }

        // 1) Inject new spheres while below time limit AND below phi target (if set)
        if (!g_phi_reached && t < g_fill_time && injected_total < g_natoms_max && g_flux>0){
            // expected total injected by now:
            int should = int(std::floor(t * g_flux));
            int want = std::min(should - injected_total, g_natoms_max - injected_total);
            if (it==0) want = std::max(want,1); // ensure at least one
            if (want > 0){
                int oldN = (int)P.size();
                P.resize(oldN + want);
                // Per-thread RNGs to avoid races and XY bias
#if defined(_OPENMP)
                #pragma omp parallel
                {
                    const int tid = omp_get_thread_num();
                    std::mt19937 trng((unsigned)g_seed + 0x9E3779B9u * (unsigned)tid + (unsigned)it*101u);
                    #pragma omp for
                    for (int k=0;k<want;++k){
                        // resample until the diameter fits the annulus
                        real d; real r; real x; real y;
                        int  t_id;
                        for(;;){
                            t_id = sample_type(trng);
                            d = sample_diameter_for_type(trng, t_id);
                            r = 0.5*d;
                            sample_injection_xy(trng, r, x, y);
                            if (x==x && y==y) { // not NaN => valid
                                break;
                            }
                        }
                        real z = g_L + 2*r;
                        Particle a;
                        a.p = {x,y,z};
                        a.v = {g_inject_vx, g_inject_vy, g_inject_vz}; // injection initial velocity
                        a.a = {0,0,0};
                        a.r = r;
                        // mass from per-type density if specified, else global
                        {
                            real rho = g_density;
                            if (t_id>=0 && t_id<4 && g_density_type[t_id] > 0.0) {
                                rho = g_density_type[t_id];
                            }
                            a.m = rho * (4.0/3.0) * real(M_PI) * r*r*r;
                        }
                        a.type_id = (uint8_t)(t_id & 0xFF);
                        a.asleep = false;
                        a.calm_steps = 0;
                        a.prev_p = a.p;
                        P[oldN+k] = a;
                    }
                } // omp parallel
#else
                // Serial fallback when OpenMP is not enabled
                std::mt19937 trng((unsigned)g_seed + (unsigned)it*101u);
                for (int k=0;k<want;++k){
                    real d; real r; real x; real y;
                    int  t_id;
                    for(;;){
                        t_id = sample_type(trng);
                        d = sample_diameter_for_type(trng, t_id);
                        r = 0.5*d;
                        sample_injection_xy(trng, r, x, y);
                        if (x==x && y==y) break;
                    }
                    real z = g_L + 2*r;
                    Particle a;
                    a.p = {x,y,z};
                    a.v = {g_inject_vx, g_inject_vy, g_inject_vz};
                    a.a = {0,0,0};
                    a.r = r;
                    {
                        real rho = g_density;
                        if (t_id>=0 && t_id<4 && g_density_type[t_id] > 0.0) rho = g_density_type[t_id];
                        a.m = rho * (4.0/3.0) * real(M_PI) * r*r*r;
                    }
                    a.type_id = (uint8_t)(t_id & 0xFF);
                    a.asleep = false;
                    a.calm_steps = 0;
                    a.prev_p = a.p;
                    P[oldN+k] = a;
                }
#endif
                injected_total += want;
                // update running volume cheaply for new range
                for (int k = oldN; k < oldN + want; ++k) {
                    const real r = P[k].r;
                    totalVol_running += (4.0/3.0) * real(M_PI) * r*r*r;
                }
            }
        }
        {
            // Gentle local relaxation for recently injected particles (avoid float suffixes)
            real max_diam = 1e-6 * max_diam_from_types_um();
            max_diam += 2.0 * g_cushion;
            // account for effective collision diameter inflation from cushion
            real z_thresh = g_L + 0.5 * max_diam; // consider particles placed above top plane
            std::vector<int> recent;
            recent.reserve(64);
            for (int i = 0; i < (int)P.size(); ++i) {
                if (P[i].p.z > z_thresh) recent.push_back(i);
            }

            if (!recent.empty()) {
                const int relax_iters = 3;
                for (int iter_rel = 0; iter_rel < relax_iters; ++iter_rel) {
                    for (int ii = 0; ii < (int)recent.size(); ++ii) {
                        int i = recent[ii];
                        // ensure walls/piston constraints satisfied for i
                        const real no_top = g_L + 1000.0 * max_diam;
                        collide_walls(P[i], no_top, g_e_pw);
                    }
                }
            }
        }

        // 2) rebuild neighbor grid
        rebuild_grid(P, head, next);

        // ------------------------------------------------------------------
        // (A) — wake logic BEFORE collisions so sleepers can participate
        // ------------------------------------------------------------------
        // Compute current top plane position (reuse later)
        real zTop_now = topPlane(t);
        g_piston_wall_vn = ((t >= g_ram_t0) && (t <= (g_ram_t0 + g_ram_dt)) ? g_ram_speed : 0.0);
        {
            // Wake everyone during the ram window (simple, robust)
            if (t >= g_ram_t0 && t <= (g_ram_t0 + g_ram_dt)) {
                #pragma omp parallel for
                for (int i=0; i<(int)P.size(); ++i) {
                    P[i].asleep = false;
                    P[i].calm_steps = 0;
                }
            }

            // Band-wake when piston moves downward (local wake near top)
            static real prevTop = -1;
            bool piston_down = false;
            if (prevTop < 0) prevTop = zTop_now;
            if (zTop_now < prevTop - 1e-15) piston_down = true;
            prevTop = zTop_now;
            if (piston_down) {
                #pragma omp parallel for
                for (int i=0; i<(int)P.size(); ++i) {
                    if (P[i].asleep && P[i].p.z > (zTop_now - g_wake_band * P[i].r)) {
                        P[i].asleep = false; P[i].calm_steps = 0;
                    }
                }
            }
        }

        // 4) pairwise collisions via neighbor search (particle-driven, not cell-scan)
        #pragma omp parallel for schedule(static)
        for (int i=0; i<(int)P.size(); ++i){
            // compute i's cell
            int ix = clampi(int(std::floor((P[i].p.x - g_minB.x)/g_cell_h)), 0, gx-1);
            int iy = clampi(int(std::floor((P[i].p.y - g_minB.y)/g_cell_h)), 0, gy-1);
            int iz = clampi(int(std::floor((P[i].p.z - g_minB.z)/g_cell_h)), 0, gz-1);

            for (int dz=-1; dz<=1; ++dz){
                int jz = clampi(iz+dz,0,gz-1);
                for (int dy=-1; dy<=1; ++dy){
                    int jy = clampi(iy+dy,0,gy-1);
                    for (int dx=-1; dx<=1; ++dx){
                        int jx = clampi(ix+dx,0,gx-1);
                        int h = idx3(jx,jy,jz);
                        for (int j=head[h]; j!=-1; j=next[j]){
                            if (j<=i) continue;
                            if (P[i].asleep && P[j].asleep) continue;
                            collide_pair(P[i], P[j], g_e_pp, g_tangent_damp);
                        }
                    }
                }
            }
        }

        // 4.9) set accelerations for this step: gravity + shaking-frame accel
        real a_shake_x=0.0, a_shake_y=0.0, a_shake_z=0.0;
        shake_accels(t, a_shake_x, a_shake_y, a_shake_z);
#if defined(_OPENMP)
        #pragma omp parallel for
        for (int i=0; i<(int)P.size(); ++i){
            auto &a = P[i];
            // reset
            a.a.x = 0.0; a.a.y = 0.0; a.a.z = 0.0;
            // do NOT push on sleepers
            if (!a.asleep) {
                // gravity (+ optional shaking frame accel)
                a.a.x += a_shake_x;
                a.a.y += a_shake_y;
                a.a.z -= (g_g_run + a_shake_z);
            }
        }
#else
        struct AccelCtx { std::vector<Particle>* P; real a_shake_x; real a_shake_y; real a_shake_z; };
        auto accel_worker = [](int i, void* vctx){
            AccelCtx* c = (AccelCtx*)vctx;
            auto &a = (*(c->P))[i];
            a.a.x = 0.0; a.a.y = 0.0; a.a.z = 0.0;
            if (!a.asleep) {
                a.a.x += c->a_shake_x;
                a.a.y += c->a_shake_y;
                a.a.z -= (g_g_run + c->a_shake_z);
            }
        };
        AccelCtx ctxA{&P, a_shake_x, a_shake_y, a_shake_z};
        for (int i = 0; i < (int)P.size(); ++i) accel_worker(i, &ctxA);
#endif

        // 4.95) optional short-range repulsive accelerations (pp + walls)
        // Note: uses the already-built neighbor grid.
        apply_repulsion_accel(P, head, next, zTop_now);

        // 5) integrate (Velocity Verlet: here simple symplectic Euler style for brevity)
        real zTop = zTop_now;
        const real v_damp_fac = (g_lin_damp > 0.0 ? std::exp(-g_lin_damp * g_dt) : 1.0);

#if defined(_OPENMP)
        #pragma omp parallel for
        for (int i=0;i<(int)P.size();++i){
            auto &a = P[i];
            if (!a.asleep) {
                // v_{n+1/2}
                a.v.x += a.a.x * g_dt;
                a.v.y += a.a.y * g_dt;
                a.v.z += a.a.z * g_dt;
                a.v.x *= v_damp_fac;
                a.v.y *= v_damp_fac;
                a.v.z *= v_damp_fac;
                // x_{n+1}
                a.p.x += a.v.x * g_dt;
                a.p.y += a.v.y * g_dt;
                a.p.z += a.v.z * g_dt;
            }

            // Always enforce walls (incl. bottom) to prevent “sleep drift”
            collide_walls(a, zTop, g_e_pw);

            // Extra floor safety clamp (belt-and-suspenders)
            const real zbot_guard = a.r + g_cushion;
            if (a.p.z < zbot_guard) {
                a.p.z = zbot_guard;
                if (a.v.z < 0.0) a.v.z = 0.0;
            }

            // Sleep logic
            bool slow = (std::abs(a.v.z) < g_sleep_v) &&
                        (a.v.x*a.v.x + a.v.y*a.v.y < g_sleep_v*g_sleep_v);

            Vec3 dp = {a.p.x - a.prev_p.x, a.p.y - a.prev_p.y, a.p.z - a.prev_p.z};
            bool on_floor = (a.p.z <= ((a.r + g_cushion) + 1e-12));
            bool not_moving = (std::abs(dp.x) < g_sleep_dz &&
                               std::abs(dp.y) < g_sleep_dz &&
                               std::abs(dp.z) < g_sleep_dz);

            if (!a.asleep) {
                if (on_floor && slow && not_moving) {
                    a.calm_steps++;
                    if (a.calm_steps >= g_sleep_N) {
                        a.asleep = true;
                        a.v = {0,0,0};  // lock it
                    }
                } else {
                    a.calm_steps = 0;
                }
            }

            a.prev_p = a.p; // update history every step
        }
#else
        struct IntegrateCtx { std::vector<Particle>* P; real zTop; };
        auto integ_worker = [](int i, void* vctx){
            IntegrateCtx* c = (IntegrateCtx*)vctx;
            auto &a = (*(c->P))[i];

            if (!a.asleep) {
                a.v.x += a.a.x * g_dt;
                a.v.y += a.a.y * g_dt;
                a.v.z += a.a.z * g_dt;
                a.v.x *= v_damp_fac;
                a.v.y *= v_damp_fac;
                a.v.z *= v_damp_fac;
                a.p.x += a.v.x * g_dt;
                a.p.y += a.v.y * g_dt;
                a.p.z += a.v.z * g_dt;
            }

            collide_walls(a, c->zTop, g_e_pw);

            const real zbot_guard = a.r + g_cushion;
            if (a.p.z < zbot_guard) {
                a.p.z = zbot_guard;
                if (a.v.z < 0.0) a.v.z = 0.0;
            }

            bool slow = (std::abs(a.v.z) < g_sleep_v) &&
                        (a.v.x*a.v.x + a.v.y*a.v.y < g_sleep_v*g_sleep_v);

            Vec3 dp = {a.p.x - a.prev_p.x, a.p.y - a.prev_p.y, a.p.z - a.prev_p.z};
            bool on_floor = (a.p.z <= ((a.r + g_cushion) + 1e-12));
            bool not_moving = (std::abs(dp.x) < g_sleep_dz &&
                               std::abs(dp.y) < g_sleep_dz &&
                               std::abs(dp.z) < g_sleep_dz);

            if (!a.asleep) {
                if (on_floor && slow && not_moving) {
                    a.calm_steps++;
                    if (a.calm_steps >= g_sleep_N) {
                        a.asleep = true;
                        a.v = {0,0,0};
                    }
                } else {
                    a.calm_steps = 0;
                }
            }

            a.prev_p = a.p;
        };
        IntegrateCtx ctxI{&P, zTop};
        for (int i = 0; i < (int)P.size(); ++i) integ_worker(i, &ctxI);
#endif


        // 5.5) successive damping sweeps (fast settling)
        if (g_phi_reached) {
            successive_damping(P, zTop);
        }

        // 6) optional early stop: once phi_target reached and the ensemble settles
        {
            static int ok_checks = 0;

            const bool want_stop = g_phi_reached && (
                (g_stop_vrms > 0.0) || (g_stop_vmax > 0.0) || (g_stop_sleep_frac > 0.0));

            if (want_stop && g_stop_check_interval > 0 && (it % g_stop_check_interval) == 0 && !P.empty()) {
                real sumv2 = 0.0;
                real vmax = 0.0;
                int asleep_count = 0;

                #pragma omp parallel for reduction(+:sumv2,asleep_count) reduction(max:vmax)
                for (int i=0; i<(int)P.size(); ++i) {
                    const Particle& a = P[i];
                    if (a.asleep) asleep_count++;
                    const real v2 = dot(a.v,a.v);
                    sumv2 += v2;
                    const real v = std::sqrt(v2);
                    if (v > vmax) vmax = v;
                }

                const real vrms = std::sqrt(sumv2 / (real)P.size());
                const real sleep_frac = (real)asleep_count / (real)P.size();

                bool ok = true;
                if (g_stop_vrms > 0.0) ok = ok && (vrms <= g_stop_vrms);
                if (g_stop_vmax > 0.0) ok = ok && (vmax <= g_stop_vmax);
                if (g_stop_sleep_frac > 0.0) ok = ok && (sleep_frac >= g_stop_sleep_frac);

                ok_checks = ok ? (ok_checks + 1) : 0;

                if (g_debug) {
                    std::printf("[stop-check] it=%d t=%.3g vrms=%.3g vmax=%.3g asleep=%.3f ok=%d (%d/%d)\n",
                                it, t, vrms, vmax, sleep_frac, ok?1:0, ok_checks, g_stop_checks_required);
                }

                if (ok_checks >= std::max(1, g_stop_checks_required)) {
                    std::printf("[early-stop] settled after phi_target: it=%d t=%.6g vrms=%.3g vmax=%.3g asleep=%.3f\n",
                                it, t, vrms, vmax, sleep_frac);
                    break;
                }
            }
        }

        // 7) dumps & monitor
        if ( (g_dump_interval>0 && (it % g_dump_interval) == 0) ||
             (g_xyz_interval>0 && (it % g_xyz_interval) == 0) ||
             (g_vtk_interval>0 && (it % g_vtk_interval) == 0) ||
             (g_vtk_domain_interval>0 && (it % g_vtk_domain_interval) == 0) ){
            if (!P.empty()){
                if (g_xyz_interval > 0 && (it % g_xyz_interval)==0) dump_xyz(P, it);
                if (g_vtk_interval > 0 && (it % g_vtk_interval)==0) dump_vtk(P, it);
                if (g_vtk_domain_interval > 0 && (it % g_vtk_domain_interval)==0) dump_domain_vtk(zTop_now, it);
            }
            real E=0.0, sumz=0.0;
            #pragma omp parallel for reduction(+:E,sumz)
            for (int i=0;i<(int)P.size();++i){
                E += 0.5*P[i].m*dot(P[i].v,P[i].v);
                sumz += P[i].p.z;
            }
            // phi from running accumulator (cheap)
            phi = (domVol > 0.0 ? totalVol_running / domVol : 0.0);
            std::printf("it=%d t=%.3g N=%zu topZ=%.4g avgZ=%.4g phi=%6.3f\n",
                it, t, P.size(), zTop, P.empty()?0.0:real(sumz/P.size()),phi);
        }
    }

    std::printf("Done. Injected=%d, final N=%zu, phi=%6.3f\n", injected_total, P.size(), phi);
    return 0;
}
