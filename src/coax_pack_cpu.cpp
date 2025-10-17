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
 *     dt                 (real) time step [s]
 *     niter              (int)   number of steps
 *     dump_interval      (int)   dump every N steps (XYZ/VTK)
 *     debug              (int)   0/1
 *     seed               (int)   RNG seed (<=0 -> time-based)
 *     Rin                (real) inner radius  [m]
 *     Rout               (real) outer radius  [m]
 *     L                  (real) length (z)    [m]
 *     flux               (int)   spheres per second injected while filling
 *     g                  (real) gravity magnitude [m/s^2] (points -z)
 *     shake_hz           (real) shake frequency [Hz]; 0 => off
 *     shake_amp          (real) shake amplitude [m] (vertical)
 *     fill_time          (real) seconds to inject (afterwards stop injecting)
 *     ram_start_t        (real) when to start ram (s from t=0)
 *     ram_duration       (real) duration of ram [s]
 *     ram_speed          (real) downward piston speed [m/s] (top plate)
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
#include <omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

int   g_natoms_max;
real g_dt;
int   g_niter;
int   g_dump_interval;

real g_Rin, g_Rout, g_L;
int   g_flux;                  // spheres/sec
real g_g = 9.81f;             // m/s^2 downward

// shaking (container vertical acceleration via frame motion)
real g_shake_hz = 0.0f;
real g_shake_amp = 0.0f;      // meters

// injection
real g_fill_time = 0.0f;      // seconds to inject

// ram (top piston)
real g_ram_t0 = 0.0f;
real g_ram_dt = 0.0f;
real g_ram_speed = 0.0f;      // downward (positive)

// materials / contacts
real g_density = 2000.0f;     // kg/m^3 (default)
real g_e_pp = 0.015f;          // restitution particle-particle
real g_e_pw = 0.015f;          // restitution particle-wall
real g_tangent_damp = 0.015f;  // simple tangential vel damping

// neighbor grid
real g_cell_h;                // cell size ~ max_diam
int   gx, gy, gz;
Vec3  g_minB, g_maxB;

// distribution (diameters) — mirrors your cumulative sampling style
static std::vector<real> diam_list_um = {
    1.25,1.46,1.67,1.88,2.09,2.30,2.51,2.72,2.93,3.14,
    3.35,3.56,3.76,3.97,4.18,4.39,4.60,4.81,5.02,5.23
};
static std::vector<real> cumu = {
    0.0250f, 0.0328f, 0.1197f, 0.3573f, 0.6326f,
    0.7687f, 0.8431f, 0.8929f, 0.9206f, 0.9420f,
    0.9609f, 0.9666f, 0.9748f, 0.9824f, 0.9868f,
    0.9905f, 0.9943f, 0.9965f, 0.9980f, 1.0000f
};

// --------------------------- Utilities -----------------------------
static inline int clampi(int v,int lo,int hi){ return v<lo?lo:(v>hi?hi:v); }

void dump_xyz(const std::vector<Particle>& P, int iter){
    char fn[64]; std::snprintf(fn,sizeof(fn),"atoms_%d.xyz",iter);
    FILE* f = std::fopen(fn,"w");
    if(!f){ std::perror("xyz"); return; }
    std::fprintf(f,"%zu\n", P.size());
    std::fprintf(f,"iter %d\n", iter);
    for(const auto& p: P) std::fprintf(f,"C %.9g %.9g %.9g %.9g\n", p.p.x, p.p.y, p.p.z, p.r);
    std::fclose(f);
}
void dump_vtk(const std::vector<Particle>& P, int iter){
    char fn[64]; std::snprintf(fn,sizeof(fn),"atom.%d.vtk",iter);
    FILE* f = std::fopen(fn,"w");
    if(!f){ std::perror("vtk"); return; }
    std::fprintf(f,"# vtk DataFile Version 2.0\nCoax packing\nASCII\nDATASET POLYDATA\n");
    std::fprintf(f,"POINTS %zu real\n", P.size());
    for(const auto& a: P) std::fprintf(f,"%g %g %g\n", a.p.x, a.p.y, a.p.z);
    std::fprintf(f,"POINT_DATA %zu\nVECTORS velocity real\n", P.size());
    for(const auto& a: P) std::fprintf(f,"%g %g %g\n", a.v.x, a.v.y, a.v.z);
    std::fclose(f);
}

// sample diameter from cumulative list (meters)
real sample_diameter(std::mt19937& rng){
    std::uniform_real_distribution<real> U(0.0f,1.0f);
    real r = U(rng);
    size_t i=0; for(; i<cumu.size(); ++i) if(r < cumu[i]) break;
    if(i>=diam_list_um.size()) i = diam_list_um.size()-1;
    return 1e-6f * diam_list_um[i];
}

// uniform point in annulus cross-section for injection (x,y), at z ~ L + margin
void sample_injection_xy(std::mt19937& rng, real& x, real& y){
    std::uniform_real_distribution<real> U(0.0f,1.0f);
    // area-uniform in annulus: r^2 uniform between Rin^2 and Rout^2
    real Rin2 = g_Rin*g_Rin, Rout2 = g_Rout*g_Rout;
    real r = std::sqrt( Rin2 + (Rout2-Rin2)*U(rng) );
    real th = 2.0f*real(M_PI)*U(rng);
    x = r*std::cos(th);
    y = r*std::sin(th);
}

// build grid indexing helpers
inline Vec3 worldMin(){ return g_minB; }
inline int  idx3(int ix,int iy,int iz){ return (iz*gy + iy)*gx + ix; }

void rebuild_grid(const std::vector<Particle>& P,
                  std::vector<int>& head,
                  std::vector<int>& next)
{
    std::fill(head.begin(), head.end(), -1);
    #pragma omp parallel for
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

// -------------------------- Collisions -----------------------------
void collide_pair(Particle& A, Particle& B, real e, real tdamp){
    Vec3 d = B.p - A.p;
    real dist = norm(d);
    real R = A.r + B.r;
    if (dist<=0.0f || dist >= R) return;

    // normalize
    Vec3 n = d*(1.0f/(dist>1e-16f?dist:1e-16f));
    // separate overlap equally
    real pen = (R - dist);
    A.p = A.p - n*(0.5f*pen);
    B.p = B.p + n*(0.5f*pen);

    // relative velocity
    Vec3 vrel = B.v - A.v;
    real vn = dot(vrel,n);
    if (vn < 0.0f){
        real invm = 1.0f/A.m;
        real invn = 1.0f/B.m;
        real j = -(1.0f + e)*vn / (invm + invn);
        Vec3 J = n*j;
        A.v = A.v - J*invm;
        B.v = B.v + J*invn;

        // crude tangential damping
        Vec3 vt = vrel - n*vn;
        A.v = A.v + vt*(-tdamp*invm);
        B.v = B.v - vt*(-tdamp*invn);
    }
}

void collide_walls(Particle& a, real top_z, real e){
    // radial walls
    real rxy = std::sqrt(a.p.x*a.p.x + a.p.y*a.p.y);
    // inner
    real rmin = g_Rin + a.r;
    if(rxy < rmin){
        // push out
        real nx = (rxy>1e-16f? a.p.x/rxy : 1.0f);
        real ny = (rxy>1e-16f? a.p.y/rxy : 0.0f);
        real pen = (rmin - rxy);
        a.p.x -= nx*pen; a.p.y -= ny*pen;
        real vn = a.v.x*nx + a.v.y*ny;
        if (vn<0.0f){
            a.v.x -= (1.0f+e)*vn*nx;
            a.v.y -= (1.0f+e)*vn*ny;
            // tangential damping
            a.v.x *= (1.0f - g_tangent_damp*0.5f);
            a.v.y *= (1.0f - g_tangent_damp*0.5f);
        }
    }
    // outer
    real rmax = g_Rout - a.r;
    if(rxy > rmax){
        real nx = (rxy>1e-16f? a.p.x/rxy : 1.0f);
        real ny = (rxy>1e-16f? a.p.y/rxy : 0.0f);
        real pen = (rxy - rmax);
        a.p.x -= nx*pen; a.p.y -= ny*pen;
        real vn = a.v.x*nx + a.v.y*ny;
        if (vn>0.0f){
            a.v.x -= (1.0f+e)*vn*nx;
            a.v.y -= (1.0f+e)*vn*ny;
            a.v.x *= (1.0f - g_tangent_damp*0.5f);
            a.v.y *= (1.0f - g_tangent_damp*0.5f);
        }
    }
    // bottom plane at z=0
    real zmin = a.r;
    if (a.p.z < zmin){
        real pen = (zmin - a.p.z);
        a.p.z += pen;
        a.v.z = 0;
        a.v.x = 0;
        a.v.y = 0;
    }
    // moving top piston plane at z = top_z
    real zmax = top_z - a.r;
    if (a.p.z > zmax){
        real pen = (a.p.z - zmax);
        a.p.z -= pen;
        if (a.v.z > 0.0f) a.v.z = -a.v.z * e; // reflect downward
        a.v.x *= (1.0f - g_tangent_damp);
        a.v.y *= (1.0f - g_tangent_damp);
    }
}

// ----------------------- Integration Loop -------------------------
int main(int argc, char** argv){
    if (argc < 18){
        std::fprintf(stderr,
            "Usage:\n"
            "%s natoms_max dt niter dump_interval debug seed Rin Rout L flux g shake_hz shake_amp fill_time ram_start ram_duration ram_speed\n", argv[0]);
        return 1;
    }
    int iarg=1;
    g_natoms_max   = std::atoi(argv[iarg++]);
    g_dt           = std::atof(argv[iarg++]);
    g_niter        = std::atoi(argv[iarg++]);
    g_dump_interval= std::atoi(argv[iarg++]);
    g_debug        = std::atoi(argv[iarg++]);
    g_seed         = std::atoi(argv[iarg++]);
    g_Rin          = std::atof(argv[iarg++]);
    g_Rout         = std::atof(argv[iarg++]);
    g_L            = std::atof(argv[iarg++]);
    g_flux         = std::atoi(argv[iarg++]);
    g_g            = std::atof(argv[iarg++]);
    g_shake_hz     = std::atof(argv[iarg++]);
    g_shake_amp    = std::atof(argv[iarg++]);
    g_fill_time    = std::atof(argv[iarg++]);
    g_ram_t0       = std::atof(argv[iarg++]);
    g_ram_dt       = std::atof(argv[iarg++]);
    g_ram_speed    = std::atof(argv[iarg++]);

    if (g_seed<=0) g_seed = int(std::time(nullptr));
    std::mt19937 rng(g_seed);

    // domain bbox for neighbor grid
    real max_diam = 1e-6f * *std::max_element(diam_list_um.begin(), diam_list_um.end());
    g_cell_h = std::max(max_diam, 0.5f*max_diam) * 1.05f; // conservative
    g_minB = {-g_Rout - max_diam, -g_Rout - max_diam, -max_diam};
    g_maxB = { g_Rout + max_diam,  g_Rout + max_diam,  g_L + max_diam};
    gx = std::max(4, int(std::ceil((g_maxB.x - g_minB.x)/g_cell_h)));
    gy = std::max(4, int(std::ceil((g_maxB.y - g_minB.y)/g_cell_h)));
    gz = std::max(8, int(std::ceil((g_maxB.z - g_minB.z)/g_cell_h)));

    std::vector<Particle> P;
    P.reserve(std::min(g_natoms_max, 200000));

    // neighbor grid buffers
    std::vector<int> head(gx*gy*gz,-1);
    std::vector<int> next(g_natoms_max, -1);

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
    real t = 0.0f;

    std::printf("Coax pack: Rin=%.4g m Rout=%.4g m L=%.4g m, fill_time=%.2fs, flux=%d 1/s\n",
                g_Rin, g_Rout, g_L, g_fill_time, g_flux);

    for (int it=0; it<=g_niter; ++it, t += g_dt){
        // 1) Inject new spheres while t < fill_time
        if (t < g_fill_time && injected_total < g_natoms_max && g_flux>0){
            // expected total injected by now:
            int should = int(std::floor(t * g_flux));
            int want = std::min(should - injected_total, g_natoms_max - injected_total);
            if (it==0) want = std::max(want,1); // ensure at least one
            if (want > 0){
                int oldN = (int)P.size();
                P.resize(oldN + want);
                #pragma omp parallel for
                for (int k=0;k<want;++k){
                    real d = sample_diameter(rng);
                    real r = 0.5f*d;
                    real x,y; sample_injection_xy(rng,x,y);
                    // place slightly above top plane
                    real z = g_L + 2*r;
                    Particle a;
                    a.p = {x,y,z};
                    a.v = {0,0,0};
                    a.a = {0,0,0};
                    a.r = r;
                    a.m = g_density * (4.0f/3.0f) * real(M_PI) * r*r*r;
                    P[oldN+k] = a;
                }
                injected_total += want;
            }
        }

        // 2) rebuild neighbor grid
        rebuild_grid(P, head, next);

        // 3) reset acceleration to gravity + shaking (vertical frame acc)
        real ashake = 0.0f;
        if (g_shake_hz > 0.0f && g_shake_amp > 0.0f){
            real w = 2.0f*real(M_PI)*g_shake_hz;
            // z_frame = A * sin(w t) → a_frame = -A w^2 * sin(w t)
            ashake = - g_shake_amp * w*w * std::sin(w*t);
        }
        #pragma omp parallel for
        for (int i=0;i<(int)P.size();++i){
            P[i].a = {0.0f, 0.0f, -(g_g) + ashake};
        }

        // 4) pairwise collisions via neighbor search
        #pragma omp parallel for schedule(static)
        for (int iz=0; iz<gz; ++iz){
            for (int iy=0; iy<gy; ++iy){
                for (int ix=0; ix<gx; ++ix){
                    int h = idx3(ix,iy,iz);
                    for (int i=head[h]; i!=-1; i=next[i]){
                        // visit neighbor cells (3x3x3)
                        for (int dz=-1; dz<=1; ++dz){
                            int jz = clampi(iz+dz,0,gz-1);
                            for (int dy=-1; dy<=1; ++dy){
                                int jy = clampi(iy+dy,0,gy-1);
                                for (int dx=-1; dx<=1; ++dx){
                                    int jx = clampi(ix+dx,0,gx-1);
                                    int hh = idx3(jx,jy,jz);
                                    for (int j=head[hh]; j!=-1; j=next[j]){
                                        if (j<=i) continue; // avoid double
                                        collide_pair(P[i], P[j], g_e_pp, g_tangent_damp);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // 5) integrate (Velocity Verlet: here simple symplectic Euler style for brevity)
        real zTop = topPlane(t);
        #pragma omp parallel for
        for (int i=0;i<(int)P.size();++i){
            auto &a = P[i];
            // v_{n+1/2}
            a.v.x += a.a.x * g_dt;
            a.v.y += a.a.y * g_dt;
            a.v.z += a.a.z * g_dt;
            // x_{n+1}
            a.p.x += a.v.x * g_dt;
            a.p.y += a.v.y * g_dt;
            a.p.z += a.v.z * g_dt;

            // walls & piston
            collide_walls(a, zTop, g_e_pw);
        }

        // 6) basic limiter to avoid blowups
        #pragma omp parallel for
        for (int i=0;i<(int)P.size();++i){
            real vmax = 10.0f;
            real s = std::sqrt(dot(P[i].v,P[i].v));
            if (s>vmax){
                real k = vmax/s;
                P[i].v = P[i].v * k;
            }
        }

        // 7) dumps & monitor
        if ( (it % g_dump_interval) == 0 ){
            if (!P.empty()){
                dump_xyz(P, it);
                dump_vtk(P, it);
            }
            double E=0.0, sumz=0.0;
            #pragma omp parallel for reduction(+:E,sumz)
            for (int i=0;i<(int)P.size();++i){
                E += 0.5*P[i].m*dot(P[i].v,P[i].v);
                sumz += P[i].p.z;
            }
            std::printf("it=%d t=%.3f N=%zu topZ=%.4f avgZ=%.4f\n",
                it, t, P.size(), zTop, P.empty()?0.0:real(sumz/P.size()));
        }
    }

    std::printf("Done. Injected=%d, final N=%zu\n", injected_total, P.size());
    return 0;
}
