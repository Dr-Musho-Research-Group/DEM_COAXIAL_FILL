#!/usr/bin/env python3
"""
Render four views from atoms_*.xyz frames written by coax_pack_cpu:
  1) Top view (x-y)     -> atoms_xy.mp4
  2) X-view  (y-z)      -> atoms_yz.mp4
  3) Y-view  (x-z)      -> atoms_xz.mp4
  4) Isometric (3D)     -> atoms_iso.mp4

Batch auto-parse (optional): if _run_win_cpu.bat exists, read set RIN, ROUT, LENGTH, DT.

XYZ row format supported (per particle):
   'SYMBOL x y z d'  (preferred; SYMBOL e.g., La/Sr/Fe/Co)
   'x y z d'         (fallback; symbol -> 'X')

line 1: N
line 2: "iter <int>"
"""
import os, re, glob, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle, Patch
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# -------------------- Helper: parse batch variables --------------------
def parse_batch_vars(filename):
    """Parse numeric 'set VAR=VALUE' assignments from a Windows .bat file.
    Accepts both:   set VAR=123   and   set "VAR=123"
    Silently ignores non-numeric values (paths, '..', %variables%, etc.).
    """
    vals = {}
    if not os.path.exists(filename):
        return vals
    # Match: set VAR=VALUE   or   set "VAR=VALUE"   (VALUE may be unquoted)
    pat = re.compile(r"^\s*set\s+(?:\"?)([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)(?:\"?)\s*$", re.IGNORECASE)
    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.lower().startswith("rem"):
                continue
            m = pat.match(line)
            if not m:
                continue
            key = m.group(1).strip().upper()
            val_s = m.group(2).strip().strip('"')
            # Skip values that obviously are not plain numbers
            if any(c in val_s for c in ["%", "\\", ":", "(", ")"]):
                continue
            try:
                vals[key] = float(val_s)
            except Exception:
                # e.g., '..' or other non-numeric tokens
                continue
    return vals
batch_vals = parse_batch_vars("_run_win_cpu.bat")

# -------------------- Args & env --------------------
ap = argparse.ArgumentParser()
ap.add_argument("--pattern", default="atoms_*.xyz", help="Glob for frame files (default: atoms_*.xyz)")
ap.add_argument("--rin", type=float, default=batch_vals.get("RIN"))
ap.add_argument("--rout", type=float, default=batch_vals.get("ROUT"))
ap.add_argument("--length", type=float, default=batch_vals.get("LENGTH"))
ap.add_argument("--dt", type=float, default=batch_vals.get("DT"))
ap.add_argument("--fps", type=int, default=1, help="Frames per second")
ap.add_argument("--ms", type=float, default=16.0,
               help="Marker area (points^2) for a *median* diameter particle")
ap.add_argument("--limit", type=int, default=None, help="Limit number of frames")
# Back-compat: --cmap=z behaves like --color=z; otherwise ignored
ap.add_argument("--cmap", choices=["none", "z"], default="none",
               help=argparse.SUPPRESS)
ap.add_argument("--color", choices=["element", "z", "none"], default="element",
               help="Coloring mode: element (La/Sr/Fe/Co), z (height), or none (single color)")
ap.add_argument("--dpi", type=int, default=150)
args = ap.parse_args()

# Back-compat shim
color_mode = args.color
if args.cmap == "z" and args.color == "element":
    color_mode = "z"

RIN, ROUT, LENGTH, DT = args.rin, args.rout, args.length, args.dt
FPS, MS, LIMIT, DPI = args.fps, args.ms, args.limit, args.dpi

print(f"\nEffective parameters:")
print(f"  RIN={RIN}, ROUT={ROUT}, LENGTH={LENGTH}, DT={DT}")
print(f"  Color mode: {color_mode}\n")

# -------------------- Discover frames --------------------
xyz_files = sorted(
    glob.glob(args.pattern),
    key=lambda x: int(re.search(r"atoms_(\d+)\.xyz$", x).group(1)) if re.search(r"atoms_(\d+)\.xyz$", x) else 0
)
if not xyz_files:
    print(f"No files found for pattern: {args.pattern}")
    raise SystemExit(1)
if LIMIT:
    xyz_files = xyz_files[:LIMIT]
print(f"Found {len(xyz_files)} frames.")

# -------------------- Loading --------------------
def read_xyz_positions_symbols(fn):
    """Return (iter_idx, pos[N,3], diam[N], symbols[N]) from XYZ.
       Supported per-particle row formats (extra columns are ignored):
         1) SYMBOL x y z d ...        (SYMBOL like La/Sr/Fe/Co)
         2) type_id x y z d ...       (type_id numeric; maps 0..3 -> La/Sr/Fe/Co)
         3) x y z d ...               (fallback; symbol='X')
    """
    type_map = {0: 'La', 1: 'Sr', 2: 'Fe', 3: 'Co'}
    with open(fn, 'r', encoding='utf-8', errors='ignore') as f:
        first = f.readline()
        if not first:
            raise ValueError(f'Empty file: {fn}')
        N = int(first.strip())
        hdr = f.readline().strip()
        m = re.search(r'iter\s+(\d+)', hdr, flags=re.IGNORECASE)
        it = int(m.group(1)) if m else 0

        pos = np.zeros((N, 3), dtype=float)
        diam = np.zeros(N, dtype=float)
        symbols = np.empty(N, dtype=object)

        j = 0
        for line in f:
            if j >= N:
                break
            if not line.strip():
                continue
            toks = line.split()
            if not toks:
                continue

            if not is_float(toks[0]):
                # SYMBOL x y z d ...
                sym = toks[0]
                x, y, z, d = map(float, toks[1:5])
            else:
                # Numeric first token: could be type_id x y z d ... OR x y z d ...
                if len(toks) >= 5 and is_int_like(toks[0]) and all(is_float(t) for t in toks[1:5]):
                    tid = int(float(toks[0]))
                    sym = type_map.get(tid, str(tid))
                    x, y, z, d = map(float, toks[1:5])
                else:
                    sym = 'X'
                    x, y, z, d = map(float, toks[0:4])

            pos[j, :] = (x, y, z)
            diam[j] = d
            symbols[j] = sym
            j += 1

        if j != N:
            # Tolerate short reads but resize to what we actually got
            pos = pos[:j, :]
            diam = diam[:j]
            symbols = symbols[:j]

    return it, pos, diam, symbols

def is_float(s):
    try:
        float(s)
        return True
    except Exception:
        return False


def is_int_like(s):
    """True if token looks like an integer (no decimal point/exponent)."""
    s = s.strip()
    if not s:
        return False
    if any(c in s.lower() for c in ['.', 'e']):
        return False
    try:
        int(s)
        return True
    except Exception:
        return False
iters, frames_pos, frames_d, frames_sym = [], [], [], []
print("Scanning frames for bounds and diameters...")
for fn in xyz_files:
    it, pos, diam, sym = read_xyz_positions_symbols(fn)
    iters.append(it)
    frames_pos.append(pos)
    frames_d.append(diam)
    frames_sym.append(sym)

iters = np.array(iters, dtype=int)

# -------------------- Size normalization (median-diameter anchored) --------------------
all_diam = np.concatenate(frames_d) if frames_d else np.array([1.0])
D_med = float(np.median(all_diam)) if all_diam.size else 1.0
if D_med <= 0:
    D_med = 1.0
print(f"Median diameter across frames: {D_med:.6g} m")
def sizes_for_frame(dvec):
    # marker area ∝ d^2; median d -> area = MS (points^2)
    return ((dvec / D_med) ** 2) * MS

# -------------------- Element → color mapping --------------------
ELEMENT_COLORS = {
    'La': 'tab:blue',
    'Sr': 'tab:orange',
    'Fe': 'tab:red',
    'Co': 'tab:purple',
    'X':  'tab:gray'  # fallback/unknown
}
def colors_for_symbols(symbols):
    # Map each symbol to a matplotlib color string
    return np.array([ELEMENT_COLORS.get(s, ELEMENT_COLORS['X']) for s in symbols], dtype=object)

# -------------------- Bounds --------------------
all_xyz = np.vstack(frames_pos)
rxy = np.sqrt(all_xyz[:, 0]**2 + all_xyz[:, 1]**2)
data_rmax = float(np.max(rxy)) if rxy.size else 1.0
zmin_d, zmax_d = float(np.min(all_xyz[:, 2])), float(np.max(all_xyz[:, 2]))

ROUT_plot = (ROUT if ROUT else data_rmax) * 1.05
RIN_plot = RIN if RIN else 0.0
zmin_plot, zmax_plot = 0.0, (LENGTH if LENGTH else zmax_d)

# -------------------- Utilities --------------------
def try_save_mp4(ani, outfile):
    try:
        ani.save(outfile, writer='ffmpeg', fps=FPS, dpi=DPI)
        print(f"Saved MP4: {outfile}")
        return True
    except Exception as e:
        print(f"FFMPEG not available or failed for {outfile} ({e}).")
        return False

def label_text(iter_val, npts):
    return f"iter: {iter_val}   t={iter_val*DT:.4f}s   N={npts}" if DT else f"iter: {iter_val}   N={npts}"

def add_element_legend(ax):
    handles = [Patch(facecolor=ELEMENT_COLORS[k], edgecolor='none', label=k) for k in ['La','Sr','Fe','Co']]
    ax.legend(handles=handles, loc='upper right', framealpha=0.7)

# ===================================================
# 1) Top View (XY)
# ===================================================
fig_xy, ax_xy = plt.subplots(figsize=(8, 8))
ax_xy.set_aspect('equal', adjustable='box')
ax_xy.set_xlim(-ROUT_plot, ROUT_plot)
ax_xy.set_ylim(-ROUT_plot, ROUT_plot)
ax_xy.set_xlabel('x [m]')
ax_xy.set_ylabel('y [m]')
ax_xy.set_title('Top view (XY)')

# Inner + outer coax walls
annulus_art_xy = []
if RIN_plot > 0:
    inner = Circle((0, 0), RIN_plot, fill=False, lw=1.5, ls='--', color='gray')
    ax_xy.add_patch(inner); annulus_art_xy.append(inner)
outer = Circle((0, 0), ROUT_plot / 1.05, fill=False, lw=1.8, ls='-', color='black')
ax_xy.add_patch(outer); annulus_art_xy.append(outer)

if color_mode == "z":
    sc_xy = ax_xy.scatter([], [], s=[], cmap='viridis', vmin=zmin_plot, vmax=zmax_plot)
else:
    sc_xy = ax_xy.scatter([], [], s=[], alpha=0.85)

txt_xy = ax_xy.text(0.02, 0.98, "", transform=ax_xy.transAxes,
                    ha='left', va='top', fontsize=10,
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
if color_mode == "element":
    add_element_legend(ax_xy)

def update_xy(i):
    pos = frames_pos[i]
    diam = frames_d[i]
    sc_xy.set_offsets(pos[:, :2])
    sc_xy.set_sizes(sizes_for_frame(diam))
    if color_mode == "z":
        sc_xy.set_array(pos[:, 2])
    elif color_mode == "element":
        sc_xy.set_color(colors_for_symbols(frames_sym[i]))
    txt_xy.set_text(label_text(iters[i], pos.shape[0]))
    ax_xy.set_title(f"Top view (XY) frame {i+1}/{len(frames_pos)}")
    return sc_xy, txt_xy, *annulus_art_xy

ani_xy = FuncAnimation(fig_xy, update_xy, frames=len(frames_pos), interval=1000.0 / FPS, blit=False)
try_save_mp4(ani_xy, "atoms_xy.mp4")
fig_xy.savefig("atoms_xy_last.png", dpi=200)

# ===================================================
# 2) Side Views (YZ, XZ)
# ===================================================
def side_view(coord1, coord2, name, xlabel, ylabel):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xlim(-ROUT_plot, ROUT_plot)
    ax.set_ylim(zmin_plot, zmax_plot)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{name} view")

    # Inner/outer lines
    if RIN_plot > 0:
        ax.axvline(+RIN_plot, ls='--', lw=1.0, color='gray')
        ax.axvline(-RIN_plot, ls='--', lw=1.0, color='gray')
    ax.axvline(+ROUT_plot / 1.05, ls='-', lw=1.2, color='black')
    ax.axvline(-ROUT_plot / 1.05, ls='-', lw=1.2, color='black')

    if color_mode == "z":
        sc = ax.scatter([], [], s=[], cmap='viridis', vmin=zmin_plot, vmax=zmax_plot)
    else:
        sc = ax.scatter([], [], s=[], alpha=0.85)
    txt = ax.text(0.02, 0.98, "", transform=ax.transAxes, ha='left', va='top',
                  fontsize=10, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    if color_mode == "element":
        add_element_legend(ax)

    def update(i):
        pos = frames_pos[i]
        diam = frames_d[i]
        p1, p2 = pos[:, coord1], pos[:, coord2]
        sc.set_offsets(np.c_[p1, p2])
        sc.set_sizes(sizes_for_frame(diam))
        if color_mode == "z":
            sc.set_array(pos[:, 2])
        elif color_mode == "element":
            sc.set_color(colors_for_symbols(frames_sym[i]))
        txt.set_text(label_text(iters[i], len(p1)))
        ax.set_title(f"{name} view frame {i+1}/{len(frames_pos)}")
        return sc, txt

    ani = FuncAnimation(fig, update, frames=len(frames_pos), interval=1000.0 / FPS, blit=False)
    outbase = name.lower().replace(' ', '_')
    try_save_mp4(ani, f"atoms_{outbase}.mp4")
    fig.savefig(f"atoms_{outbase}_last.png", dpi=200)

side_view(1, 2, "YZ (X-view)", "y [m]", "z [m]")
side_view(0, 2, "XZ (Y-view)", "x [m]", "z [m]")

# ===================================================
# 3) Isometric 3D
# ===================================================
fig_iso = plt.figure(figsize=(8, 7))
ax_iso = fig_iso.add_subplot(111, projection='3d')
ax_iso.set_xlim(-ROUT_plot, ROUT_plot)
ax_iso.set_ylim(-ROUT_plot, ROUT_plot)
ax_iso.set_zlim(zmin_plot, zmax_plot)
ax_iso.set_xlabel('x [m]'); ax_iso.set_ylabel('y [m]'); ax_iso.set_zlabel('z [m]')
ax_iso.set_title('Isometric 3D')

x0, y0, z0 = frames_pos[0][:, 0], frames_pos[0][:, 1], frames_pos[0][:, 2]
s0 = sizes_for_frame(frames_d[0])
if color_mode == "z":
    scat = ax_iso.scatter(x0, y0, z0, s=s0, c=z0, cmap='viridis', vmin=zmin_plot, vmax=zmax_plot, depthshade=False)
else:
    # initialize with element colors or single color
    c0 = colors_for_symbols(frames_sym[0]) if color_mode == "element" else None
    scat = ax_iso.scatter(x0, y0, z0, s=s0, c=c0, alpha=0.85, depthshade=False)
txt_iso = ax_iso.text2D(0.02, 0.98, label_text(iters[0], len(x0)), transform=ax_iso.transAxes,
                        fontsize=10, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
if color_mode == "element":
    add_element_legend(ax_iso)

def update_iso(i):
    pos = frames_pos[i]
    x, y, z = pos[:,0], pos[:,1], pos[:,2]
    scat._offsets3d = (x, y, z)
    scat.set_sizes(sizes_for_frame(frames_d[i]))
    if color_mode == "z":
        scat.set_array(z)
    elif color_mode == "element":
        colors = colors_for_symbols(frames_sym[i])
        # set_color works for 3D collections in recent matplotlib
        try:
            scat.set_color(colors)
        except Exception:
            # fallback for older mpl
            scat._facecolor3d = colors
            scat._edgecolor3d = colors
    txt_iso.set_text(label_text(iters[i], len(x)))
    ax_iso.view_init(elev=20, azim=(30 + 0.5*i) % 360)
    return scat, txt_iso

ani_iso = FuncAnimation(fig_iso, update_iso, frames=len(frames_pos), interval=1000.0 / FPS, blit=False)
try_save_mp4(ani_iso, "atoms_iso.mp4")
fig_iso.savefig("atoms_iso_last.png", dpi=200)

print("\nAll views complete with diameter-scaled markers:")
print("  atoms_xy.mp4, atoms_yz.mp4, atoms_xz.mp4, atoms_iso.mp4")
print("  (plus *_last.png snapshots)")
