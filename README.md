# Coaxial Fill DEM (CPU)

3D hard-sphere (DEM-style) packing in a **coaxial/annular column** with gravity-driven filling, optional vertical shaking, and a top piston/ram stage.  
The simulation uses a cell-linked neighbor search, impulse-based collisions, and energy damping for realistic settling behavior.

---

## 🧩 Features

- Annular domain: inner/outer cylindrical walls with finite height `L`; gravity acts in −z.  
- Time-resolved injection by particle **flux** until a desired **fill_time**.  
- Multi-type particle size distributions with optional per-type densities.  
- Contact model: normal restitution + tangential friction (static/kinetic) with damping.  
- Walls: inner & outer cylinders, bottom plane (stick/slip), and optional top piston (ram).  
- Optional **vertical shaking** via sinusoidal acceleration.  
- Outputs `.xyz` and `.vtk` files for visualization in ParaView.  
- CLI prints energy, mean height, and packing fraction φ during runtime.

---

## ⚙️ Build Instructions

### Linux / macOS (Clang or GCC)
```bash
clang++ -O3 -march=native -fopenmp -std=c++17 -o coax_pack_cpu coax_pack_cpu.cpp
# or
g++ -O3 -march=native -fopenmp -std=c++17 -o coax_pack_cpu coax_pack_cpu.cpp
```

### Windows (MSVC)
```bat
cl /O2 /openmp /std:c++17 coax_pack_cpu.cpp /Fe:coax_pack_cpu.exe
```

> If needed, add `-lm` on Linux to link the math library.

---

## ▶️ Run Example

```bash
./coax_pack_cpu   20000      1e-6  200000 1000   0   12345   0.00010    0.00030 0.00038   2000       9.81   0.0    0.0    0.50   1.0        0.25   0.0 0.5
```

This example injects up to **20,000** spheres at **2,000 s⁻¹** into an annulus  
(Rin = 100 µm, Rout = 300 µm, L = 380 µm) under gravity, with no shaking or ram compression.

---

## 🧮 Command-Line Parameters

| # | Parameter | Description |
|:-:|:-----------|:-------------|
| 1 | `natoms_max` | Total number of particles to inject |
| 2 | `dt` | Time step (s) |
| 3 | `niter` | Number of integration steps |
| 4 | `dump_interval` | Write output every N steps |
| 5 | `debug` | Debug flag (0/1) |
| 6 | `seed` | RNG seed (≤0 = time-based) |
| 7 | `Rin` | Inner cylinder radius (m) |
| 8 | `Rout` | Outer cylinder radius (m) |
| 9 | `L` | Column height (m) |
| 10 | `flux` | Injection rate (particles/s) |
| 11 | `g` | Gravity (m/s²) |
| 12 | `shake_hz` | Shake frequency (Hz, 0 = off) |
| 13 | `shake_amp` | Shake amplitude (m) |
| 14 | `fill_time` | Duration of injection (s) |
| 15 | `ram_start` | Ram start time (s) |
| 16 | `ram_duration` | Ram duration (s) |
| 17 | `ram_speed` | Downward piston speed (m/s) |
| 18 | `VF` | Volume Fraction Solid (0-1) |

---

## 🧠 Physics Overview

- **Integration:** Semi-implicit Euler scheme (`v += a*dt; p += v*dt`).  
- **Contacts:**  
  - Normal restitution and frictional damping.  
  - Tangential friction (static/kinetic) with small-slip “snap.”  
  - Position correction capped to avoid particle teleportation.  
- **Walls:**  
  - Cylindrical inner/outer boundaries.  
  - Bottom plane with adjustable sticking/sliding.  
  - Optional moving ram at the top boundary.  
- **Neighbor Search:**  
  Uniform 3D cell-linked list with 27-cell checks per particle.  
- **Sleeping/Waking:**  
  Slow, settled particles near the bottom are marked asleep and stop integrating.  
  Collisions or ram motion can wake them.

---

## ⚡ Random Vertical Shaking

The container’s vertical acceleration is controlled by:
```cpp
a_z = g_shake_amp * sin(2π * g_shake_hz * t);
```
You can randomize it or modulate it for better packing efficiency.

---

## 📊 Output Files

| File | Description |
|:------|:-------------|
| `atoms_<iter>.xyz` | Text file with `element x y z r` |
| `atom.<iter>.vtk` | ParaView-compatible format with velocity vectors |
| stdout | Prints iteration, time, number of particles, mean height, and packing fraction φ |

---

## 🧱 Particle Distributions

- Supports **four types** (e.g., La, Sr, Fe, Co) with tabulated cumulative radii distributions in micrometers.  
- Fractions are defined in the `g_type_frac[]` array and normalized at runtime.  
- Densities can be customized per type for mass-weighted dynamics.

---

## 🧪 Performance Tips

- Compile with `-O3 -march=native -fopenmp`.  
- Increase `dump_interval` to reduce I/O load.  
- Use a small `dt` for high-speed impacts or small particles.  
- Enable *successive damping* for faster settling.  
- Fix RNG seed for reproducible packing.

---

## 🎨 Visualization

Open `.vtk` outputs in **ParaView**:  
- Color by `radius` or `type_id`.  
- Add velocity glyphs to view flow during injection or shaking.

For lightweight inspection, `.xyz` files can plotted using the included python script.

---

## 📜 License
- Free as in beer
  
---


