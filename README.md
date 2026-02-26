
# Coaxial Fill DEM (CPU)

3D hard-sphere DEM-style packing in a coaxial (annular) column with gravity-driven filling, optional multi-axis shaking, cushion radius support, enhanced wall spring contacts, and an optional top ram stage.

The solver uses a cell-linked neighbor search, impulse-based sphere-sphere collisions, optional spring-damper wall contacts, and damping for realistic settling behavior.

---

## Features

- Annular domain: inner and outer cylindrical walls with finite height L.
- Gravity acts in negative z direction.
- Time-resolved particle injection via flux until fill_time.
- Multi-type particle size distributions with per-type densities.
- Normal restitution plus tangential friction model.
- Cushion radius option for soft-contact buffering.
- Effective wall spring-damper model (independent from particle collisions).
- Independent shaking amplitudes in x, y, and z.
- Optional legacy XY shake compatibility mode.
- Configurable injection velocity via CLI.
- Optional top piston (ram) stage.
- Outputs XYZ and VTK files for ParaView.
- Runtime prints energy, mean height, and packing fraction.

---

## Build Instructions

### Linux or macOS (Clang or GCC)

clang++ -O3 -march=native -fopenmp -std=c++17 -o coax_pack_cpu coax_pack_cpu.cpp

or

g++ -O3 -march=native -fopenmp -std=c++17 -o coax_pack_cpu coax_pack_cpu.cpp

If needed on Linux:
-lm

### Windows (MSVC)

cl /O2 /openmp /std:c++17 coax_pack_cpu.cpp /Fe:coax_pack_cpu.exe

---

## Core Positional Command-Line Parameters

1  natoms_max     Total particles to inject
2  dt             Time step (s)
3  niter          Number of integration steps
4  dump_interval  Output interval
5  debug          Debug flag (0 or 1)
6  seed           RNG seed (<=0 time-based)
7  Rin            Inner radius (m)
8  Rout           Outer radius (m)
9  L              Column height (m)
10 flux           Injection rate (1/s)
11 g              Gravity (m/s^2)
12 shake_hz       Shake frequency (Hz)
13 shake_amp      Legacy shake amplitude (m)
14 fill_time      Injection duration (s)
15 ram_start      Ram start time (s)
16 ram_duration   Ram duration (s)
17 ram_speed      Downward ram speed (m/s)
18 VF             Target solid volume fraction

---

## New Optional CLI Flags

### Cushion Radius

--cushion <meters>

Adds cushion thickness to effective collision radius for particle-particle and particle-wall contacts.
Does not change particle mass or physical radius.
Default: 0.0

---

### Effective Wall Spring-Damper Model

--wall_k <N_per_m>
--wall_zeta <0_to_1>
--wall_dvmax <m_per_s>

wall_k: effective wall stiffness
wall_zeta: damping ratio
wall_dvmax: max per-step velocity correction clamp

Set wall_k = 0 to disable.

---

### Injection Velocity Control

--inject_vx <m_per_s>
--inject_vy <m_per_s>
--inject_vz <m_per_s>

Defaults:
vx = 0
vy = 0
vz = -0.5

---

### Multi-Axis Shaking

--shake_amp_x <meters>
--shake_amp_y <meters>
--shake_amp_z <meters>

Acceleration applied as:
a_i = A_i * sin(2*pi*shake_hz*t)

---

### Legacy XY Shake Compatibility

--shake_xy_legacy 1

Preserves historical behavior where vertical shake acceleration was also applied in x and y.
Default: 0

---

## Output Files

atoms_<iter>.xyz  element x y z r
atom.<iter>.vtk   VTK with velocity vectors
stdout            Iter, time, particle count, mean height, phi

---

## License

Free as in beer.
