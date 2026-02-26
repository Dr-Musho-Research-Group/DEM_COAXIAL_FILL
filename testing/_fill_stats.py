#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def read_xyz(path: Path) -> pd.DataFrame:
    with open(path, "r") as f:
        n = int(f.readline().strip())
        _comment = f.readline().strip()
        rows = []
        for _ in range(n):
            parts = f.readline().split()
            if len(parts) < 5:
                continue
            elem = parts[0]
            x, y, z, a = map(float, parts[1:5])  # a is particle radius
            rows.append((elem, x, y, z, a))
    df = pd.DataFrame(rows, columns=["elem", "x", "y", "z", "a"])
    df["r"] = np.sqrt(df["x"]**2 + df["y"]**2)
    return df

def binned_stats(x, y, nbins):
    edges = np.linspace(np.min(x), np.max(x), nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mean = np.full(nbins, np.nan)
    std = np.full(nbins, np.nan)
    count = np.zeros(nbins, dtype=int)

    for i in range(nbins):
        m = (x >= edges[i]) & (x < edges[i + 1])
        yi = y[m]
        count[i] = yi.size
        if yi.size >= 5:
            mean[i] = float(np.mean(yi))
            std[i] = float(np.std(yi, ddof=1))
    cv = std / mean
    return centers, mean, std, cv, count

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xyz", type=str, help="Path to .xyz file")
    ap.add_argument("--rin", type=float, default=22e-6, help="Inner radius [m]")
    ap.add_argument("--rout", type=float, default=40e-6, help="Outer radius [m]")
    ap.add_argument("--length", type=float, default=380e-6, help="Bed length [m] (for reference)")
    ap.add_argument("--rbins", type=int, default=30, help="Number of radial bins")
    ap.add_argument("--zbins", type=int, default=30, help="Number of axial bins")
    ap.add_argument("--out", type=str, default="fill_stats", help="Output prefix (no extension)")
    args = ap.parse_args()

    xyz_path = Path(args.xyz)
    df = read_xyz(xyz_path)

    # Global stats
    a_mean = df["a"].mean()
    a_std  = df["a"].std(ddof=1)

    # Basic “segregation index” vs radius: (local_mean - global_mean)/global_std
    r_cent, r_mean, r_std, r_cv, r_count = binned_stats(df["r"].to_numpy(), df["a"].to_numpy(), args.rbins)
    seg_idx = (r_mean - a_mean) / a_std

    # Axial stats
    z_cent, z_mean, z_std, z_cv, z_count = binned_stats(df["z"].to_numpy(), df["a"].to_numpy(), args.zbins)

    # --- Plot 1: radial profiles ---
    plt.figure()
    plt.plot(r_cent * 1e6, r_mean * 1e6, marker="o", linewidth=1)
    plt.xlabel("r [µm]")
    plt.ylabel("Mean particle radius [µm]")
    plt.grid(True, which="both", linestyle=":")
    plt.tight_layout()
    plt.savefig(f"{args.out}_radial_mean.png", dpi=200)

    plt.figure()
    plt.plot(r_cent * 1e6, r_cv, marker="o", linewidth=1)
    plt.xlabel("r [µm]")
    plt.ylabel("CV of particle radius in bin [-]")
    plt.grid(True, which="both", linestyle=":")
    plt.tight_layout()
    plt.savefig(f"{args.out}_radial_cv.png", dpi=200)

    plt.figure()
    plt.plot(r_cent * 1e6, seg_idx, marker="o", linewidth=1)
    plt.xlabel("r [µm]")
    plt.ylabel("Segregation index (mean-a_global)/std_global [-]")
    plt.grid(True, which="both", linestyle=":")
    plt.tight_layout()
    plt.savefig(f"{args.out}_radial_segregation_index.png", dpi=200)

    # --- Plot 2: axial profiles ---
    plt.figure()
    plt.plot(z_cent * 1e6, z_mean * 1e6, marker="o", linewidth=1)
    plt.xlabel("z [µm]")
    plt.ylabel("Mean particle radius [µm]")
    plt.grid(True, which="both", linestyle=":")
    plt.tight_layout()
    plt.savefig(f"{args.out}_axial_mean.png", dpi=200)

    plt.figure()
    plt.plot(z_cent * 1e6, z_cv, marker="o", linewidth=1)
    plt.xlabel("z [µm]")
    plt.ylabel("CV of particle radius in bin [-]")
    plt.grid(True, which="both", linestyle=":")
    plt.tight_layout()
    plt.savefig(f"{args.out}_axial_cv.png", dpi=200)

    # --- Plot 3: 2D heatmap mean a(r,z) ---
    r_edges = np.linspace(df["r"].min(), df["r"].max(), args.rbins + 1)
    z_edges = np.linspace(df["z"].min(), df["z"].max(), args.zbins + 1)
    df["rbin"] = np.digitize(df["r"], r_edges) - 1
    df["zbin"] = np.digitize(df["z"], z_edges) - 1

    A = np.full((args.zbins, args.rbins), np.nan)
    C = np.zeros((args.zbins, args.rbins), dtype=int)

    for (zb, rb), g in df.groupby(["zbin", "rbin"]):
        if 0 <= zb < args.zbins and 0 <= rb < args.rbins and len(g) >= 5:
            A[zb, rb] = g["a"].mean()
            C[zb, rb] = len(g)

    plt.figure()
    plt.imshow(A * 1e6, origin="lower", aspect="auto",
               extent=[r_edges[0]*1e6, r_edges[-1]*1e6, z_edges[0]*1e6, z_edges[-1]*1e6])
    plt.colorbar(label="Mean particle radius [µm]")
    plt.xlabel("r [µm]")
    plt.ylabel("z [µm]")
    plt.tight_layout()
    plt.savefig(f"{args.out}_heatmap_mean_a.png", dpi=200)

    # --- Save a CSV summary of the binned stats ---
    out_csv = pd.DataFrame({
        "r_center_m": r_cent,
        "mean_a_m": r_mean,
        "std_a_m": r_std,
        "cv_a": r_cv,
        "count": r_count,
        "seg_index": seg_idx
    })
    out_csv.to_csv(f"{args.out}_radial_stats.csv", index=False)

    print("Wrote:")
    for fn in [
        f"{args.out}_radial_mean.png",
        f"{args.out}_radial_cv.png",
        f"{args.out}_radial_segregation_index.png",
        f"{args.out}_axial_mean.png",
        f"{args.out}_axial_cv.png",
        f"{args.out}_heatmap_mean_a.png",
        f"{args.out}_radial_stats.csv",
    ]:
        print("  ", fn)

if __name__ == "__main__":
    main()