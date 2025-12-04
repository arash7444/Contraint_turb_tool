#!/usr/bin/env python
"""
Multi-height LiDAR/SCADA-constrained Mann turbulence using Hipersim.

- Reads a 10 min MAT file (H2A_YYYY-MM-DD_HH-MM-SS.mat)
- Extracts several L_WS_* signals (heights)
- Builds u'(t) per height
- Generates a Mann field and constrains it at all heights simultaneously
- Compares LiDAR vs constrained field
- Exports constrained / unconstrained boxes to Bladed + HAWC2
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import numpy.core.numerictypes as ntypes
import pandas as pd
from hipersim.mann_turbulence import MannTurbulenceField


# ================================================================
# 0. NumPy 2.x compatibility patch  (Hipersim uses 'cfloat')
# ================================================================
if "cfloat" not in ntypes.sctypeDict:
    ntypes.sctypeDict["cfloat"] = np.complex128



# ================================================================
# 2. Helper functions
# ================================================================

def read_matfile(matfile: str) -> pd.DataFrame:
    """
    Read your H2A_YYYY-MM-DD_HH-MM-SS.mat file, return DataFrame with
    a Time index in seconds (10-min segment).
    """
    from scipy.io import loadmat

    mat = loadmat(matfile, squeeze_me=True)
    data = mat["DATA"]

    all_fields = data.dtype.names
    df = pd.DataFrame({name: data[name].item() for name in all_fields})

    # Sampling
    sample_length = df.shape[0]
    sample_freq = sample_length / (60 * 10)  # assume 10 min
    dt = 1.0 / sample_freq

    # Build absolute time from filename
    filename = Path(matfile).stem  # e.g. H2A_2025-07-07_16-10-00
    _, date_str, time_str = filename.split("_")
    time_str_clean = time_str.replace("-", ":")  # 16:10:00

    start_time = pd.to_datetime(
        f"{date_str} {time_str_clean}",
        format="%Y-%m-%d %H:%M:%S",
    )

    df["Time"] = start_time + pd.to_timedelta(df.index * dt, unit="s")
    df = df.set_index("Time")
    return df


def find_ws_columns_and_heights(df: pd.DataFrame):
    """
    Find all L_WS_* columns and parse their heights as integers.
    Returns:
        ws_cols_all: list of column names
        heights_all: list of heights [m] (same length)
    """
    ws_cols_all = [c for c in df.columns if c.startswith("L_WS_")]
    heights_all = []
    for col in ws_cols_all:
        m = re.search(r"L_WS_(\d+)", col)
        if m:
            heights_all.append(int(m.group(1)))
        else:
            heights_all.append(None)
    return ws_cols_all, heights_all


def prepare_multilevel_timeseries(df: pd.DataFrame,
                                  ws_cols_all,
                                  heights_all,
                                  col_indices,
                                  block_duration_min):
    """
    Selects a 10-min block and returns:
    - times (common to all selected heights)
    - ws_cols: selected column names
    - heights: corresponding heights [m]
    - U_mean_h: dict[height] = mean
    - u_prime_h: dict[height] = fluctuations
    """
    # select which WS columns & heights to use
    ws_cols = [ws_cols_all[i] for i in col_indices]
    heights = [heights_all[i] for i in col_indices]

    # 10-min block from start of file
    t0 = df.index.min()
    t1 = t0 + pd.Timedelta(minutes=block_duration_min)
    block = df[(df.index >= t0) & (df.index < t1)].copy()

    sub = block[ws_cols].copy()
    sub["Time"] = sub.index

    # require valid data at all selected heights
    mask_valid = np.ones(len(sub), dtype=bool)
    for c in ws_cols:
        mask_valid &= sub[c].notna().values
    sub = sub[mask_valid].reset_index(drop=True)

    times = sub["Time"].copy()

    print(f"10-min block: {t0} – {t1}")
    print(f"N samples (common to all heights) = {len(sub)}")

    U_mean_h = {}
    u_prime_h = {}
    for col, h in zip(ws_cols, heights):
        ws = sub[col].values
        U_mean_h[h] = ws.mean()
        u_prime_h[h] = ws - U_mean_h[h]
        print(
            f"  Height {h:3d} m: U_mean = {U_mean_h[h]:.3f} m/s, "
            f"std(u') = {u_prime_h[h].std():.3f} m/s"
        )

    return times, ws_cols, heights, U_mean_h, u_prime_h


def generate_mann_fields(NX, NY, NZ,DX, DY, DZ):
    """
    Generate unconstrained and constrained Mann boxes with the same parameters.
    """
    Nxyz = (NX, NY, NZ)
    dxyz = (DX, DY, DZ)
    print("\nMann box:")
    print("  Nxyz =", Nxyz)
    print("  dxyz =", dxyz)
    print(f"  Lx={NX*DX:.1f} m, Ly={NY*DY:.1f} m, Lz={NZ*DZ:.1f} m")

    mtf_nocon = MannTurbulenceField.generate(
        alphaepsilon=1.0,
        L=33.6,
        Gamma=3.9,
        Nxyz=Nxyz,
        dxyz=dxyz,
        seed=1,
        HighFreqComp=0,
        double_xyz=(False, True, True),
        n_cpu=1,
        verbose=0,
        cache_spectral_tensor=False,
    )

    mtf_con = MannTurbulenceField.generate(
        alphaepsilon=1.0,
        L=33.6,
        Gamma=3.9,
        Nxyz=Nxyz,
        dxyz=dxyz,
        seed=1,
        HighFreqComp=0,
        double_xyz=(False, True, True),
        n_cpu=1,
        verbose=0,
        cache_spectral_tensor=False,
    )

    return mtf_nocon, mtf_con


def build_constraints(mtf_con,
                      times,
                      heights,
                      U_mean_h,
                      u_prime_h,
                      NX, NY, NZ,DX, DY, DZ,
                      hub_height):
    """
    Build constraint array for all heights, using Taylor's hypothesis with
    U_mean_ref = U_mean(hub_height).
    Returns:
        x_pos, t_rel, y_pos, constraints_all
    """
    Lx = NX * DX
    Ly = NY * DY
    Lz = NZ * DZ

    # time → x positions
    t_rel = (times - times.iloc[0]).dt.total_seconds().values
    U_mean_ref = U_mean_h[hub_height]
    print(f"\nUsing U_mean_ref = U_mean({hub_height} m) = {U_mean_ref:.3f} m/s for advection")

    x_pos = U_mean_ref * t_rel
    mask = x_pos < Lx
    x_pos = x_pos[mask]
    t_rel = t_rel[mask]
    for h in heights:
        u_prime_h[h] = u_prime_h[h][mask]

    print(f"\nMax x_pos from LiDAR/SCADA = {x_pos.max():.1f} m (Lx = {Lx:.1f} m)")
    print(f"N samples after x clipping = {len(x_pos)}")

    # fixed lateral line
    y0 = 0.5 * Ly
    y_pos = np.full_like(x_pos, y0)
    print(f"Constraint line in y: y0 = {y0:.1f} m")

    constraints_list = []
    for h in heights:
        if h > Lz:
            raise RuntimeError(f"Height {h} m exceeds box vertical extent Lz={Lz} m")

        z_pos_h = np.full_like(x_pos, h, dtype=float)

        # Background v,w before constraining
        uvw_bg = mtf_con(x_pos, y_pos, z_pos_h)
        v_bg = uvw_bg[:, 1]
        w_bg = uvw_bg[:, 2]

        cons_h = np.column_stack(
            [
                x_pos,        # x [m]
                y_pos,        # y [m]
                z_pos_h,      # z [m] (height)
                u_prime_h[h], # u' target
                v_bg,         # keep original v
                w_bg,         # keep original w
            ]
        )
        constraints_list.append(cons_h)
        print(f"  Added constraints for height {h} m: shape {cons_h.shape}")

    constraints_all = np.vstack(constraints_list)
    print("\nTotal constraints shape:", constraints_all.shape)

    return x_pos, t_rel, y_pos, constraints_all


def verify_per_height(mtf_con, x_pos, y_pos, heights, u_prime_h):
    """
    Print std, corr, and max error per height after constraining.
    """
    results = {}
    for h in heights:
        z_pos_h = np.full_like(x_pos, h, dtype=float)
        uvw_con_h = mtf_con(x_pos, y_pos, z_pos_h)
        u_con_h = uvw_con_h[:, 0]

        std_lidar = u_prime_h[h].std()
        std_field = u_con_h.std()
        corr = np.corrcoef(u_prime_h[h], u_con_h)[0, 1]
        max_err = np.max(np.abs(u_con_h - u_prime_h[h]))

        results[h] = (std_lidar, std_field, corr, max_err)

        print(f"\nHeight {h} m:")
        print(f"  std(u'_LiDAR/SCADA)  = {std_lidar:.3f} m/s")
        print(f"  std(u'_constrained)  = {std_field:.3f} m/s")
        print(f"  corr(LiDAR, field)   = {corr:.3f}")
        print(f"  max |u_field - u_LiDAR| = {max_err:.3f} m/s")

    return results


def plot_multilevel(mtf_con, x_pos, t_rel, y_pos, heights, u_prime_h,NX, NY, NZ,DX, DY, DZ,
                    upsample_for_plot=True, upsample_factor=10):
    """
    Quick visual check: LiDAR/SCADA u'(t) vs constrained u'(t) at each height.
    Optional: also eval field on a finer time grid for smoother lines.
    """
    Lx = NX * DX
    plt.figure(figsize=(10, 6))

    if upsample_for_plot:
        t_rel_new = np.linspace(t_rel[0], t_rel[-1], num=len(t_rel) * upsample_factor)
        x_pos_new = np.interp(t_rel_new, t_rel, x_pos)
        y_pos_new = np.full_like(x_pos_new, y_pos[0])
    else:
        t_rel_new = None
        x_pos_new = None
        y_pos_new = None

    for i, h in enumerate(heights):
        z_pos_h = np.full_like(x_pos, h, dtype=float)
        uvw_con_h = mtf_con(x_pos, y_pos, z_pos_h)
        u_con_h = uvw_con_h[:, 0]

        plt.subplot(len(heights), 1, i + 1)
        plt.plot(
            t_rel,
            u_prime_h[h],
            label=f"LiDAR/SCADA u'(t) at {h} m",
            linewidth=1.2,
            marker="o",
            markersize=2,
        )
        plt.plot(
            t_rel,
            u_con_h,
            label=f"Constrained u'(t) at {h} m",
            alpha=0.7,
            marker="x",
            markersize=2,
        )

        if upsample_for_plot:
            z_pos_h_new = np.full_like(x_pos_new, h, dtype=float)
            uvw_con_h_new = mtf_con(x_pos_new, y_pos_new, z_pos_h_new)
            u_con_h_new = uvw_con_h_new[:, 0]
            plt.plot(
                t_rel_new,
                u_con_h_new,
                label=f"Constrained u'(t) at {h} m (interp)",
                alpha=0.7,
            )

        plt.ylabel("u' [m/s]")
        plt.grid(True, linestyle="--", alpha=0.3)
        plt.legend(loc="upper right")
        if i == len(heights) - 1:
            plt.xlabel("Time since block start [s]")

    plt.suptitle("Multi-height LiDAR/SCADA-constrained Mann turbulence")
    plt.tight_layout()
    plt.show()


def export_boxes(mtf_con, mtf_nocon, U_mean_ref):
    """
    Export constrained and unconstrained boxes to Bladed and HAWC2.
    """
    # Bladed-style export
    mtf_con.to_bladed(
        U=BLADED_U_MEAN_FOR_EXPORT,
        folder=EXPORT_FOLDER_BLADED_CON,
        basename=None,
    )
    mtf_nocon.to_bladed(
        U=BLADED_U_MEAN_FOR_EXPORT,
        folder=EXPORT_FOLDER_BLADED_NOCON,
        basename=None,
    )

    # HAWC2-style export (turbulence boxes)
    htc_wind_section_con = mtf_con.to_hawc2(
        folder=EXPORT_FOLDER_HAWC2_CON,
        basename=None,
        htc_dict={
            "wind.wsp": HAWC2_WSP_FOR_EXPORT,
            "wind.center_pos0": HAWC2_CENTER_POS0,
        },
    )
    print("\nConstrained HAWC2 wind section:")
    print(htc_wind_section_con)

    htc_wind_section_nocon = mtf_nocon.to_hawc2(
        folder=EXPORT_FOLDER_HAWC2_NOCON,
        basename=None,
        htc_dict={
            "wind.wsp": HAWC2_WSP_FOR_EXPORT,
            "wind.center_pos0": HAWC2_CENTER_POS0,
        },
    )
    print("\nUnconstrained HAWC2 wind section:")
    print(htc_wind_section_nocon)


# ================================================================
# 3. MAIN SCRIPT
# ================================================================
if __name__ == "__main__":


    # ================================================================
    # 1. USER INPUT SECTION  (edit for your case)
    # ================================================================

    MAT_FILE = r"C:\Users\12650009\ArashData\Projects\47_Measurement_LoadValidation\Measurement_data\Mat_files\New folder\H2A_2025-07-07_16-10-00.mat"

    # indices of L_WS_ columns to use (after sorting)
    # e.g. if df columns contain: L_WS_44, L_WS_75, L_WS_100, L_WS_125, ...
    # and you want 75, 125, 175 m, you can use [1, 3, 5]
    WS_COL_INDICES = [1, 3, 5]

    BLOCK_DURATION_MIN = 10       # length of segment to use
    HUB_HEIGHT = 125              # height used for Taylor advection U_mean_ref

    # Mann box parameters
    NX, NY, NZ = 4096, 32, 32
    DX, DY, DZ = 2.0, 10.0, 5.0   # DZ chosen so Lz >= max(height)

    # Output options
    EXPORT_FOLDER_BLADED_CON   = "./constraint_wind/bladed_con"
    EXPORT_FOLDER_BLADED_NOCON = "./constraint_wind/bladed_NoCon"
    EXPORT_FOLDER_HAWC2_CON    = "./constraint_wind/HAWC2_con"
    EXPORT_FOLDER_HAWC2_NOCON  = "./constraint_wind/HAWC2_NoCon"

    BLADED_U_MEAN_FOR_EXPORT = 9.0      # for .wnd
    HAWC2_WSP_FOR_EXPORT     = 6.8      # wind.wsp
    HAWC2_CENTER_POS0        = (0, 0, -125)








    # --- Read MAT file ---
    df_raw = read_matfile(MAT_FILE)
    ws_cols_all, heights_all = find_ws_columns_and_heights(df_raw)

    # --- Prepare multi-height time series ---
    times, ws_cols, heights, U_mean_h, u_prime_h = prepare_multilevel_timeseries(
        df_raw,
        ws_cols_all,
        heights_all,
        WS_COL_INDICES,
        BLOCK_DURATION_MIN,
    )

    # Ensure hub_height is in the selected heights
    hub_height = HUB_HEIGHT
    if hub_height not in heights:
        hub_height = heights[len(heights) // 2]
        print(f"\n[Warning] HUB_HEIGHT not in selected heights, "
              f"falling back to {hub_height} m")

    # --- Generate Mann boxes ---
    mtf_nocon, mtf_con = generate_mann_fields(NX, NY, NZ,DX, DY, DZ)

    # --- Build & apply constraints ---
    x_pos, t_rel, y_pos, constraints_all = build_constraints(
        mtf_con,
        times,
        heights,
        U_mean_h,
        u_prime_h,
        NX, NY, NZ,DX, DY, DZ,
        hub_height,
    )

    mtf_con.constrain(constraints_all)
    print("\nconstrain() finished.")

    # --- Verification ---
    results = verify_per_height(mtf_con, x_pos, y_pos, heights, u_prime_h)

    # --- Plots ---
    plot_multilevel(mtf_con, x_pos, t_rel, y_pos, heights, u_prime_h,
                    upsample_for_plot=True, upsample_factor=10)

    # --- Export to Bladed / HAWC2 ---
    export_boxes(mtf_con, mtf_nocon, U_mean_h[hub_height])

    input("\nPress Enter to exit...")
