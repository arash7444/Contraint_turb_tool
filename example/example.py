
import constraint_turb_tool as ctt


# ================================================================
# 1. USER INPUT SECTION  (edit for your case)
# ================================================================

MAT_FILE = r"C:\Users\12650009\ArashData\Projects\47_Measurement_LoadValidation\Measurement_data\Mat_files\New folder\H2A_2025-07-07_16-10-00.mat"

# indices of L_WS_ columns to use (after sorting)
# e.g. if df columns contain: L_WS_44, L_WS_75, L_WS_100, L_WS_125, ...
# and you want 75, 125, 175 m, you can use [1, 3, 5]
# WS_COL_INDICES = [1, 3, 5]
WS_COL_INDICES = [3]

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
df_raw = ctt.read_matfile(MAT_FILE)
ws_cols_all, heights_all = ctt.find_ws_columns_and_heights(df_raw)

# --- Prepare multi-height time series ---
times, ws_cols, heights, U_mean_h, u_prime_h = ctt.prepare_multilevel_timeseries(
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
mtf_nocon, mtf_con = ctt.generate_mann_fields(NX, NY, NZ,DX, DY, DZ)

# --- Build & apply constraints ---
x_pos, t_rel, y_pos, constraints_all = ctt.build_constraints(
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
results = ctt.verify_per_height(mtf_con, x_pos, y_pos, heights, u_prime_h)

# --- Plots ---
ctt.plot_multilevel(mtf_con, x_pos, t_rel, y_pos, heights, u_prime_h,NX, NY, NZ,DX, DY, DZ,
                upsample_for_plot=True, upsample_factor=10)

# --- Export to Bladed / HAWC2 ---
ctt.export_boxes(mtf_con, mtf_nocon, U_mean_h[hub_height])

input("\nPress Enter to exit...")
