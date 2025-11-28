
---

## 5. Add an `examples/` script

Create `example/run_from_mat.py`:

```python
#!/usr/bin/env python
"""
Minimal example of building a constrained Mann box from one .mat file.
"""

from constraint_turb_tool import (
    read_matfile,
    find_ws_columns_and_heights,
    prepare_multilevel_timeseries,
    generate_mann_fields,
    build_constraints,
    verify_per_height,
    plot_multilevel,
    export_boxes,
)

# --- user inputs (you can edit these) ----------------------------------------
MAT_FILE = r"C:\path\to\H2A_2025-07-07_16-10-00.mat"
WS_COL_INDICES = [1, 3, 5]   # which heights to use
BLOCK_DURATION_MIN = 10
HUB_HEIGHT = 125
# -----------------------------------------------------------------------------


def main():
    df_raw = read_matfile(MAT_FILE)
    ws_cols_all, heights_all = find_ws_columns_and_heights(df_raw)

    times, ws_cols, heights, U_mean_h, u_prime_h = prepare_multilevel_timeseries(
        df_raw,
        ws_cols_all,
        heights_all,
        WS_COL_INDICES,
        BLOCK_DURATION_MIN,
    )

    hub_height = HUB_HEIGHT if HUB_HEIGHT in heights else heights[len(heights) // 2]

    mtf_nocon, mtf_con = generate_mann_fields()

    x_pos, t_rel, y_pos, constraints_all = build_constraints(
        mtf_con, times, heights, U_mean_h, u_prime_h, hub_height
    )

    mtf_con.constrain(constraints_all)

    verify_per_height(mtf_con, x_pos, y_pos, heights, u_prime_h)
    plot_multilevel(mtf_con, x_pos, t_rel, y_pos, heights, u_prime_h)
    export_boxes(mtf_con, mtf_nocon, U_mean_h[hub_height])


if __name__ == "__main__":
    main()
