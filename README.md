# Contraint Turbulence Tool using HiperSim

Small Python package to generate **LiDAR/SCADA-constrained Mann turbulence boxes**
using [Hipersim](https://gitlab.windenergy.dtu.dk/HiperSim/hipersim).

The tool:
- reads 10-min high-frequency measurement data from `.mat` files, LiDAR CSV files
- extracts wind speed time series at selected heights,
- converts time → spatial coordinates using Taylor’s hypothesis,
- calls `MannTurbulenceField.constrain()` at multiple heights, and
- exports Mann boxes to **HAWC2** and **Bladed** formats.

## Install (local, development mode)

```bash
git clone https://github.com/arash7444/Contraint_turb_tool.git
cd Contraint_turb_tool
pip install -e .    # installs constraint_turb_tool in editable mode
