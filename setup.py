from setuptools import setup, find_packages

setup(
    name="constraint_turb_tool",
    version="0.1.0",
    description="Generate constrained Mann turbulence boxes from LiDAR/SCADA using Hipersim",
    author="Arash Atasen",
    packages=find_packages("src"),
    package_dir={"": "src"},
    python_requires=">=3.9",
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "scipy",
        "hipersim @ git+https://gitlab.windenergy.dtu.dk/HiperSim/hipersim.git",
    ],
)
