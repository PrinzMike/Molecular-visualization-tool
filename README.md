# Molecular Visualization Tool

This project is an interactive Streamlit app to visualize and correlate properties of electrolyte molecules relevant to lithium-ion battery research. It uses RDKit for molecular processing, Py3Dmol for 3D visualization, and Streamlit for the web interface.

## Features
- Computes molecular weight of common electrolytes.
- Generates 2D and 3D visualizations.
- Plots molecular weight trends.
- Allows adding new molecules interactively.

## Files
- `molecular_visualization.py`: The main Python script with Streamlit interface.
- `electrolyte_properties.csv`: Dataset with molecular weights.
- `molecule_2d_grid.png`: 2D grid image of molecules.
- `mw_plot.png`: Plot of molecular weight vs. index.
- `3d_molecules.html`: Interactive 3D visualization (open in browser).

## Usage
Run with `streamlit run molecular_visualization.py` to launch the interactive app. Requires Python 3.x and libraries: `rdkit`, `py3dmol`, `numpy`, `matplotlib`, `pandas`, `streamlit`.

