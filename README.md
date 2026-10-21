# Drone Flight Plan GUI

Generates KMZ (DJI waypoints with turns) + PDF flight plan.
Methane measurement – DJI M400 – EU A1/A3 Open Category.

## Requirements

- Python 3.11
- `fpdf2`
- `pyproj`

## Setup (Conda)

```bash
conda env create -f environment.yml
conda activate flightplan