# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Thermoembo is a research codebase for **thermal and perfusion analysis in biological tissue with embedded vessels**, targeting oncology thermal ablation workflows. It combines:

1. **Medical image processing** — CT/MRI vessel segmentation using Hessian-based vesselness filters
2. **Mesh generation** — image-to-FEM mesh pipeline (tetgen → Exodus format)
3. **Coupled FEM simulation** — C++ solver (PETSc + VTK) for pressure, saturation, temperature, phase, and damage fields

## Build

### FEM Solver (C++)

```bash
cd couplemesh/
make thermoembovessel
```

Requires PETSc and VTK 5.10. The resulting binary is named `thermoembo-${PETSC_ARCH}`.

### Image Processing Pipelines

Driven by GNU Make with `c3d`, ANTS, and ITK tools:

```bash
# Root: vessel segmentation from raw CT
make vessel

# Venogram workflow
cd venogram/ && make
```

## Docker

Build the image once:

```bash
docker build -t thermoembo .
```

Run the full pipeline — always pass `--user $(id -u):$(id -g)` so output files are owned by the host user, not root:

```bash
docker run --rm --user $(id -u):$(id -g) \
  -v $(pwd):/data thermoembo \
  /data/004_vessel_phantom.nii.gz \
  --out-dir /data/thermoembo_run --steps 5
```

## Running the Solver

```bash
mpirun -n 12 ./thermoembo-${PETSC_ARCH} -dim 3 \
  -ts_max_steps 11150 -ts_dt 1.e-1 \
  -temp_petscspace_degree 1 -pres_petscspace_degree 1 \
  -mesh ./mytetmesh.1.exo -vtk ./output.vtk
```

## Mesh Generation

```bash
# Image → tetrahedral mesh → Exodus format
python couplemesh/phasemesh.py --file_name vessel.nii.gz --output mytetmesh
```

## Parameter Sweeps (Validation)

There is no formal test suite. Validation is done via parameter sweeps comparing vessel segmentation to ground truth:

```bash
python parametersweep.py   # runs pipeline variants, reports Dice overlap
```

## Architecture

### Data Flow

```
DICOM CT/MRI
    ↓  (c3d, ANTS)
vessel.nii.gz  (Hessian vesselness → Otsu threshold → binary mask)
    ↓  (phasemesh.py, tetgen)
mytetmesh.1.exo  (Exodus II mesh with vessel node sets)
    ↓  (thermoembovessel.c + PETSc)
output.vtk  (pressure / saturation / temperature / phase / damage fields)
```

### Key Files

| File | Role |
|------|------|
| `couplemesh/thermoembovessel.c` | Main FEM solver (~4100 lines): IMEX time-stepping, SNES Newton-LS, field-split preconditioner |
| `couplemesh/phasemesh.py` | Mesh generation orchestration: image → distance transform → tetgen → Exodus |
| `examplemesh/combinemesh.py` | Sphere geometry mesh example (Cubit-based) |
| `makefile` | Root image processing pipeline |
| `venogram/makefile` | Automated vessel detection and segmentation |
| `Params.yaml` | PyRadiomics configuration |
| `builddb.sql` | DICOM database queries for study organization |

### FEM Solver Fields

The coupled system solves five fields simultaneously:
- **u** — pressure
- **s** — saturation
- **temp** — temperature
- **phas** — phase
- **damg** — damage

Time integration uses backward Euler (implicit); the weak formulation covers advection-diffusion-reaction.

## Dependencies

- **C++**: PETSc, VTK 5.10
- **Python**: `nibabel`, `numpy`, `vtk`, `SimpleITK`, `radiomics` (PyRadiomics)
- **CLI tools**: `c3d`, `tetgen`, ANTS, ITK-Snap
- **Mesh tools**: Cubit (for example meshes)

Note: Python scripts use Python 2 syntax in some places (`.iteritems()`, bare `print` statements).
