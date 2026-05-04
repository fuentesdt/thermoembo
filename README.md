# thermoembo

Coupled 1D–3D finite-element solver for thermal ablation simulation in
biological tissue with embedded vessels.  The Docker image bundles the C++ FEM
solver (built against PETSc 3.10.2 and VTK 5.10) together with the Python
pipeline that prepares the mesh inputs from a labeled NIfTI volume.

## Prerequisites

- Docker (any recent version)
- A labeled NIfTI volume with the following label convention:

| Label | Tissue |
|-------|--------|
| 1 | Liver parenchyma |
| 2 | Inflow vessel lumen (embolic injection site) |
| 3 | Outflow vessel lumen (optional) |
| 4 | Tumor (target — temperature evaluated here) |

## Building the image

```bash
docker build -t thermoembo .
```

The build compiles **four binaries** inside the image:

| Binary | Flags | Purpose |
|--------|-------|---------|
| `thermoembo-3.10.2-xenial-gcc-5.4.0-opt` | `-O2` | Production solver (3-D) |
| `thermoembo1d-3.10.2-xenial-gcc-5.4.0-opt` | `-O2` | Production 1-D/3-D coupled solver |
| `thermoembo-3.10.2-xenial-gcc-5.4.0-opt-debug` | `-O0 -g3` | Debug solver (3-D) |
| `thermoembo1d-3.10.2-xenial-gcc-5.4.0-opt-debug` | `-O0 -g3` | Debug 1-D/3-D coupled solver |

All four live in `/work/tutorials/` inside the container.

---

## Running the full pipeline

Mount a directory containing the labeled NIfTI and pass it as the first argument:

```bash
docker run --rm \
  -v $(pwd):/data \
  thermoembo \
  /data/PreTxArtLoRes.vessellabel.nii.gz \
  --out-dir /data/thermoembo_run \
  --steps 5
```

**All options:**

```
usage: run_thermoembo.py [-h] [--phantom-dir DIR] [--phantom-ids IDS]
                         [--out-dir DIR] [--decimate R] [--maxvol MM3]
                         [--tet-dihedral DEG] [--steps N] [--dt S]
                         [--vesselcoupling G] [--mesh-only]
                         [nii]

  nii                 Labeled NIfTI (label 1=liver, 2=inflow vessel, 3=outflow vessel, 4=tumor)
  --phantom-dir DIR   Directory of phantom NIfTI files; process all matching *_vessel_phantom.nii.gz
  --phantom-ids IDS   Comma-separated phantom IDs to process (e.g. 001,004); requires --phantom-dir
  --out-dir DIR       Output directory (default: thermoembo_run)
  --decimate R        Surface decimation ratio, 0–1 (default 0.90 = keep 10%)
  --maxvol MM3        TetGen max tet volume mm³ (default 5000)
  --tet-dihedral DEG  TetGen min dihedral angle (default 15)
  --steps N           Time steps (default 5)
  --dt S              Time-step size in seconds (default 60)
  --vesselcoupling G  Gaussian 1D–3D coupling amplitude (default 1e-3)
  --mesh-only         Prepare meshes only; skip simulation
```

Prepare meshes only (no simulation):

```bash
docker run --rm -v $(pwd):/data thermoembo \
  /data/label.nii.gz --out-dir /data/out --mesh-only
```

**Outputs written to `--out-dir`:**

| File | Description |
|------|-------------|
| `vessel_skeleton.nii.gz` | Skeletonised label-2 volume |
| `centerline.vtp` | 1-D vessel centerline mesh |
| `surface.vtp` | Decimated liver surface |
| `mesh.vtu` | Tetrahedral volume mesh (mm) |
| `mesh.exo` | Exodus II mesh (metres, for PETSc) |
| `vessel.vtk` | Zero-filled phase-field initial image |
| `resultsolution000.NNNN.vtu` | Solution snapshots per time step |

### Key solver parameters

The pipeline passes the following fixed PETSc flags to `thermoembo1d`:

| Flag | Value | Effect |
|------|-------|--------|
| `-gamma` | `0.087` | Heat-of-reaction rate constant (10× default); drives tumour heating |
| `-vesselcoupling` | `1e-3` | Gaussian 1D–3D coupling amplitude; set via `--vesselcoupling` |
| `-baselinepressure` | `1.0` | Tissue pressure boundary condition (atm) |
| `-snes_type` | `ksponly` | Single linearisation per time step (backward Euler) |
| `-pc_type` | `ilu` | ILU(2) preconditioner with NONZERO diagonal shift |

**Expected output:** with default settings (5 steps, dt=60 s) the mean temperature
in the label==4 tumour region increases by ≥10 % above the 37 °C body temperature
baseline (0.37 hK → ≥ 0.407 hK).

---

## Interactive shell

Override the entrypoint to drop into a shell:

```bash
docker run --rm -it \
  --entrypoint /bin/bash \
  -v $(pwd):/data \
  thermoembo
```

---

## Debugging with gdb

### 1 — Start an interactive container

```bash
docker run --rm -it \
  --entrypoint /bin/bash \
  --cap-add=SYS_PTRACE \       # required for gdb / valgrind
  -v $(pwd):/data \
  thermoembo
```

> `--cap-add=SYS_PTRACE` grants the ptrace capability gdb needs to attach to
> processes and read core files.

### 2 — Run the mesh preparation step first

Inside the container, run the Python pipeline with `--mesh-only` to generate
the Exodus and VTP inputs the solver needs:

```bash
python /work/pipeline/run_thermoembo.py \
  /data/PreTxArtLoRes.vessellabel.nii.gz \
  --out-dir /data/dbg_run \
  --mesh-only
```

### 3 — Start the debug build under gdb

The debug binary (`-O0 -g3`, no inlining) lets gdb step through source lines
without optimisation artefacts:

```bash
# set the output directory so the solver can find mesh.exo / centerline.vtp
cd /data/dbg_run

gdb /work/tutorials/thermoembo1d-3.10.2-xenial-gcc-5.4.0-opt-debug
```

Inside gdb:

```gdb
(gdb) set args -dim 3 \
               -mesh /data/dbg_run/mesh.exo \
               -vtp1d_in /data/dbg_run/centerline.vtp \
               -vesselcoupling 1e-3 \
               -vtk  /data/dbg_run/vessel.vtk \
               -o    /data/dbg_run/result \
               -temp_petscspace_degree 1 -pres_petscspace_degree 1 \
               -ts_type beuler -ts_max_steps 3 -ts_dt 60.0 \
               -gamma 0.087 \
               -snes_type ksponly -ksp_type gmres \
               -pc_type ilu -pc_factor_levels 2 \
               -pc_factor_shift_type NONZERO -pc_factor_shift_amount 1e-10

(gdb) break ThermoEmbodiesRHSFunction   # break on a function in thermoembo1d.c
(gdb) run
(gdb) backtrace
(gdb) list
(gdb) next
(gdb) print pressure
```

### 4 — cgdb (vi-keybinding TUI)

`cgdb` splits the terminal: source view on top, gdb prompt below.
Vim motion keys navigate the source, `i` switches to gdb input:

```bash
cgdb /work/tutorials/thermoembo1d-3.10.2-xenial-gcc-5.4.0-opt-debug
```

Key bindings inside cgdb:

| Key | Action |
|-----|--------|
| `Esc` | Move focus to source window |
| `i` | Move focus to gdb window |
| `k` / `j` | Scroll source up / down |
| `Ctrl+]` (gdb) | Follow tag (if ctags loaded) |
| `/` | Search in source window |
| `spacebar` | Toggle breakpoint on current line |

### 5 — tmux layout for simultaneous editing and debugging

```bash
tmux new-session -s dbg
# split horizontally: Ctrl+b then %
# left pane: cgdb; right pane: vim
```

---

## Vim + ctags / cscope navigation

### Generate tags and cscope database

Inside the container, from `/work/tutorials`:

```bash
cd /work/tutorials

# ctags — symbol index (functions, variables, structs)
ctags -R --languages=c,c++ --langmap=c++:+.txx \
      ${PETSC_DIR} /usr/include/vtk-5.10/ .

# cscope — cross-reference database (callers, callees, text search)
find . ${PETSC_DIR}/include /usr/include/vtk-5.10 \
     -name "*.c" -o -name "*.h" -o -name "*.cpp" \
     > cscope.files
cscope -b -q -k
```

### vim commands for C++ navigation

| Command | Action |
|---------|--------|
| `Ctrl+]` | Jump to definition under cursor |
| `Ctrl+t` | Jump back |
| `:tag funcname` | Jump to `funcname` |
| `:tnext` / `:tprev` | Cycle through multiple matches |
| `Ctrl+\`s` | cscope: find all uses of symbol |
| `Ctrl+\`g` | cscope: find definition |
| `Ctrl+\`c` | cscope: find callers |
| `Ctrl+\`f` | cscope: find file by name |
| `:e .` | Open directory browser |

### Open the solver source in vim

```bash
vim /work/tutorials/thermoembo1d.c
```

---

## Memory debugging with valgrind

```bash
valgrind --leak-check=full --track-origins=yes \
  /work/tutorials/thermoembo1d-3.10.2-xenial-gcc-5.4.0-opt-debug \
  -dim 3 -mesh /data/dbg_run/mesh.exo \
  -vtp1d_in /data/dbg_run/centerline.vtp \
  -vtk /data/dbg_run/vessel.vtk \
  -o /data/dbg_run/result \
  -temp_petscspace_degree 1 -pres_petscspace_degree 1 \
  -ts_type beuler -ts_max_steps 2 -ts_dt 60.0 \
  -gamma 0.087 \
  -snes_type ksponly -ksp_type gmres \
  -pc_type ilu -pc_factor_levels 2
```

> Use the debug binary with valgrind — the `-O0 -g3` flags give valgrind
> accurate file/line information.  Valgrind is very slow; use a tiny mesh and
> a small number of steps.

---

## Running the production solver directly

Use `mpirun` for multi-rank production runs:

```bash
mpirun --allow-run-as-root -n 4 \
  /work/tutorials/thermoembo1d-3.10.2-xenial-gcc-5.4.0-opt \
  -dim 3 \
  -mesh /data/mytetmesh.exo \
  -vtp1d_in /data/centerline.vtp \
  -vesselcoupling 1e-3 \
  -vtk /data/vessel.vtk \
  -o /data/result \
  -temp_petscspace_degree 1 -pres_petscspace_degree 1 \
  -ts_type beuler -ts_max_steps 11150 -ts_dt 60.0 \
  -gamma 0.087 \
  -snes_type ksponly -ksp_type gmres \
  -pc_type ilu -pc_factor_levels 2 \
  -pc_factor_shift_type NONZERO -pc_factor_shift_amount 1e-10 \
  -ts_monitor -modulowrite 100
```

---

## Architecture

```
docker run thermoembo <label.nii.gz> [options]
       │
       └─ pipeline/run_thermoembo.py
              │
              ├─ centerline.skel_to_vtp()   → centerline.vtp
              ├─ meshing.*                  → surface.vtp + mesh.vtu
              ├─ vtu_to_exodus()            → mesh.exo (metres)
              ├─ make_vessel_vtk()          → vessel.vtk (binary label-2 mask)
              └─ run_solver() ─────────────→ /work/tutorials/thermoembo1d-...-opt
                   single-pass coupled solve → resultsolution000.NNNN.vtu
```
