"""
run_thermoembo.py — Prepare and launch a thermoembo1d coupled 1D-3D simulation.

Pipeline
--------
  1a. Extract label==2 (inflow) mask  → skeletonize → {id}_inflow_skel.nii.gz
  1b. Extract label==3 (outflow) mask → skeletonize → {id}_outflow_skel.nii.gz
      (outflow skeleton is optional; step is skipped if label==3 is absent)
  2a. skel_to_vtp() on inflow  → raw inflow VTP graph
  2b. resistance_lumping.solve() on inflow  → {id}_inflow_centerline.vtp
      Node pressures set by Hagen-Poiseuille resistance lumping (see below).
  2c. skel_to_vtp() on outflow → raw outflow VTP graph
  2d. resistance_lumping.solve() on outflow → {id}_outflow_centerline.vtp
  3.  meshing.*  → surface.vtp + mesh.vtu   (3-D liver tet mesh, mm)
  4.  vtu_to_exodus()           → mesh.exo  (mm → m, for PETSc)
  5.  make_vessel_vtk()         → vessel.vtk (binary label-2 phase-field IC)
  6.  run_solver()              → resultsolution000.NNNN.vtu per time step

Inputs
------
  <nii>  — labeled NIfTI (.nii or .nii.gz)
    label == 1 : liver parenchyma  (3-D mesh domain)
    label == 2 : inflow vessel lumen  (hepatic artery / embolic injection site)
    label == 3 : outflow vessel lumen (portal vein / venous drainage), optional
    label == 4 : tumor target  (temperature and concentration evaluated here)

Resistance lumping (1-D solve)
------------------------------
  Node pressures along each vessel centerline are pre-computed by solving a
  sparse Hagen-Poiseuille conductance system  K·p = q  before the 3-D FEM
  solve begins.  This replaces the uniform-pressure Dirichlet BC with a
  physiologically consistent pressure gradient along the vessel tree.

  Method (resistance_lumping.py):
    1. Read the raw skeleton VTP graph (nodes = centerline points, edges = line
       cells).  Vessel radii at each node are derived from the distance
       transform of the binary label mask (scipy distance_transform_edt with
       true anisotropic voxel spacing).
    2. Disconnected sub-graphs are bridged with minimum-spanning-tree (MST)
       phantom edges connecting inter-component leaf pairs, ensuring a single
       connected conductance network with the fewest spurious parallel paths.
    3. Conductance of each edge: G_ij = π r⁴ / (8 μ L), where r = mean radius
       of the two endpoints, L = Euclidean edge length, μ = blood viscosity.
    4. Boundary conditions:
         Inflow  — root node (max-radius leaf): P = 836 mmHg (1.10 atm, HA)
                   injection-site seed (from .fcsv): P = 882 mmHg (1.16 atm)
                   terminal branches: P = 760 mmHg (1.00 atm, tissue)
         Outflow — max-radius leaf: P = 752 mmHg (0.99 atm, IVC)
                   terminal branches: P = 760 mmHg (1.00 atm, tissue)
    5. The linear system K·p = q is solved with scipy.sparse.linalg.spsolve.
       The resulting nodal pressures (mmHg) are written as a point-data array
       "Pressure_mmHg" to the output VTP and converted to atm when read by the
       3-D solver via -vtp1d_in / -vtp1d_out.

Outputs (written to --out-dir, default: thermoembo_run/)
  {id}_inflow_centerline.vtp    1-D inflow centerline with lumped pressures
  {id}_outflow_centerline.vtp   1-D outflow centerline with lumped pressures
  surface.vtp                   decimated liver surface
  mesh.vtu                      tetrahedral volume mesh (mm)
  mesh.exo                      Exodus II mesh (m, for PETSc/thermoembo1d)
  vessel.vtk                    binary vessel mask (phase-field IC, mm)
  resultsolution000.NNNN.vtu    solution snapshots per time step
  dashboard.html                Chart.js viewer: tumor T and concentration vs t

Usage
-----
  # Single NIfTI
  python run_thermoembo.py label.nii.gz --out-dir thermoembo_run --steps 5

  # Batch phantom directory
  python run_thermoembo.py --phantom-dir /data --phantom-ids 001,002,003,004 \\
      --out-dir thermoembo_run --steps 5

Dependencies
------------
  pip install nibabel pyvista pymeshfix tetgen scipy scikit-image meshio netCDF4
  solver binary: /work/tutorials/thermoembo1d-3.10.2-xenial-gcc-5.4.0-opt
"""

import argparse
import glob
import os
import subprocess
import sys

import numpy as np
import nibabel as nib
import vtk
from vtk.util.numpy_support import numpy_to_vtk
import pyvista as pv
import meshio
from scipy.ndimage import binary_fill_holes
from skimage.morphology import skeletonize

from centerline import skel_to_vtp
from meshing import load_nii, extract_surface, process_surface, generate_tet_mesh, write_outputs

import resistance_lumping


THERMOEMBO_EXE = "/work/tutorials/thermoembo1d-3.10.2-xenial-gcc-5.4.0-opt"


def log(msg):
    print(f"[run_thermoembo] {msg}", flush=True)


# ── Step 1: skeletonize vessel by label ───────────────────────────────────────

def make_label_skeleton(nii_path, out_nii, label_val=2):
    """Extract voxels where label == label_val, skeletonise, write NIfTI."""
    log(f"Skeletonising vessel mask (label=={label_val}) ...")
    nii    = nib.load(nii_path)
    data   = np.asarray(nii.dataobj)
    binary = (data == label_val)
    log(f"  Vessel voxels  : {int(binary.sum())}")
    if not binary.any():
        log(f"Warning: no label=={label_val} voxels found — skipping")
        return False

    skel = skeletonize(binary)          # Lee 1994 thinning (3-D)
    log(f"  Skeleton voxels: {int(skel.sum())}")

    out = nib.Nifti1Image(skel.astype(np.uint8), nii.affine, nii.header)
    nib.save(out, out_nii)
    log(f"  Written        : {out_nii}")
    return True


# ── Step 3: write legacy VTK image (material file) ────────────────────────────

def make_vessel_vtk(nii_path, out_vtk):
    """
    Write a zero-filled VTK legacy structured-points image covering the NIfTI
    extent.

    thermoembo1d uses this image via -vtk to initialise the phase field via
    analytic_phas(): phase(x) = image_value(x).  A zero image means
    "no embolic material deposited at t=0," which puts all tissue nodes in the
    free-DOF set (phase=0 < phasethresh=0.5).  Using the tissue label image
    would set phase=1 or 2 everywhere, making isnotdirichlet empty and leaving
    the coupled solver with no free degrees of freedom.

    vtkImageData requires positive axis-aligned spacing.  For NIfTI affines
    with negative diagonal elements (common RAS/LPS flip) we shift the origin
    to the minimum world corner so that ComputeStructuredCoordinates maps FEM
    world coordinates to correct voxel indices.
    """
    log("Writing VTK material image (vessel label==2 → phase=1.0) ...")
    nii    = nib.load(nii_path)
    affine = nii.affine.astype(np.float64)

    signs   = np.sign(np.diag(affine[:3, :3]))
    spacing = np.abs(np.diag(affine[:3, :3]))   # always positive

    # Origin = minimum world corner (shift for negative-diagonal axes)
    origin = affine[:3, 3].copy()
    for ax in range(3):
        if signs[ax] < 0:
            n = nii.shape[ax]
            origin[ax] = affine[ax, 3] + signs[ax] * spacing[ax] * (n - 1)
    origin = np.minimum(origin, affine[:3, 3])

    log(f"  Origin  (mm)   : {origin.round(2).tolist()}")
    log(f"  Spacing (mm)   : {spacing.round(4).tolist()}")
    log(f"  Dimensions     : {list(nii.shape)}")

    # Vessel mask: label==2 → phase=1.0 (Dirichlet BC nodes), else 0.0 (free DOFs).
    # Flip axes with negative diagonal so spacing is always positive in the VTK image.
    data = np.asarray(nii.dataobj).astype(np.float32)
    vessel_mask = (data == 2).astype(np.float32)
    for ax in range(3):
        if signs[ax] < 0:
            vessel_mask = np.flip(vessel_mask, axis=ax)

    img = vtk.vtkImageData()
    img.SetDimensions(nii.shape)
    img.SetSpacing(spacing.tolist())
    img.SetOrigin(origin.tolist())

    arr = numpy_to_vtk(vessel_mask.ravel(order='F'), deep=True, array_type=vtk.VTK_FLOAT)
    arr.SetName("phase")
    img.GetPointData().SetScalars(arr)

    writer = vtk.vtkDataSetWriter()
    writer.SetFileName(out_vtk)
    writer.SetInputData(img)
    writer.SetFileTypeToBinary()
    writer.Write()
    log(f"  Written        : {out_vtk}")


# ── Step 4: VTU → Exodus conversion (mm → m) ─────────────────────────────────

def vtu_to_exodus(vtu_path, exo_path):
    """
    Convert a pyvista VTU (tetrahedral, coordinates in mm) to Exodus II
    (coordinates in metres), as required by PETSc/thermoembo1d.
    """
    log("Converting VTU → Exodus (mm → m) ...")
    grid = pv.read(vtu_path)

    VTK_TETRA = 10
    if VTK_TETRA not in grid.cells_dict:
        sys.exit(f"Error: VTU mesh does not contain tetrahedral cells "
                 f"(found types: {list(grid.cells_dict.keys())})")

    points_m = grid.points * 1.0e-3               # mm → m
    cells    = [meshio.CellBlock('tetra', grid.cells_dict[VTK_TETRA])]
    mesh     = meshio.Mesh(points_m, cells)
    meshio.write(exo_path, mesh)

    log(f"  Nodes          : {grid.n_points}")
    log(f"  Elements       : {grid.n_cells}")
    bb = grid.bounds
    log(f"  BBox (m)       : x [{bb[0]*1e-3:.4f}, {bb[1]*1e-3:.4f}]  "
        f"y [{bb[2]*1e-3:.4f}, {bb[3]*1e-3:.4f}]  "
        f"z [{bb[4]*1e-3:.4f}, {bb[5]*1e-3:.4f}]")
    log(f"  Written        : {exo_path}")


# ── Step 5: run solver ────────────────────────────────────────────────────────

def _solver_cmd(abs_out, steps, dt, vesselcoupling,
                vtp1d_in=None, vtp1d_out=None):
    """Return the argument list for thermoembo1d (both passes use the same args)."""
    vtp_in = vtp1d_in or os.path.join(abs_out, "centerline.vtp")
    cmd = [
        THERMOEMBO_EXE,
        # --- mesh ---
        "-dim",  "3",
        "-mesh", os.path.join(abs_out, "mesh.exo"),
        # --- 1-D inflow vessel centerline ---
        "-vtp1d_in",       vtp_in,
        "-vesselcoupling", str(vesselcoupling),
        # --- material image ---
        "-vtk", os.path.join(abs_out, "vessel.vtk"),
        # --- output prefix (also defines phase-solution cache filename) ---
        "-o",   os.path.join(abs_out, "result"),
        # --- FE degree (P1 for all fields) ---
        "-temp_petscspace_degree", "1",
        "-pres_petscspace_degree", "1",
        "-damg_petscspace_degree", "1",
        "-conc_petscspace_degree", "1",
        "-phas_petscspace_degree", "1",
        # --- time integration ---
        "-ts_type",      "beuler",
        "-ts_max_steps", str(steps),
        "-ts_dt",        str(dt),
        "-ts_monitor",
        "-modulowrite",  "1",           # write VTU every time step
        # --- main solve ---
        # ILU(2) handles the large-dt (backward Euler) Jacobian of this
        # 5-field coupled system; without sufficient fill the factorisation
        # is too inaccurate and GMRES stalls.  The nonzero shift protects
        # against near-zero pivots on MatZeroRowsIS Dirichlet BC rows.
        "-snes_type",               "ksponly",
        "-ksp_type",                "gmres",
        "-ksp_gmres_restart",       "200",
        "-ksp_max_it",              "1000",
        "-ksp_rtol",                "1e-4",
        "-pc_type",                 "ilu",
        "-pc_factor_levels",        "2",
        "-pc_factor_shift_type",    "NONZERO",
        "-pc_factor_shift_amount",  "1e-10",
        "-ksp_converged_reason",
        "-snes_converged_reason",
        "-ts_max_snes_failures", "-1",
        # --- misc ---
        "-disppressure",     "0.0",
        "-baselinepressure", "1.0",
        "-gamma",            "0.087", # 10× default gammaconst → stronger heat source
    ]
    if vtp1d_out:
        cmd += ["-vtp1d_out", vtp1d_out]
    return cmd


def _rewrite_vtu_base64(path):
    """
    Re-write a VTU file written by VTK 5.10 (raw appended, may contain null
    bytes that are invalid in XML 1.0) to inline base64 binary so ParaView
    can open it without XML parse errors.  Overwrites the file.

    SetDataModeToBinary() writes each DataArray's data inline as base64
    (format="binary") rather than in an AppendedData section.  This avoids
    the VTK underscore-marker / Expat interaction that causes parse failures
    in ParaView 5.11 even with encoding="base64" appended mode.
    """
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(path)
    r.Update()
    w = vtk.vtkXMLUnstructuredGridWriter()
    w.SetFileName(path)
    w.SetInputData(r.GetOutput())
    w.SetDataModeToBinary()       # inline base64 per DataArray, no AppendedData
    w.SetCompressorTypeToNone()   # skip zlib — keeps XML structure simple
    w.Write()


def _rewrite_vtu_files(out_dir, pattern):
    files = sorted(glob.glob(os.path.join(out_dir, pattern)))
    if not files:
        return
    log(f"Converting {len(files)} VTU file(s) to base64 for ParaView ...")
    for f in files:
        _rewrite_vtu_base64(f)
    log(f"  Done.")


def _collect_label_timeseries(out_dir, nii_path, label_val, dt):
    """
    Read each resultsolution*.vtu once and collect temperature and concentration
    statistics for mesh nodes whose NIfTI voxel has label==label_val.

    Returns a dict with lists: times, temp_mean/min/max, conc_mean/min/max.
    Returns None if no VTU files exist or the label is absent.
    """
    vtu_files = sorted(glob.glob(os.path.join(out_dir, "resultsolution*.vtu")))
    if not vtu_files:
        return None

    nii     = nib.load(nii_path)
    labels  = np.asarray(nii.dataobj)
    if not (labels == label_val).any():
        log(f"Label=={label_val} not present in {os.path.basename(nii_path)}; "
            f"skipping evaluation")
        return None

    inv_aff = np.linalg.inv(nii.affine)   # affine is in mm; mesh coords converted to mm
    shape   = np.array(labels.shape)

    times = []
    temp_mean, temp_min, temp_max = [], [], []
    conc_mean, conc_min, conc_max = [], [], []

    for vtu_path in vtu_files:
        # parse step index from "resultsolution000.NNNN.vtu" → NNNN
        stem     = os.path.basename(vtu_path)
        step_idx = int(stem.split(".")[-2])
        t        = step_idx * dt

        m   = pv.read(vtu_path)
        pts = m.points * 1e3                  # m → mm
        ph  = np.hstack([pts, np.ones((len(pts), 1))])
        vox = (inv_aff @ ph.T).T[:, :3]
        idx = np.clip(np.round(vox).astype(int), 0, shape - 1)
        mask = labels[idx[:, 0], idx[:, 1], idx[:, 2]] == label_val
        if not mask.any():
            log(f"  {stem}: no label=={label_val} nodes in mesh")
            continue

        temp = m.point_data["solutiontemperature.0"][mask]
        conc = m.point_data["solutionconcentration.0"][mask]
        temp_all = m.point_data["solutiontemperature.0"]

        times.append(t)
        temp_mean.append(float(np.mean(temp)))
        temp_min.append(float(np.min(temp)))
        temp_max.append(float(np.max(temp)))
        conc_mean.append(float(np.mean(conc)))
        conc_min.append(float(np.min(conc)))
        conc_max.append(float(np.max(conc)))

        step_tag = stem.replace("resultsolution", "step").replace(".vtu", "")
        log(f"  {step_tag}: label{label_val} n={mask.sum():5d}  "
            f"temp mean={np.mean(temp):.4f} min={np.min(temp):.4f} max={np.max(temp):.4f}  "
            f"conc mean={np.mean(conc):.4f} min={np.min(conc):.4f} max={np.max(conc):.4f}  "
            f"| all_nodes temp max={np.max(temp_all):.4f} mean={np.mean(temp_all):.4f}")

    if not times:
        return None

    return dict(
        times=times,
        temp_mean=temp_mean, temp_min=temp_min, temp_max=temp_max,
        conc_mean=conc_mean, conc_min=conc_min, conc_max=conc_max,
    )


def _write_dashboard(out_path, named_series, label_val):
    """
    Write a self-contained HTML dashboard with two Chart.js line graphs per entry:
    temperature over time and concentration over time for the given label region.

    named_series: list of (name_str, series_dict) tuples.
      - Single run: pass [("", series)] — charts rendered without a phantom prefix.
      - Phantom-dir run: pass [("001", s1), ("002", s2), ...] — each pair is headed
        with the phantom ID, producing 2 × len(named_series) total charts.

    Chart.js is loaded from jsDelivr CDN — requires internet on first open.
    """
    import json

    entries_js = json.dumps([{"id": name, "series": s} for name, s in named_series])

    # Build canvas elements in Python so they exist before JS runs.
    canvas_html = ""
    for name, _ in named_series:
        safe_id     = (name or "run").replace(" ", "_")
        prefix      = f"Phantom {name} — " if name else ""
        canvas_html += (
            f'\n<h2>{prefix}Temperature over time (hK)</h2>\n'
            f'<canvas id="tempChart_{safe_id}"></canvas>\n'
            f'<h2>{prefix}Concentration over time</h2>\n'
            f'<canvas id="concChart_{safe_id}"></canvas>\n'
        )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Thermoembo Dashboard — Label {label_val}</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4/dist/chart.umd.min.js"></script>
<style>
  body {{
    font-family: sans-serif;
    max-width: 900px;
    margin: 2rem auto;
    padding: 0 1rem;
    color: #222;
  }}
  h1 {{ font-size: 1.4rem; margin-bottom: 0.25rem; }}
  h2 {{ font-size: 1.1rem; color: #444; margin-top: 2.5rem; }}
  canvas {{ margin-bottom: 1rem; }}
</style>
</head>
<body>
<h1>Thermoembo — Label {label_val} (tumor) region</h1>
<p style="color:#666;font-size:0.9rem">
  Shaded band = min/max across label nodes. Line = mean.
</p>
{canvas_html}
<script>
const entries = {entries_js};

function makeChart(canvasId, yLabel, times, meanData, minData, maxData, rgbColor) {{
  const alpha = 'rgba(' + rgbColor + ',0.15)';
  const solid = 'rgba(' + rgbColor + ',1)';
  new Chart(document.getElementById(canvasId), {{
    type: 'line',
    data: {{
      labels: times,
      datasets: [
        {{
          label: yLabel + ' max',
          data: maxData,
          borderColor: 'transparent',
          backgroundColor: alpha,
          fill: '+1',
          pointRadius: 0,
          order: 3,
        }},
        {{
          label: yLabel + ' mean',
          data: meanData,
          borderColor: solid,
          backgroundColor: 'transparent',
          borderWidth: 2,
          pointRadius: 4,
          fill: false,
          order: 1,
        }},
        {{
          label: yLabel + ' min',
          data: minData,
          borderColor: 'transparent',
          backgroundColor: alpha,
          fill: '-1',
          pointRadius: 0,
          order: 3,
        }},
      ],
    }},
    options: {{
      responsive: true,
      plugins: {{
        legend: {{
          labels: {{ filter: item => item.text.includes('mean') }},
        }},
        tooltip: {{
          callbacks: {{ title: ctx => 'Time: ' + ctx[0].label + ' s' }},
        }},
      }},
      scales: {{
        x: {{
          title: {{ display: true, text: 'Time (s)' }},
          type: 'linear',
          ticks: {{ callback: v => v + ' s' }},
        }},
        y: {{ title: {{ display: true, text: yLabel }} }},
      }},
    }},
  }});
}}

entries.forEach(e => {{
  const safeId = (e.id || 'run').replace(/\\s+/g, '_');
  const d = e.series;
  makeChart('tempChart_' + safeId, 'Temperature (hK)',
            d.times, d.temp_mean, d.temp_min, d.temp_max, '220,50,50');
  makeChart('concChart_' + safeId, 'Concentration',
            d.times, d.conc_mean, d.conc_min, d.conc_max, '50,100,220');
}});
</script>
</body>
</html>"""
    with open(out_path, "w") as f:
        f.write(html)
    log(f"Dashboard written: {out_path} ({len(named_series) * 2} charts)")


def evaluate_label_temperature(out_dir, nii_path, label_val=5, dt=60.0):
    """
    Log temperature and concentration statistics per time step for the label region,
    and return the collected time-series dict (or None if unavailable).
    """
    log(f"Mean temperature in label=={label_val} region:")
    return _collect_label_timeseries(out_dir, nii_path, label_val, dt)


def _fix_ownership(out_dir, reference_path):
    """When running as root inside Docker, chown out_dir to match the input file owner."""
    if os.getuid() != 0:
        return
    try:
        st = os.stat(reference_path)
        if st.st_uid == 0:
            return  # input also root-owned — nothing to fix
        for fname in os.listdir(out_dir):
            os.chown(os.path.join(out_dir, fname), st.st_uid, st.st_gid)
        log(f"Fixed output ownership → {st.st_uid}:{st.st_gid}")
    except Exception as e:
        log(f"Warning: could not fix output ownership: {e}")


def run_solver(out_dir, steps, dt, vesselcoupling,
               vtp1d_in=None, vtp1d_out=None):
    """
    Launch thermoembo1d once.  Phase field is initialised directly from the
    vessel-label image (vessel=1, tissue=0); no binary pre-solve cache.
    Solution VTU snapshots: resultsolution000.NNNN.vtu
    """
    abs_out = os.path.abspath(out_dir)
    cmd     = _solver_cmd(abs_out, steps, dt, vesselcoupling,
                          vtp1d_in=vtp1d_in, vtp1d_out=vtp1d_out)

    log("Running thermoembo1d ...")
    log("  " + " ".join(cmd))
    r = subprocess.run(cmd)
    if r.returncode != 0:
        log(f"Warning: thermoembo1d exited with code {r.returncode}")
    else:
        log("thermoembo1d completed successfully")
        log(f"Solution VTU files written to: {out_dir}/resultsolution*.vtu")
        _rewrite_vtu_files(out_dir, "result*.vtu")


# ── Phantom-dir mode: process all 00?_vessel_phantom.nii.gz ──────────────────

def process_phantom(sample_id, phantom_dir, out_dir, args):
    """
    Full pipeline for one phantom sample.
    phantom_dir: directory containing {id}_vessel_phantom.nii.gz and {id}_seed.fcsv
    """
    phantom_nii = os.path.join(phantom_dir, f"{sample_id}_vessel_phantom.nii.gz")
    seed_fcsv   = os.path.join(phantom_dir, f"{sample_id}_seed.fcsv")
    if not os.path.exists(phantom_nii):
        log(f"  {sample_id}: {phantom_nii} not found — skipping")
        return

    sample_out = os.path.join(out_dir, sample_id)
    os.makedirs(sample_out, exist_ok=True)

    def p(name): return os.path.join(sample_out, name)

    log(f"{'='*60}")
    log(f"Processing sample {sample_id}")

    # ── 1a. Inflow skeleton (label==2) ───────────────────────────────────────
    inflow_skel_nii = p(f"{sample_id}_inflow_skel.nii.gz")
    inflow_raw_vtp  = p(f"{sample_id}_inflow_raw.vtp")
    inflow_vtp      = p(f"{sample_id}_inflow_centerline.vtp")
    if make_label_skeleton(phantom_nii, inflow_skel_nii, label_val=2):
        log("Extracting inflow centerline ...")
        skel_to_vtp(inflow_skel_nii, inflow_raw_vtp)
        seed_mm = (resistance_lumping.read_fcsv_seed(seed_fcsv)
                   if os.path.exists(seed_fcsv) else None)
        if seed_mm is None:
            log("  No seed FCSV found — omitting injection BC")
        log("Running resistance lumping on inflow ...")
        resistance_lumping.solve(
            vtp_in=inflow_raw_vtp,
            vtp_out=inflow_vtp,
            opts={
                'label_path':    phantom_nii,
                'label_val':     2,
                'p_in_mmhg':     836.0,    # 1.1 atm — HA root (max-radius leaf)
                'p_out_mmhg':    760.0,    # 1.0 atm — terminal branches
                'p_seed_mmhg':   882.0,    # 1.16 atm — injection site
                'seed_coord_mm': seed_mm,
            })
    else:
        inflow_vtp = None

    # ── 1b. Outflow skeleton (label==3) ─────────────────────────────────────
    outflow_skel_nii = p(f"{sample_id}_outflow_skel.nii.gz")
    outflow_raw_vtp  = p(f"{sample_id}_outflow_raw.vtp")
    outflow_vtp      = p(f"{sample_id}_outflow_centerline.vtp")
    if make_label_skeleton(phantom_nii, outflow_skel_nii, label_val=3):
        log("Extracting outflow centerline ...")
        skel_to_vtp(outflow_skel_nii, outflow_raw_vtp)
        log("Running resistance lumping on outflow ...")
        resistance_lumping.solve(
            vtp_in=outflow_raw_vtp,
            vtp_out=outflow_vtp,
            opts={
                'label_path':  phantom_nii,
                'label_val':   3,
                'p_in_mmhg':   752.0,    # 0.99 atm — IVC end (max-radius leaf)
                'p_out_mmhg':  760.0,    # 1.0 atm — tissue terminal branches
            })
    else:
        outflow_vtp = None

    # ── 2. 3-D liver tet mesh ─────────────────────────────────────────────────
    mesh_vtu = p("mesh.vtu")
    mesh_exo = p("mesh.exo")
    vessel_vtk = p("vessel.vtk")
    log("Generating 3-D liver tet mesh ...")
    data, affine = load_nii(phantom_nii)
    surface_raw  = extract_surface(data, affine)
    surface      = process_surface(surface_raw, decimate=args.decimate)
    vol          = generate_tet_mesh(surface,
                                     min_dihedral=args.tet_dihedral,
                                     max_vol=args.maxvol)
    write_outputs(surface, vol, sample_out)

    if not os.path.exists(mesh_vtu):
        log(f"  Warning: meshing failed for {sample_id} — skipping solver")
        return

    # ── 3. VTU → Exodus ───────────────────────────────────────────────────────
    vtu_to_exodus(mesh_vtu, mesh_exo)

    # ── 4. Material VTK image ─────────────────────────────────────────────────
    make_vessel_vtk(phantom_nii, vessel_vtk)

    log(f"Mesh preparation complete for {sample_id}:")
    for f in [f for f in [inflow_vtp, outflow_vtp, mesh_vtu, mesh_exo, vessel_vtk] if f and os.path.exists(f)]:
        size = os.path.getsize(f) / 1024
        log(f"  {os.path.basename(f):35s}  {size:8.1f} kB")

    if args.mesh_only:
        log("--mesh-only specified; skipping simulation run.")
        return

    # ── 5. Run thermoembo1d ───────────────────────────────────────────────────
    run_solver(sample_out, args.steps, args.dt, args.vesselcoupling,
               vtp1d_in=inflow_vtp, vtp1d_out=outflow_vtp)

    # ── 6. Collect tumor temperature/concentration time-series ────────────────
    series = evaluate_label_temperature(sample_out, phantom_nii, label_val=4, dt=args.dt)
    return series


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Prepare meshes and run thermoembo1d",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    # Single-file mode (original)
    ap.add_argument("nii", nargs='?', default=None,
                    help="Labeled NIfTI (.nii or .nii.gz) with label 1=liver, 2=vessel")

    # Phantom-dir mode (new)
    ap.add_argument("--phantom-dir", default=None, metavar="DIR",
                    help="Directory containing 00?_vessel_phantom.nii.gz and 00?_seed.fcsv")
    ap.add_argument("--phantom-ids", default=None, metavar="IDS",
                    help="Comma-separated IDs to process (default: all found in --phantom-dir)")

    ap.add_argument("--out-dir",       default="thermoembo_run",  metavar="DIR",
                    help="Output directory (default: thermoembo_run)")
    ap.add_argument("--decimate",      type=float, default=0.90,  metavar="R",
                    help="Surface decimation reduction ratio (default 0.90 = keep 10%%)")
    ap.add_argument("--maxvol",        type=float, default=5000,  metavar="MM3",
                    help="TetGen max tet volume mm³ (default 5000)")
    ap.add_argument("--tet-dihedral",  type=float, default=15.0,  metavar="DEG",
                    help="TetGen min dihedral angle (default 15)")
    ap.add_argument("--steps",         type=int,   default=5,     metavar="N",
                    help="thermoembo1d time steps (default 5)")
    ap.add_argument("--dt",            type=float, default=60.0,  metavar="S",
                    help="thermoembo1d time-step size in seconds (default 60)")
    ap.add_argument("--vesselcoupling",type=float, default=1.0e-3,metavar="G",
                    help="Vessel-tissue coupling conductance [1/s/atm] (default 1e-3)")
    ap.add_argument("--mesh-only",     action="store_true",
                    help="Prepare meshes only — skip simulation run")

    args = ap.parse_args()

    # ── Phantom-dir mode ────────────────────────────────────────────────────
    if args.phantom_dir is not None:
        phantom_dir = os.path.abspath(args.phantom_dir)
        out_dir     = os.path.abspath(args.out_dir)
        os.makedirs(out_dir, exist_ok=True)

        if args.phantom_ids:
            ids = [s.strip() for s in args.phantom_ids.split(',')]
        else:
            found = glob.glob(os.path.join(phantom_dir, '*_vessel_phantom.nii.gz'))
            ids   = sorted(os.path.basename(f).split('_vessel')[0] for f in found)
        if not ids:
            sys.exit(f"Error: no *_vessel_phantom.nii.gz found in {phantom_dir}")
        log(f"Phantom IDs to process: {ids}")
        all_series = []
        for sid in ids:
            series = process_phantom(sid, phantom_dir, out_dir, args)
            if series:
                all_series.append((sid, series))
        if all_series:
            _write_dashboard(os.path.join(out_dir, "dashboard.html"), all_series, label_val=4)
        return

    # ── Single-file mode (original) ─────────────────────────────────────────
    if args.nii is None:
        ap.error("Provide either a NIfTI file or --phantom-dir")

    nii_path = os.path.abspath(args.nii)
    out_dir  = os.path.abspath(args.out_dir)

    if not os.path.exists(nii_path):
        sys.exit(f"Error: {nii_path} not found")

    os.makedirs(out_dir, exist_ok=True)

    # Paths inside out_dir
    skel_nii       = os.path.join(out_dir, "vessel_skeleton.nii.gz")
    centerline_vtp = os.path.join(out_dir, "centerline.vtp")
    mesh_vtu       = os.path.join(out_dir, "mesh.vtu")
    mesh_exo       = os.path.join(out_dir, "mesh.exo")
    vessel_vtk     = os.path.join(out_dir, "vessel.vtk")

    # ── 1. Skeletonise vessel ────────────────────────────────────────────────
    make_label_skeleton(nii_path, skel_nii, label_val=2)

    # ── 2. 1-D centerline mesh ───────────────────────────────────────────────
    log("Extracting centerline ...")
    skel_to_vtp(skel_nii, centerline_vtp)

    # ── 3. 3-D liver tet mesh ────────────────────────────────────────────────
    log("Generating 3-D liver tet mesh ...")
    data, affine = load_nii(nii_path)
    surface_raw  = extract_surface(data, affine)
    surface      = process_surface(surface_raw, decimate=args.decimate)
    vol          = generate_tet_mesh(surface,
                                     min_dihedral=args.tet_dihedral,
                                     max_vol=args.maxvol)
    write_outputs(surface, vol, out_dir)

    if not os.path.exists(mesh_vtu):
        sys.exit(f"Error: expected {mesh_vtu} — meshing may have failed")

    # ── 4. VTU → Exodus ─────────────────────────────────────────────────────
    vtu_to_exodus(mesh_vtu, mesh_exo)

    # ── 5. Material VTK image ────────────────────────────────────────────────
    make_vessel_vtk(nii_path, vessel_vtk)

    log("=" * 60)
    log("Mesh preparation complete:")
    for f in [centerline_vtp, mesh_vtu, mesh_exo, vessel_vtk]:
        size = os.path.getsize(f) / 1024
        log(f"  {os.path.basename(f):25s}  {size:8.1f} kB")
    log("=" * 60)

    if args.mesh_only:
        log("--mesh-only specified; skipping simulation run.")
        return

    # ── 6. Run thermoembo1d ──────────────────────────────────────────────────
    run_solver(out_dir, args.steps, args.dt, args.vesselcoupling,
               vtp1d_in=centerline_vtp)

    # ── 7. Evaluate mean temperature/concentration in tumor region ──────────
    series = evaluate_label_temperature(out_dir, nii_path, label_val=4, dt=args.dt)
    if series:
        _write_dashboard(os.path.join(out_dir, "dashboard.html"), [("", series)], label_val=4)

    _fix_ownership(out_dir, nii_path)


if __name__ == "__main__":
    main()
