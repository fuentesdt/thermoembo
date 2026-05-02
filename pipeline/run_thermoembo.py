"""
run_thermoembo.py — Prepare and launch a thermoembo1d coupled 1D-3D simulation.

Pipeline
--------
  1. Extract label==2 vessel mask → skeletonize → vessel_skeleton.nii.gz
  2. centerline.skel_to_vtp()  → centerline.vtp  (1-D vessel graph)
  3. meshing.*                 → mesh.vtu         (3-D liver tet mesh)
  4. Convert mesh.vtu → mesh.exo                  (mm → m, meshio/netCDF4)
  5. Write vessel.vtk  (legacy VTK image from NIfTI label, mm coords)
  6. Run thermoembo1d directly (two-pass: phase presolve + main solve)

Inputs
------
  <nii>  — labeled NIfTI (.nii or .nii.gz)
    label == 1 : liver parenchyma
    label == 2 : vessel lumen

Outputs (written to --out-dir, default: thermoembo_run/)
  vessel_skeleton.nii.gz   skeletonised label-2 volume
  centerline.vtp           1-D vessel centerline mesh
  surface.vtp              decimated liver surface
  mesh.vtu                 tetrahedral volume mesh  (mm)
  mesh.exo                 Exodus mesh              (m, for PETSc/thermoembo1d)
  vessel.vtk               legacy VTK image (mm, material field)
  result*.vtu              thermoembo1d solution snapshots

Usage
-----
  python run_thermoembo.py prepost/PreTxArtLoRes.vessellabel.nii.gz
  python run_thermoembo.py prepost/PreTxArtLoRes.vessellabel.nii.gz \\
      --out-dir thermoembo_run --decimate 0.90 --maxvol 5000 --steps 5

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


def evaluate_label_temperature(out_dir, nii_path, label_val=5):
    """
    For each resultsolution*.vtu, report mean temperature at mesh nodes whose
    corresponding NIfTI voxel has label==label_val.
    Mesh coordinates are in metres; NIfTI affine is in mm — multiply by 1e3.
    """
    vtu_files = sorted(glob.glob(os.path.join(out_dir, "resultsolution*.vtu")))
    if not vtu_files:
        return

    nii    = nib.load(nii_path)
    labels = np.asarray(nii.dataobj)
    if not (labels == label_val).any():
        log(f"Label=={label_val} not present in {os.path.basename(nii_path)}; "
            f"skipping temperature evaluation")
        return

    inv_aff = np.linalg.inv(nii.affine)        # affine is mm; coords must be mm
    shape   = np.array(labels.shape)

    log(f"Mean temperature in label=={label_val} region:")
    for vtu_path in vtu_files:
        m   = pv.read(vtu_path)
        pts = m.points * 1e3                   # m → mm
        ph  = np.hstack([pts, np.ones((len(pts), 1))])
        vox = (inv_aff @ ph.T).T[:, :3]
        idx = np.clip(np.round(vox).astype(int), 0, shape - 1)
        mask = labels[idx[:, 0], idx[:, 1], idx[:, 2]] == label_val
        if not mask.any():
            log(f"  {os.path.basename(vtu_path)}: no label=={label_val} nodes in mesh")
            continue
        temp = m.point_data["solutiontemperature.0"][mask]
        step = os.path.basename(vtu_path).replace("resultsolution", "step").replace(".vtu", "")
        log(f"  {step}: n={mask.sum():5d}  mean={np.mean(temp):.4f}  "
            f"min={np.min(temp):.4f}  max={np.max(temp):.4f}")


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
    ap.add_argument("--vesselcoupling",type=float, default=1.0e-4,metavar="G",
                    help="Vessel-tissue coupling conductance [1/s/atm] (default 1e-4)")
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
        for sid in ids:
            process_phantom(sid, phantom_dir, out_dir, args)
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

    # ── 7. Evaluate mean temperature in tumor region ─────────────────────────
    evaluate_label_temperature(out_dir, nii_path, label_val=4)

    _fix_ownership(out_dir, nii_path)


if __name__ == "__main__":
    main()
