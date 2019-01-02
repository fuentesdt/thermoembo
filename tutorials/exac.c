// $ c3d segmentation.nii.gz -info
// Image #1: dim = [512, 512, 37];  bb = {[-224.8 -200 -420], [175.2 200 -235]};  vox = [0.78125, 0.78125, 5];  range = [0, 5];  orient = RAI
// c3d -verbose segmentation.nii.gz -region 0x0x15vox 512x512x17vox -type uchar -o liver.vtk -replace 3 1 2 1 5 0  -canny 1mm 0 2 -o canny.vtk -as A -dilate 1 1x1x0 -push A  -scale -1 -add -o cannyedge.vtk
// ./exac -dim 3 -simplex -temp_petscspace_degree 1 -dm_view -ts_type beuler -ts_max_steps 20 -ts_dt 1.e0 -pc_type bjacobi -ksp_monitor_short -ksp_rtol 1.e-12 -ksp_converged_reason -snes_type ksponly -snes_monitor_short -snes_lag_jacobian 1  -snes_converged_reason -ts_monitor  -vtk liver.vtk -edge cannyedge.vtk -log_summary 
static char help[] = "Heat Equation in 2d and 3d with finite elements.\n\
We solve the heat equation in a rectangular\n\
domain, using a parallel unstructured mesh (DMPLEX) to discretize it.\n\
Contributed by: Julian Andrej <juan@tf.uni-kiel.de>\n\n\n";

#include <petscdmplex.h>
#include <petscds.h>
#include <petscts.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkStructuredPoints.h>
#include <vtkDataSetReader.h>

/*
  Heat equation:

    du/dt - \Delta u = -1 * dim

  Exact 2D solution:

    u = 2t + x^2 + y^2

    2 - (2 + 2) + 2 = 0

  Exact 3D solution:

    u = 3t + x^2 + y^2 + z^2

    3 - (2 + 2 + 2) + 3 = 0
*/

typedef struct {
  PetscInt          dim;
  PetscBool         simplex;
  char          imagefile[2048];   /* The vtk Image file */
  char           edgefile[2048];   /* The vtk Image file */
  PetscErrorCode (**exactFuncs)(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx);
  vtkSmartPointer<vtkImageData> ImageData ; 
  vtkSmartPointer<vtkImageData> EdgeData ; 
  double bounds[6];
} AppCtx;

// FIXME - need interface update
AppCtx         *_global_HACK_ctx;

static PetscErrorCode analytic_temp(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  PetscInt d;

  *u = dim*time;
  for (d = 0; d < dim; ++d) *u += x[d]*x[d];
  return 0;
}

static void f0_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  f0[0] = u_t[0] + (PetscScalar) dim;
}

static void f1_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) {
    f1[d] = u_x[d];
  }
}

static void g3_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) {
    g3[d*dim+d] = 1.0;
  }
}

static void g0_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  g0[0] = u_tShift*1.0;
}

static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options)
{
  PetscErrorCode ierr;
  PetscBool      flg;

  PetscFunctionBeginUser;
  options->dim     = 2;
  options->simplex = PETSC_TRUE;

  ierr = PetscOptionsBegin(comm, "", "Heat Equation Options", "DMPLEX");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-dim", "The topological mesh dimension", "ex45.c", options->dim, &options->dim, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-simplex", "Simplicial (true) or tensor (false) mesh", "ex45.c", options->simplex, &options->simplex, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-vtk", "vtk material filename to read", "exac.c", options->imagefile, options->imagefile, sizeof(options->imagefile), &flg);CHKERRQ(ierr);
  if (flg)
     {
       ierr = PetscPrintf(PETSC_COMM_WORLD, "opening file...\n");CHKERRQ(ierr);
       vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
       reader->SetFileName(options->imagefile);
       reader->Update();
       // initilize the bounding data structure
       options->ImageData = vtkImageData::SafeDownCast( reader->GetOutput()) ;
       options->ImageData->PrintSelf(std::cout,vtkIndent());
       options->ImageData->GetBounds(options->bounds);
       ierr = PetscPrintf(PETSC_COMM_WORLD, "ZBounds [%10.3e,%10.3e] \n",
                            options->bounds[4],options->bounds[5]);
     }
  else 
     {
       options->ImageData = 0;
     }
  ierr = PetscOptionsString("-edge", "vtk edge enhancement file to read", "exac.c", options->edgefile, options->edgefile, sizeof(options->edgefile), &flg);CHKERRQ(ierr);
  if (flg)
     {
       ierr = PetscPrintf(PETSC_COMM_WORLD, "opening file...\n");CHKERRQ(ierr);
       vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
       reader->SetFileName(options->edgefile);
       reader->Update();
       // initilize the bounding data structure
       options->EdgeData = vtkImageData::SafeDownCast( reader->GetOutput()) ;
       options->EdgeData->PrintSelf(std::cout,vtkIndent());
       // FIXME - add error checking that the bounds are the same
       options->EdgeData->GetBounds(options->bounds);
       ierr = PetscPrintf(PETSC_COMM_WORLD, "ZBounds [%10.3e,%10.3e] \n",
                            options->bounds[4],options->bounds[5]);
     }
  else 
     {
       options->EdgeData = 0;
     }
  ierr = PetscOptionsEnd();
  PetscFunctionReturn(0);
}

static PetscErrorCode CreateBCLabel(DM dm, const char name[])
{
  DMLabel        label;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = DMCreateLabel(dm, name);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, name, &label);CHKERRQ(ierr);
  ierr = DMPlexMarkBoundaryFaces(dm, 1, label);CHKERRQ(ierr);
  ierr = DMPlexLabelComplete(dm, label);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode edgerefinement(const PetscReal xloc[], PetscReal *limit)
{
  PetscErrorCode ierr;
  // FIXME - interface
  AppCtx *user= _global_HACK_ctx;

  double coord[3]= {xloc[0],xloc[1],xloc[2]};
  double pcoord[3];
  int    index[3];
  //transform the point and return the intensity value
  *limit = 50.;
  if ( user->EdgeData->ComputeStructuredCoordinates(coord,index,pcoord) )
   {
     // get material property
     PetscInt edgeid = static_cast<PetscInt>( user->EdgeData->GetScalarComponentAsDouble(index[0],index[1],index[2],0) );
     ierr = PetscPrintf(PETSC_COMM_WORLD, "(%10.3e,%10.3e,%10.3e) %d \n",coord[0],coord[1],coord[2],edgeid );CHKERRQ(ierr);
     // refine at the edges
     if (edgeid == 1 )
      {
        *limit = 5.;
      }
   }
   
  PetscFunctionReturn(0);
}

static PetscErrorCode CreateMesh(MPI_Comm comm, DM *dm, AppCtx *ctx)
{
  DM              pdm = NULL;
  const PetscInt  dim = ctx->dim;
  const PetscReal lower[3]= {ctx->bounds[0],ctx->bounds[2],ctx->bounds[4]};
  const PetscReal upper[3]= {ctx->bounds[1],ctx->bounds[3],ctx->bounds[5]};
  PetscBool       hasLabel;
  DM              refinedm = NULL;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  ierr = DMPlexCreateBoxMesh(comm, dim, ctx->simplex, NULL, lower, upper, NULL, PETSC_TRUE, dm);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) *dm, "Mesh");CHKERRQ(ierr);
  /* If no boundary marker exists, mark the whole boundary */
  ierr = DMHasLabel(*dm, "marker", &hasLabel);CHKERRQ(ierr);
  if (!hasLabel) {ierr = CreateBCLabel(*dm, "marker");CHKERRQ(ierr);}
  /* Distribute mesh over processes */
  ierr = DMPlexDistribute(*dm, 0, NULL, &pdm);CHKERRQ(ierr);
  if (pdm) {
    ierr = DMDestroy(dm);CHKERRQ(ierr);
    *dm  = pdm;
  }
  ierr = DMSetFromOptions(*dm);CHKERRQ(ierr);
  ierr = DMViewFromOptions(*dm, NULL, "-dm_view");CHKERRQ(ierr);


  ierr =  DMPlexSetRefinementFunction(*dm, edgerefinement);CHKERRQ(ierr);
  ierr = DMRefine(*dm, PetscObjectComm((PetscObject) dm), &refinedm);CHKERRQ(ierr);
  if (refinedm) {
    ierr = DMDestroy(dm);CHKERRQ(ierr);
    *dm  = refinedm;
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode SetupProblem(PetscDS prob, AppCtx *ctx)
{
  const PetscInt id = 1;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscDSSetResidual(prob, 0, f0_temp, f1_temp);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(prob, 0, 0, g0_temp, NULL, NULL, g3_temp);CHKERRQ(ierr);
  ctx->exactFuncs[0] = analytic_temp;
  ierr = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, "wall", "marker", 0, 0, NULL, (void (*)(void)) ctx->exactFuncs[0], 1, &id, ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SetupDiscretization(DM dm, AppCtx* ctx)
{
  DM             cdm = dm;
  const PetscInt dim = ctx->dim;
  PetscDS        prob;
  PetscFE        fe;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* Create finite element */
  ierr = PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, ctx->simplex, "temp_", -1, &fe);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fe, "temperature");CHKERRQ(ierr);
  /* Set discretization and boundary conditions for each mesh */
  ierr = DMGetDS(dm, &prob);CHKERRQ(ierr);
  ierr = PetscDSSetDiscretization(prob, 0, (PetscObject) fe);CHKERRQ(ierr);
  ierr = SetupProblem(prob, ctx);CHKERRQ(ierr);
  while (cdm) {
    PetscBool hasLabel;

    ierr = DMSetDS(cdm, prob);CHKERRQ(ierr);
    ierr = DMHasLabel(cdm, "marker", &hasLabel);CHKERRQ(ierr);
    if (!hasLabel) {ierr = CreateBCLabel(cdm, "marker");CHKERRQ(ierr);}
    ierr = DMGetCoarseDM(cdm, &cdm);CHKERRQ(ierr);
  }
  ierr = PetscFEDestroy(&fe);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc, char **argv)
{
  AppCtx         ctx;
  _global_HACK_ctx =  &ctx; // FIXME  - interface
  DM             dm;
  TS             ts;
  Vec            u, r;
  PetscReal      t       = 0.0;
  PetscReal      L2error = 0.0;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, help);if (ierr) return ierr;
  ierr = ProcessOptions(PETSC_COMM_WORLD, &ctx);CHKERRQ(ierr);
  ierr = CreateMesh(PETSC_COMM_WORLD, &dm, &ctx);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(dm, &ctx);CHKERRQ(ierr);
  ierr = PetscMalloc1(1, &ctx.exactFuncs);CHKERRQ(ierr);
  ierr = SetupDiscretization(dm, &ctx);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm, &u);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) u, "temperature");CHKERRQ(ierr);
  ierr = VecDuplicate(u, &r);CHKERRQ(ierr);

  ierr = TSCreate(PETSC_COMM_WORLD, &ts);CHKERRQ(ierr);
  ierr = TSSetDM(ts, dm);CHKERRQ(ierr);
  ierr = DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, &ctx);CHKERRQ(ierr);
  ierr = DMTSSetIFunctionLocal(dm, DMPlexTSComputeIFunctionFEM, &ctx);CHKERRQ(ierr);
  ierr = DMTSSetIJacobianLocal(dm, DMPlexTSComputeIJacobianFEM, &ctx);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  ierr = DMProjectFunction(dm, t, ctx.exactFuncs, NULL, INSERT_ALL_VALUES, u);CHKERRQ(ierr);
  //ierr = TSMonitorSet(ts,TSMonitorSolutionVTK,&ctx,(void*)&TSMonitorSolutionVTKDestroy);CHKERRQ(ierr);
  // write vtk file at every time point
  char               vtkfilenametemplate[PETSC_MAX_PATH_LEN] = "solution.%04d.vtu";
  ierr = TSMonitorSet(ts,TSMonitorSolutionVTK,&vtkfilenametemplate,NULL);CHKERRQ(ierr);
  // ierr = TSSetPostStage(ts,TSUpdateArrhenius);CHKERRQ(ierr);
  ierr = TSSolve(ts, u);CHKERRQ(ierr);

  ierr = TSGetTime(ts, &t);CHKERRQ(ierr);
  ierr = DMComputeL2Diff(dm, t, ctx.exactFuncs, NULL, u, &L2error);CHKERRQ(ierr);
  if (L2error < 1.0e-11) {ierr = PetscPrintf(PETSC_COMM_WORLD, "L_2 Error: < 1.0e-11\n");CHKERRQ(ierr);}
  else                   {ierr = PetscPrintf(PETSC_COMM_WORLD, "L_2 Error: %g\n", (double)L2error);CHKERRQ(ierr);}
  ierr = VecViewFromOptions(u, NULL, "-sol_vec_view");CHKERRQ(ierr);

  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscFree(ctx.exactFuncs);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/*TEST

  # Full solves
  test:
    suffix: 2d_p1_r1
    requires: triangle
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dm_refine 1 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_p1_r3
    requires: triangle
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dm_refine 3 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_p2_r1
    requires: triangle
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dm_refine 1 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_p2_r3
    requires: triangle
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dm_refine 3 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q1_r1
    requires: !single
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 1 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q1_r3
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 3 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q2_r1
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 1 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q2_r3
    requires: !single
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 3 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p1_r1
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 1 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p1_r2
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 2 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p2_r1
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 1 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p2_r2
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 2 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q1_r1
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 1 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q1_r2
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 2 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q2_r1
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 1 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q2_r2
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 2 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor

TEST*/
