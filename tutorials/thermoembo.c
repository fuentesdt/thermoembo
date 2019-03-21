// ./thermoembo -dim 3 -temp_petscspace_degree 1 -pres_petscspace_degree 1 -damg_petscspace_degree 1 -conc_petscspace_degree 1 -phas_petscspace_degree 1 -dm_view -ts_type beuler -pc_type fieldsplit  -ksp_monitor_short -ksp_type preonly -ksp_converged_reason -snes_type newtonls  -snes_rtol 9.e-1 -snes_monitor_short -snes_lag_jacobian 1  -snes_converged_reason -ts_monitor -log_summary -artdiff 1e-6  -ts_max_steps 40 -ts_dt 1.e-1  -snes_linesearch_monitor -info -info_exclude  null,vec,mat,pc   -pc_fieldsplit_type additive  -fieldsplit_u_pc_type bjacobi  -fieldsplit_u_ksp_converged_reason -fieldsplit_u_ksp_monitor_short -fieldsplit_u_ksp_type gmres -fieldsplit_u_ksp_rtol 1.e-4  -fieldsplit_s_pc_type bjacobi -fieldsplit_s_ksp_rtol 1.e-9 -fieldsplit_s_ksp_converged_reason -fieldsplit_s_ksp_monitor_short -fieldsplit_s_ksp_type gmres  -salttemp .57  -phasepresolve_pc_type fieldsplit -phasepresolve_ksp_type preonly  -phasepresolve_ts_type beuler -phasepresolve_ts_max_steps 20 -phasepresolve_fieldsplit_c_pc_type bjacobi -phasepresolve_fieldsplit_c_ksp_type gmres -phasepresolve_fieldsplit_1_ksp_type preonly -phasepresolve_ksp_monitor_short -phasepresolve_fieldsplit_c_ksp_monitor_short -phasepresolve_fieldsplit_1_ksp_monitor_short -phasepresolve_fieldsplit_c_ksp_rtol 1.e-12 -phasepresolve_fieldsplit_1_pc_type none -phasepresolve_ksp_converged_reason -phasepresolve_snes_type ksponly -phasepresolve_snes_monitor_short -phasepresolve_snes_lag_jacobian 1  -phasepresolve_snes_converged_reason -phasepresolve_ksp_view -phasepresolve_ts_monitor   -phasepresolve_pc_fieldsplit_type additive -vtk ../temperaturedata/Kidney1Left_04202017_Exp42/vesselregion.vtk  -log_summary  -dm_refine 2 -o test

// PCApply_FieldSplit 
// -snes_type <newtonls>: Nonlinear solver method (one of) newtonls newtontr test nrichardson ksponly vinewtonrsls vinewtonssls ngmres qn shell ngs ncg fas ms nasm anderson aspin composite (SNESSetType)
// SNESSolve_KSPONLY
// SNESSolve_NEWTONLS SNESLineSearchApply_BT  SNESLineSearchApply_CP
// SNESKSPSetUseEW
// SNESComputeJacobian

// ./thermoembo -dim 3 -temp_petscspace_degree 1 -pres_petscspace_degree 1 -damg_petscspace_degree 1 -conc_petscspace_degree 1 -ts_type beuler -ts_max_steps 20 -ts_dt 1.e0 -pc_type bjacobi -ksp_monitor -ksp_rtol 1.e-12 -ksp_converged_reason -snes_type ksponly -snes_monitor_short -snes_lag_jacobian 1  -snes_converged_reason -ts_monitor -log_summary -artdiff 1e-7 -velocity .11  -ksp_view_mat ascii:mat.m:ascii_matlab   -ksp_view_rhs ascii:rhs.m:ascii_matlab

static char help[] = "Heat Equation in 2d and 3d with finite elements.\n\
We solve the heat equation with an convection term using implicit explicit time stepping\n\
, using a parallel unstructured mesh (DMPLEX) to discretize it.\n\
Contributed by: Julian Andrej <juan@tf.uni-kiel.de>\n\n\n";

#include <petscdmplex.h>
#include <petscds.h>
#include <petscts.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkStructuredPoints.h>
#include <vtkDataSetReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>

/*
  parabolic equation:

    du/dt - \alpha \Delta u = - c * \mu_a u

*/

typedef enum {COEFF_NONE, COEFF_ANALYTIC, COEFF_FIELD, COEFF_NONLINEAR} CoeffType;
// use enum to identify fields and field deriviatives
typedef enum {FIELD_PHASE,FIELD_TEMPERATURE, FIELD_DAMAGE, FIELD_PRESSURE, FIELD_SATURATION} FieldEnumType;
const     int NUMPARAMETERS=25;
const     double _globalepsilon = 1.e-12;
// list of parameters
typedef enum {PARAM_OMEGA,
              PARAM_RHOBLOOD,
              PARAM_RHOOIL,
              PARAM_RHODCACL,
              PARAM_ALPHA,
              PARAM_SPECIFICHEATBLOOD,
              PARAM_SPECIFICHEATTISSUE,
              PARAM_USALT,
              PARAM_UARTERY,
              PARAM_POROSITY,
              PARAM_MOBILITYOIL,
              PARAM_MOBILITYBLOOD,
              PARAM_INJECTIONVELOCITY,
              PARAM_DISPLACEMENTPRESSURE,
              PARAM_BASELINEPRESSURE,
              PARAM_BOUNDARYPRESSURE,
              PARAM_EPSILON,
              PARAM_PHASEEPSILON,
              PARAM_SATURATION_SOURCE,
              PARAM_PRESSURE_SOURCE, 
              PARAM_TEMPERATURE_SOURCE, 
              PARAM_ADVECTIONTERM, 
              PARAM_ARTIFICIALDIFFUSION} ParameterType;
typedef struct {
  PetscInt          dim;
  PetscInt          refine;
  PetscReal         time_step; /* phase field timestep */
  PetscInt          max_steps; /* phase field max steps */
  PetscBool         simplex;
  PetscBool      fieldBC;
  char              imagefile[2048];   /* The vtk Image file */
  char              filenosuffix[2048] ;
  double            parameters[NUMPARAMETERS] ; //{param1, param2, ...}
  double            temperaturescaling ; //normalize temperature to [0,1]
  double bounds[6];
  double spacing[3];
  vtkSmartPointer<vtkImageData> ImageData ; 
  PetscInt numFields;
  char  **fieldNames;
  IS         *fields;
  IS   isnotpressuresaturation;
  CoeffType      variableCoefficient;
  PetscErrorCode (**exactFuncs)(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx);
} AppCtx;

PetscErrorCode TSUpdateArrhenius(TS ts, PetscReal stagetime, PetscInt stageindex, Vec *Y)
{
  AppCtx         *ctx;
  DM dm;
  Vec            temperature,saturation,work,damage;
  PetscReal deltaT; 
  PetscErrorCode ierr;

  ierr = TSGetDM(ts, &dm);CHKERRQ(ierr);
  ierr = TSGetTimeStep(ts, &deltaT);CHKERRQ(ierr);
  ierr = DMGetApplicationContext(dm, &ctx);CHKERRQ(ierr);
  double frequencyfactor=3.1e98; // 1/s
  double activationenergy=6.28e5; // J/mol
  double gasconstant=8.314; // J/mol/K
  double kelvinconversion=273; // J/mol/K

  ierr = VecGetSubVector(*Y, ctx->fields[FIELD_TEMPERATURE], &temperature);CHKERRQ(ierr);
  ierr = VecDuplicate(temperature,&work);CHKERRQ(ierr);
  ierr = VecGetSubVector(*Y, ctx->fields[FIELD_DAMAGE]     , &damage);CHKERRQ(ierr);

  // compute arrhenius update
  ierr = VecCopy(temperature,work);CHKERRQ(ierr);
  ierr = VecShift(work,kelvinconversion);CHKERRQ(ierr);
  ierr = VecReciprocal(work);CHKERRQ(ierr);
  ierr = VecScale(work,-activationenergy/gasconstant);CHKERRQ(ierr);
  ierr = VecExp(work);CHKERRQ(ierr);
  ierr = VecAXPY(damage,deltaT * frequencyfactor,work);CHKERRQ(ierr);

  // post process saturations - bound between 0 and 1
  ierr = VecGetSubVector(*Y, ctx->fields[FIELD_SATURATION], &saturation);CHKERRQ(ierr);
  int iii,nlocal;
  PetscReal    *array;
  ierr = VecGetLocalSize(saturation,&nlocal);
  ierr = VecGetArray(saturation,&array);
  for (iii=0; iii<nlocal; iii++) array[iii] = PetscMax(PetscMin(1, array[iii]),0);
  ierr = VecRestoreArray(saturation,&array);


  // cleanup
  ierr = VecRestoreSubVector(*Y, ctx->fields[FIELD_TEMPERATURE], &temperature);CHKERRQ(ierr);
  ierr = VecRestoreSubVector(*Y, ctx->fields[FIELD_DAMAGE]     , &damage     );CHKERRQ(ierr);
  ierr = VecRestoreSubVector(*Y, ctx->fields[FIELD_SATURATION] , &saturation );CHKERRQ(ierr);
  ierr = VecDestroy( &work);CHKERRQ(ierr);

  return 0;
}

static PetscErrorCode quadratic_u_2d(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  *u = x[0]*x[0] + x[1]*x[1];
  return 0;
}

static PetscErrorCode quadratic_u_3d(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  *u = 2.0*(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])/3.0;
  return 0;
}

static PetscErrorCode nu_2d(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  *u = 10.*PetscAbs(sin(x[0] + x[1]))+.01;
  const PetscScalar radius = .004;
  const PetscScalar    centroid[3]  = {0.0,0.003,0.0};
  if (
  (x[0] -  centroid[0])* (x[0] -  centroid[0])  + (x[1] -  centroid[1])* (x[1] -  centroid[1]) + (x[2] -  centroid[2])* (x[2] -  centroid[2])  < radius * radius
     ) { 
       *u = 1.e2;
       } 
  return 0;
}

static PetscErrorCode analytic_temp(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  PetscInt d;

  AppCtx *options = (AppCtx *)ctx;
  *u = options->parameters[PARAM_UARTERY];
  //*u = dim*time;
  //for (d = 0; d < dim; ++d) *u += x[d]*x[d];
  return 0;
}

static PetscErrorCode analytic_conc(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  PetscInt d;

  *u = 0.0  ;
  //*u = dim*time;
  //for (d = 0; d < dim; ++d) *u += x[d]*x[d];
  return 0;
}

static PetscErrorCode analytic_phas(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  PetscInt d; 
  PetscReal imagevalue;
  AppCtx *user = (AppCtx *)ctx;

  double coord[3]= {x[0],x[1],x[2]};
  double pcoord[3];
  int    index[3];
  if ( user->ImageData->ComputeStructuredCoordinates(coord,index,pcoord) )
   {
     // get material property
     imagevalue = static_cast<PetscReal>( user->ImageData->GetScalarComponentAsDouble(index[0],index[1],index[2],0) );
   }
  *u = imagevalue;
  
  return 0;
}

static PetscErrorCode analytic_pres(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  PetscInt d;

  AppCtx *options = (AppCtx *)ctx;
  *u = options->parameters[PARAM_BASELINEPRESSURE];
  // PetscReal  centroid[3] = {.09,.08,.01};
  // PetscReal  radiussq =    (x[0]-centroid[0])*(x[0]-centroid[0]) 
  //                         +(x[1]-centroid[1])*(x[1]-centroid[1]) 
  //                         +(x[2]-centroid[2])*(x[2]-centroid[2]) ;
  // PetscReal  variance =   .01*.01;
  // *u = 1.e-2/sqrt(2.0*3.141592653589793*variance )*exp(-radiussq/2.0/variance );
  // *u = dim*time;
  // for (d = 0; d < dim; ++d) *u += 1.e2*x[d]*x[d];
  return 0;
}

static PetscErrorCode analytic_damg(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  PetscInt d;

  *u = 0.0 ;
  //*u = dim*time;
  //for (d = 0; d < dim; ++d) *u += x[d]*x[d];
  return 0;
}

void zero(const PetscReal coords[], PetscScalar *u, void *ctx)
{
  *u = 0.0;
}
static PetscErrorCode zerotwo(PetscInt dim, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  PetscInt i;
  for (i = 0; i < dim; ++i) u[i] = -1.0;
  return 0;
}

static PetscErrorCode salttemperature(PetscInt dim, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  AppCtx *options = (AppCtx *)ctx;
  *u = options->parameters[PARAM_USALT];
  return 0;
}

static PetscErrorCode bodytemperature(PetscInt dim, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  AppCtx *options = (AppCtx *)ctx;
  *u = options->parameters[PARAM_UARTERY];
  return 0;
}


static PetscErrorCode bolusinjection(PetscInt dim, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  *u = 1.0;
  return 0;
}

static PetscErrorCode fieldzero(PetscInt dim, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  *u = 0.0;
  return 0;
}

static PetscErrorCode saltconcentration (PetscInt dim, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  *u = 1.0  ;
  return 0;
}


static PetscErrorCode tissuedamagefcn(PetscInt dim, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  *u = 1.0  ;
  return 0;
}
static void f0_p(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                 const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                 PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{ // break PetscFEIntegrateResidual_Basic
  f0[0]  = - constants[PARAM_PRESSURE_SOURCE]* u[FIELD_SATURATION] ;
}
static void f1_p(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                 const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                 PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{ // break PetscFEIntegrateResidual_Basic
  PetscInt d;
  double   totalmobility = constants[PARAM_MOBILITYOIL]+constants[PARAM_MOBILITYBLOOD];
  PetscReal  dpds    = -constants[PARAM_DISPLACEMENTPRESSURE]/2.0/(sqrt(1- PetscMin(1, u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) ;
  for (d = 0; d < dim; ++d) f1[d] = u_x[uOff_x[FIELD_PRESSURE]+d] * totalmobility +  constants[PARAM_MOBILITYOIL] * dpds * u_x[uOff_x[FIELD_SATURATION]+d];
}

static void g3_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{ // break PetscFEIntegrateJacobian_Basic
  PetscInt   d;
  double   totalmobility = constants[PARAM_MOBILITYOIL]+constants[PARAM_MOBILITYBLOOD];
  for (d = 0; d < dim; ++d) g3[d*dim+d] =  totalmobility;
}

static void g0_ps(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  g0[0] = - constants[PARAM_PRESSURE_SOURCE];
}

static void g2_ps(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{ // break PetscFEIntegrateJacobian_Basic
  PetscInt   d;
  PetscReal  dp2ds2  = -constants[PARAM_DISPLACEMENTPRESSURE]/4.0/(sqrt(1.- PetscMin(1., u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) /(sqrt(1.- PetscMin(1., u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) / sqrt(1.- PetscMin(1., u[FIELD_SATURATION])+_globalepsilon)*  ( (u[FIELD_SATURATION]<1.) ?  1. : 0.);;
  for (d = 0; d < dim; ++d) g2[d] = constants[PARAM_MOBILITYOIL] * dp2ds2 * u_x[uOff_x[FIELD_SATURATION]+d];
}


static void g3_ps(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{ // break PetscFEIntegrateJacobian_Basic
  PetscInt   d;
  PetscReal  dpds    = -constants[PARAM_DISPLACEMENTPRESSURE]/2.0/(sqrt(1- PetscMin(1, u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) ;
  for (d = 0; d < dim; ++d) g3[d*dim+d] = constants[PARAM_MOBILITYOIL] * dpds  ;
}


static PetscErrorCode bd_applicator_pres(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  AppCtx *options = (AppCtx *)ctx;
  *u = 1.e10 * PetscAbs(sin(3.141592653589793 * time/10.)); // 
  return 0;
}

static PetscErrorCode vessel_pres(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  AppCtx *options = (AppCtx *)ctx;
  *u = options->parameters[PARAM_BOUNDARYPRESSURE];
  return 0;
}

static PetscErrorCode baseline_pres(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  AppCtx *options = (AppCtx *)ctx;
  *u = options->parameters[PARAM_BASELINEPRESSURE];
  return 0;
}


//  PetscFEIntegrateBdResidual_Basic DMPlexComputeBdResidual_Internal
static void f0_bd_p(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  f0[0] = - constants[PARAM_INJECTIONVELOCITY];
}

void g0_bd_uu_3d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  PetscInt d;
  double radius = 0.0;
  for (d = 0; d < dim; ++d) radius += x[d]*x[d];
  radius = sqrt(radius);
  if ( radius > .02 ) // FIXME - HACK location dependent bc .... not invariant ... assign by vertex set ? 
    {g0[0] =   constants[2] ;}
}


static void f1_bd_zero(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscInt comp;
  for (comp = 0; comp < dim; ++comp) f1[comp] = 0.0;
}


static void f0_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt   comp;
  PetscReal  dpds    = -constants[PARAM_DISPLACEMENTPRESSURE]/2.0/(sqrt(1- PetscMin(1, u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) ;
  // buffers
  PetscReal  tmpone  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]*dpds;
  PetscReal  tmptwo  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* ((1-u[FIELD_SATURATION])*constants[PARAM_MOBILITYBLOOD] + u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]); 
  PetscReal  beta[3] = {
    tmpone*u_x[uOff_x[FIELD_SATURATION]+0] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+0] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+1] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+1] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+2] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+2]  };
  
  //PetscPrintf(PETSC_COMM_WORLD, "f0: u_t = %12.5e beta = %12.5e %12.5e %12.5e   ",u_t[FIELD_TEMPERATURE],beta[0], beta[1], beta[2] );
  double advection=0.0;
  for (comp = 0; comp < dim; ++comp) advection += u_x[uOff_x[FIELD_TEMPERATURE]+ comp] * beta[comp];
  f0[0] = u_t[FIELD_TEMPERATURE] + advection  -  constants[PARAM_TEMPERATURE_SOURCE]*u[FIELD_SATURATION];

}

static void g0_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  g0[0] = u_tShift*1.0 ;
}


static void f1_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscInt d;
  PetscReal  dpds    = -constants[PARAM_DISPLACEMENTPRESSURE]/2.0/(sqrt(1- PetscMin(1, u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) ;
  // buffers
  PetscReal  tmpone  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]*dpds;
  PetscReal  tmptwo  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* ((1-u[FIELD_SATURATION])*constants[PARAM_MOBILITYBLOOD] + u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]); 
  PetscReal  beta[3] = {
    tmpone*u_x[uOff_x[FIELD_SATURATION]+0] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+0] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+1] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+1] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+2] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+2]  };

  //PetscPrintf(PETSC_COMM_WORLD, "f1: u_x  %12.5e %12.5e %12.5e \n",u_x[uOff_x[FIELD_TEMPERATURE]+0] ,u_x[uOff_x[FIELD_TEMPERATURE]+1] ,u_x[uOff_x[FIELD_TEMPERATURE]+2] );
  double  innerprod = 0.0;
  for (d = 0; d < dim; ++d)  innerprod = innerprod + u_x[uOff_x[FIELD_TEMPERATURE]+d]  * beta[d];
  for (d = 0; d < dim; ++d) {
    f1[d] = constants[PARAM_ALPHA] * u_x[uOff_x[FIELD_TEMPERATURE]+d] 
          + constants[PARAM_ARTIFICIALDIFFUSION] * innerprod * beta[d];
  }
}

static void g1_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{
  PetscInt d;
  PetscReal  dpds    = -constants[PARAM_DISPLACEMENTPRESSURE]/2.0/(sqrt(1- PetscMin(1, u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) ;
  // buffers
  PetscReal  tmpone  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]*dpds;
  PetscReal  tmptwo  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* ((1-u[FIELD_SATURATION])*constants[PARAM_MOBILITYBLOOD] + u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]); 
  PetscReal  beta[3] = {
    tmpone*u_x[uOff_x[FIELD_SATURATION]+0] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+0] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+1] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+1] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+2] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+2]  };
  for (d = 0; d < dim; ++d) {
    g1[d] =  beta[d];
  }
}

static void g3_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscInt d,iii,jjj;
  PetscReal  dpds    = -constants[PARAM_DISPLACEMENTPRESSURE]/2.0/(sqrt(1- PetscMin(1, u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) ;
  // buffers
  PetscReal  tmpone  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]*dpds;
  PetscReal  tmptwo  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* ((1-u[FIELD_SATURATION])*constants[PARAM_MOBILITYBLOOD] + u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]); 
  PetscReal  beta[3] = {
    tmpone*u_x[uOff_x[FIELD_SATURATION]+0] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+0] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+1] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+1] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+2] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+2]  };

  //PetscPrintf(PETSC_COMM_WORLD, "%f ",conduction );
  for (d = 0; d < dim; ++d) {
    g3[d*dim+d] = constants[PARAM_ALPHA];
  }
  for (iii = 0; iii < dim; ++iii) for (jjj = 0; jjj < dim; ++jjj) {
    g3[iii*dim+jjj] = g3[iii*dim+jjj] + constants[PARAM_ARTIFICIALDIFFUSION] * beta[iii] * beta[jjj];
  }
}

static void g1_temppres(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{
  PetscInt d;
  PetscReal  dpds    = -constants[PARAM_DISPLACEMENTPRESSURE]/2.0/(sqrt(1- PetscMin(1, u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) ;
  // buffers
  PetscReal  tmptwo  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* ((1-u[FIELD_SATURATION])*constants[PARAM_MOBILITYBLOOD] + u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]); 

  for (d = 0; d < dim; ++d) {
    g1[d] = tmptwo *  u_x[uOff_x[FIELD_TEMPERATURE]+d];
  }

}

static void g3_temppres(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscInt d,iii,jjj;
  PetscReal  dpds    = -constants[PARAM_DISPLACEMENTPRESSURE]/2.0/(sqrt(1- PetscMin(1, u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) ;
  // buffers
  PetscReal  tmpone  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]*dpds;
  PetscReal  tmptwo  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* ((1-u[FIELD_SATURATION])*constants[PARAM_MOBILITYBLOOD] + u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]); 
  PetscReal  beta[3] = {
    tmpone*u_x[uOff_x[FIELD_SATURATION]+0] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+0] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+1] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+1] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+2] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+2]  };

  //PetscPrintf(PETSC_COMM_WORLD, "%f ",conduction );
  double  innerprod = 0.0;
  for (d = 0; d < dim; ++d)  innerprod = innerprod + u_x[uOff_x[FIELD_TEMPERATURE]+d]  * beta[d];

  for (iii = 0; iii < dim; ++iii) for (jjj = 0; jjj < dim; ++jjj) {
    g3[iii*dim+jjj] = u_x[uOff_x[FIELD_TEMPERATURE]+iii] * beta[jjj];
  }
  for (d = 0; d < dim; ++d) {
    g3[d*dim+d] = g3[d*dim+d] + innerprod;
  }
  for (iii = 0; iii < dim; ++iii) for (jjj = 0; jjj < dim; ++jjj) {
    g3[iii*dim+jjj] = g3[iii*dim+jjj] * tmptwo * constants[PARAM_ARTIFICIALDIFFUSION] ;
  }
}
static void g0_tempsat(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  PetscInt d;
  PetscReal  dpds    = -constants[PARAM_DISPLACEMENTPRESSURE]/2.0/(sqrt(1- PetscMin(1, u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) ;
  PetscReal  dp2ds2  = -constants[PARAM_DISPLACEMENTPRESSURE]/4.0/(sqrt(1.- PetscMin(1., u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) /(sqrt(1.- PetscMin(1., u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) / sqrt(1.- PetscMin(1., u[FIELD_SATURATION])+_globalepsilon)*  ( (u[FIELD_SATURATION]<1.) ?  1. : 0.);;
  PetscReal  dbtmpone=-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* (dpds + u[FIELD_SATURATION] * dp2ds2  ) *constants[PARAM_MOBILITYOIL];
  PetscReal  dbtmptwo=-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* (constants[PARAM_MOBILITYOIL]  - constants[PARAM_MOBILITYBLOOD]);
  PetscReal  dbetads[3] = {
      dbtmpone*u_x[uOff_x[FIELD_SATURATION]+0] + dbtmptwo*u_x[uOff_x[FIELD_PRESSURE]+0] ,
      dbtmpone*u_x[uOff_x[FIELD_SATURATION]+1] + dbtmptwo*u_x[uOff_x[FIELD_PRESSURE]+1] ,
      dbtmpone*u_x[uOff_x[FIELD_SATURATION]+2] + dbtmptwo*u_x[uOff_x[FIELD_PRESSURE]+2]  };
  double  innerprod = 0.0;
  for (d = 0; d < dim; ++d)  innerprod = innerprod + u_x[uOff_x[FIELD_TEMPERATURE]+d]  * dbetads[d];
  g0[0] = -  constants[PARAM_TEMPERATURE_SOURCE] + innerprod ;
}

static void g1_tempsat(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[])
{
  PetscInt d;
  PetscReal  dpds    = -constants[PARAM_DISPLACEMENTPRESSURE]/2.0/(sqrt(1- PetscMin(1, u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) ;
  PetscReal  tmpthree=-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]*u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]*dpds;

  for (d = 0; d < dim; ++d) {
    g1[d] =  tmpthree * u_x[uOff_x[FIELD_TEMPERATURE]+d];
  }
}

static void g2_tempsat(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[])
{
  PetscInt d,iii,jjj;


  PetscReal  dpds    = -constants[PARAM_DISPLACEMENTPRESSURE]/2.0/(sqrt(1- PetscMin(1, u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) ;
  PetscReal  dp2ds2  = -constants[PARAM_DISPLACEMENTPRESSURE]/4.0/(sqrt(1.- PetscMin(1., u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) /(sqrt(1.- PetscMin(1., u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) / sqrt(1.- PetscMin(1., u[FIELD_SATURATION])+_globalepsilon)*  ( (u[FIELD_SATURATION]<1.) ?  1. : 0.);;
  // buffers
  PetscReal  tmpone  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]*dpds;
  PetscReal  tmptwo  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* ((1-u[FIELD_SATURATION])*constants[PARAM_MOBILITYBLOOD] + u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]); 
  PetscReal  beta[3] = {
    tmpone*u_x[uOff_x[FIELD_SATURATION]+0] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+0] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+1] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+1] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+2] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+2]  };
  PetscReal  dbtmpone=-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* (dpds + u[FIELD_SATURATION] * dp2ds2  ) *constants[PARAM_MOBILITYOIL];
  PetscReal  dbtmptwo=-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* (constants[PARAM_MOBILITYOIL]  - constants[PARAM_MOBILITYBLOOD]);
  PetscReal  dbetads[3] = {
      dbtmpone*u_x[uOff_x[FIELD_SATURATION]+0] + dbtmptwo*u_x[uOff_x[FIELD_PRESSURE]+0] ,
      dbtmpone*u_x[uOff_x[FIELD_SATURATION]+1] + dbtmptwo*u_x[uOff_x[FIELD_PRESSURE]+1] ,
      dbtmpone*u_x[uOff_x[FIELD_SATURATION]+2] + dbtmptwo*u_x[uOff_x[FIELD_PRESSURE]+2]  };
  double  innerprod = 0.0;
  for (d = 0; d < dim; ++d)  innerprod = innerprod + u_x[uOff_x[FIELD_TEMPERATURE]+d]  * beta[d];

  for (iii = 0; iii < dim; ++iii)  {
    g2[iii] = innerprod * constants[PARAM_ARTIFICIALDIFFUSION] * dbetads[iii];
  }
  for (iii = 0; iii < dim; ++iii) for (jjj = 0; jjj < dim; ++jjj) {
    g2[iii] = g2[iii] + u_x[uOff_x[FIELD_TEMPERATURE]+iii] * beta[jjj] * constants[PARAM_ARTIFICIALDIFFUSION] * dbetads[jjj];
  }
}

static void g3_tempsat(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscInt d,iii,jjj;
  PetscReal  dpds    = -constants[PARAM_DISPLACEMENTPRESSURE]/2.0/(sqrt(1- PetscMin(1, u[FIELD_SATURATION])+_globalepsilon)+_globalepsilon) ;
  // buffers
  PetscReal  tmpone  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]*dpds;
  PetscReal  tmptwo  =-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]* ((1-u[FIELD_SATURATION])*constants[PARAM_MOBILITYBLOOD] + u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]); 
  PetscReal  tmpthree=-constants[PARAM_ADVECTIONTERM]*constants[PARAM_POROSITY]*u[FIELD_SATURATION]*constants[PARAM_MOBILITYOIL]*dpds;
  PetscReal  beta[3] = {
    tmpone*u_x[uOff_x[FIELD_SATURATION]+0] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+0] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+1] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+1] ,
    tmpone*u_x[uOff_x[FIELD_SATURATION]+2] + tmptwo*u_x[uOff_x[FIELD_PRESSURE]+2]  };

  //PetscPrintf(PETSC_COMM_WORLD, "%f ",conduction );
  double  innerprod = 0.0;
  for (d = 0; d < dim; ++d)  innerprod = innerprod + u_x[uOff_x[FIELD_TEMPERATURE]+d]  * beta[d];

  for (iii = 0; iii < dim; ++iii) for (jjj = 0; jjj < dim; ++jjj) {
    g3[iii*dim+jjj] = u_x[uOff_x[FIELD_TEMPERATURE]+iii] * beta[jjj];
  }
  for (d = 0; d < dim; ++d) {
    g3[d*dim+d] = g3[d*dim+d] + innerprod;
  }
  for (iii = 0; iii < dim; ++iii) for (jjj = 0; jjj < dim; ++jjj) {
    g3[iii*dim+jjj] = g3[iii*dim+jjj] * tmpthree  * constants[PARAM_ARTIFICIALDIFFUSION] ;
  }
}

static void f0_conc(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  f0[0] = - constants[PARAM_POROSITY]*u_t[FIELD_SATURATION]
          - constants[PARAM_SATURATION_SOURCE]*u[FIELD_SATURATION];
}
static void f0_bd_conc(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], const PetscReal n[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  f0[0] = - 0.5 * constants[PARAM_INJECTIONVELOCITY];
}
static void g0_conc(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  g0[0] = - u_tShift*1.0 * constants[PARAM_POROSITY]   - constants[PARAM_SATURATION_SOURCE];
}

static void f1_conc(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) f1[d] = constants[PARAM_MOBILITYBLOOD] * u_x[uOff_x[FIELD_PRESSURE]+d]  ;
}

static void g3_conc(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) g3[d*dim+d] = constants[PARAM_MOBILITYBLOOD] ;
  
}

static void f0_damg(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  f0[0] =  0.0; // steady state for no solution change during solve... post process solution
}

static void g0_damg(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  g0[0] = u_tShift*1.0 ;
}

// phase field
static void f0_phas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  // F(c) = 2c^2 (c-1)^2 - 1/8
  PetscScalar  doublewell =  4.*u[FIELD_PHASE] *(u[FIELD_PHASE] -1.)*(u[FIELD_PHASE] -1.) + 4. * u[FIELD_PHASE] * u[FIELD_PHASE] * (u[FIELD_PHASE] - 1.)  ; 
  f0[0] = u_t[FIELD_PHASE] + doublewell ;
}

static void f1_phas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscInt d;
  const PetscScalar  epsilon = constants[PARAM_PHASEEPSILON]; 
  for (d = 0; d < dim; ++d) {
    f1[d] = epsilon*epsilon * u_x[d];
  }
}

static void g3_phas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscInt d;
  const PetscScalar  epsilon = constants[PARAM_PHASEEPSILON]; 
  for (d = 0; d < dim; ++d) {
    g3[d*dim+d] = epsilon*epsilon ;
  }
}

static void g0_phas(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  // F(c)   = 2c^2 (c-1)^2 - 1/8
  // dFdc   = 4c  (c-1)^2 +  4 c^2 (c-1)
  // d2Fdc2 = 4  (c-1)^2 + 8c (c-1) +  8 c (c-1) +   4 c^2
  PetscScalar  ddoublewelldc =  4. * (u[FIELD_PHASE]-1)* (u[FIELD_PHASE]-1) + 8.*u[FIELD_PHASE]*(u[FIELD_PHASE]-1.) +  8.* u[FIELD_PHASE]*(u[FIELD_PHASE]-1) +   4*u[FIELD_PHASE]* u[FIELD_PHASE] ;
  g0[0] = u_tShift*1.0 + ddoublewelldc ;
}
/* ------------------------------------------------------------------- */
/*
   compute objective on subspace
*/
PetscErrorCode subspaceobjective(SNES snes,Vec X,PetscReal *f,void *voidctx)
{
  AppCtx *ctx = (AppCtx*)voidctx;
  Vec            complementresidual,work;
  PetscErrorCode    ierr;

  // compute residual
  ierr = VecDuplicate(X,&work);CHKERRQ(ierr);
  ierr = SNESComputeFunction(snes,X,work);CHKERRQ(ierr);

  // evaluate objective function only on pressure and saturation equations
  ierr = VecGetSubVector(work, ctx->isnotpressuresaturation, &complementresidual);CHKERRQ(ierr);
  ierr = VecSet(complementresidual,0.0);CHKERRQ(ierr);
  ierr = VecRestoreSubVector(work, ctx->isnotpressuresaturation, &complementresidual);CHKERRQ(ierr);
  ierr = VecNormBegin(work, NORM_2, f);CHKERRQ(ierr);
  ierr = VecNormEnd(  work, NORM_2, f);CHKERRQ(ierr);

  // clean up
  ierr = VecDestroy( &work);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
   PreCheck - Optional user-defined routine that checks the validity of
   candidate steps of a line search method.  Set by SNESLineSearchSetPreCheck().

   Input Parameters:
   snes - the SNES context
   xcurrent - current solution
   y - search direction and length

   Output Parameters:
   y         - proposed step (search direction and length) (possibly changed)
   changed_y - tells if the step has changed or not
 */
PetscErrorCode myprecheck(SNESLineSearch linesearch,Vec xcurrent,Vec y, PetscBool *changed_y, void * voidctx)
{
  AppCtx *ctx = (AppCtx*)voidctx;
  Vec            temperature,temperaturesearch;
  PetscErrorCode    ierr;

  PetscFunctionBeginUser;
  // zero search direction for temperature and update from linear linear solve with lambda = 1
  //  ie solve linear system with SNESSolve_KSPONLY
  ierr = VecGetSubVector(y, ctx->fields[FIELD_TEMPERATURE], &temperaturesearch);CHKERRQ(ierr);
  ierr = VecGetSubVector(xcurrent, ctx->fields[FIELD_TEMPERATURE], &temperature);CHKERRQ(ierr);
  ierr = VecAXPY(temperature,-1.0,temperaturesearch);CHKERRQ(ierr);
  ierr = VecSet(temperaturesearch,0.0);CHKERRQ(ierr);
  ierr = VecRestoreSubVector(y, ctx->fields[FIELD_TEMPERATURE], &temperaturesearch);CHKERRQ(ierr);
  ierr = VecRestoreSubVector(xcurrent, ctx->fields[FIELD_TEMPERATURE], &temperature);CHKERRQ(ierr);
  *changed_y = PETSC_TRUE;
  PetscFunctionReturn(0);
}

static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  options->dim     = 3;
  options->simplex = PETSC_TRUE;
  options->variableCoefficient = COEFF_NONE;
  options->fieldBC             = PETSC_FALSE;

  // set initial parameters
  // FIXME - need units for all parameters
  options->parameters[PARAM_OMEGA               ] = 6.0 ;      // [kg/m^3/s]  blood perfusion
  options->parameters[PARAM_RHOBLOOD            ] = 1.0e3;     // [kg/m^3] water density
  options->parameters[PARAM_RHOOIL              ] = 0.8e3;     // [kg/m^3] oil density
  // https://en.wikipedia.org/wiki/Dichloroacetyl_chloride
  options->parameters[PARAM_RHODCACL            ] = 1.531e3;     // [kg/m^3] dcacl density
  options->parameters[PARAM_SPECIFICHEATBLOOD   ] = 3600.0;    // [J/kg/K]
  options->parameters[PARAM_SPECIFICHEATTISSUE  ] = 3600.0;    // [J/kg/K]


  // normalize temperature to units on [0,1] scale 
  options->temperaturescaling      = 100.;      // [hK] [hecto Kelvin]
  options->parameters[PARAM_USALT               ] = 100./options->temperaturescaling;      // [hC] [hecto Celcius]
  options->parameters[PARAM_UARTERY             ] = 37./options->temperaturescaling;       // [hC] [hecto Celcius]
  options->parameters[PARAM_POROSITY            ] = 0.05;       // [volume fraction]

 
  double     conduction              = 5.27;      // [W/m/K]
  // saturation equation effectively a first order ODE
  // http://hyperphysics.phy-astr.gsu.edu/hbase/Math/deinhom.html
  // analytical solutions guide the appropriate time constants
  double     gammaconst              = 0.052/6.*20.;  // [1/s] FIXME - verify this time constant
  // Unconsolidated sands may have permeabilities of over 5000 md = 5 d =  5e-12 m^2
  // https://en.wikipedia.org/wiki/Permeability_(earth_sciences)
  double     tissue_permeability     = 5.e-10;   // [m^2]

  // https://www.engineersedge.com/physics/water__density_viscosity_specific_weight_13146.htm
  double     water_viscosity         = 8.9e-4;   // [Pa s]
  double     oil_viscosity           = 10.e-3;   // [Pa s]

  //For reference, atmospheric pressure is 101kPa and typical injection pressure
  // is 698 - 2792 kPa~\cite{vuvckovic2006injection}.
  // 1 atm = 101325 Pa
  double     atmosphericpressure             =  101325 ;      // [Pa]
  double     injectionpressure               = 2000.e3 ;      // [Pa]

  // https://en.wikipedia.org/wiki/Blood_pressure
  // Blood pressure is one of the vital signs, along with respiratory rate, heart rate, oxygen saturation, and body temperature. Normal resting blood pressure in an adult is approximately 120 millimetres of mercury (16 kPa) systolic, and 80 millimetres of mercury (11 kPa) diastolic, abbreviated "120/80 mmHg
  // Blood pressure is relative to atmospheric pressure.  
  double     systolicpressure               = atmosphericpressure  + 16.e3 ;      // [Pa]
  double     diastolicpressure              = atmosphericpressure  + 11.e3 ;      // [Pa]

  // An approximation to the vessel flow velocity will be obtain from Poiseuille flow solution
  // - \frac{dp}{dx} \left[\frac{Pa}{m}\right]\frac{1}{ 4 \mu} \left[\frac{1}{Pa \; s}\right] \left(  r_0^2 - r^2\right) \left[m^2\right]
  double     lumenradius              = .001 ;      // [m]
  double     pressuregradientlength   = .001 ;      // [m]
  double     poiseuillevelocity       = (systolicpressure - diastolicpressure)/pressuregradientlength/4./water_viscosity * lumenradius  * lumenradius ;     // [m/s]

  // https://www.ahajournals.org/doi/10.1161/01.CIR.40.5.603
  // verage peak and mean blood velocities were 66 and 11 cm/sec in the ascending aorta, 57 and 10 cm/sec in the pulmonary artery, 28 and 12 cm/sec in the superior vena cava, and 26 and 13 cm/sec in the inferior vena cava.
  options->parameters[PARAM_INJECTIONVELOCITY   ] = .11;     // [m/s]
  // Peaceman - capillary pressure at wetting phase saturation = 0 is ~25 inches of water
  // 30 inches of water =  7465.2 Pa
  // solver in terms of atm for scaling
  options->parameters[PARAM_DISPLACEMENTPRESSURE] = 7.465e3/atmosphericpressure;      // [Pa]
  options->parameters[PARAM_BASELINEPRESSURE    ] = atmosphericpressure/atmosphericpressure;      // [Pa]
  options->parameters[PARAM_BOUNDARYPRESSURE    ] = systolicpressure/atmosphericpressure   ;      // [Pa]
  options->parameters[PARAM_EPSILON             ] = 5.e-1;     // [volume fraction of DCACL]
  // FINITE ELEMENT METHODS FOR LINEAR HYPERBOLIC PROBLEMS  - Claes JOHNSON
  //    artificial diffusion parameter should be on the order of the mesh size, millimeter or submillimeter
  options->parameters[PARAM_ADVECTIONTERM       ] = 1.0;        // [units]
  options->parameters[PARAM_ARTIFICIALDIFFUSION ] = 1.e-7;        // [units]
  // https://en.wikipedia.org/wiki/Dichloroacetyl_chloride
  double     molecularmass           = .147;    // [kg/mole]
  // First In Vivo Test of Thermoembolization: Turning Tissue Against Itself Using Transcatheter Chemistry in a Porcine Model -  Erik N. K. Cressman  • Chunxiao Guo
  double     heatofreaction          = 93.e3;   // [J/mole]


  ierr = PetscOptionsBegin(comm, "", "Thermoembolization model parameter options", "DMPLEX");CHKERRQ(ierr);
  PetscBool      flg;
  options->refine = 0;
  ierr = PetscOptionsInt("-dim", "The topological mesh dimension", "ex45.c", options->dim, &options->dim, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-dm_refine", "The number of uniform refinements", "DMCreate", options->refine, &options->refine, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-simplex", "Simplicial (true) or tensor (false) mesh", "ex45.c", options->simplex, &options->simplex, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-artdiff", "artificial diffusion [...]", "ex45.c", options->parameters[PARAM_ARTIFICIALDIFFUSION], &options->parameters[PARAM_ARTIFICIALDIFFUSION], NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-advection", "scale temperature advection term[...]", "ex45.c", options->parameters[PARAM_ADVECTIONTERM], &options->parameters[PARAM_ADVECTIONTERM], NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-velocity", "applicator injection velocity [...]", "ex45.c", options->parameters[PARAM_INJECTIONVELOCITY], &options->parameters[PARAM_INJECTIONVELOCITY], NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-salttemp", "salt temperature [...]", "ex45.c", options->parameters[PARAM_USALT], &options->parameters[PARAM_USALT], NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-disppressure", "displacement pressure [...]", "ex45.c", options->parameters[PARAM_DISPLACEMENTPRESSURE], &options->parameters[PARAM_DISPLACEMENTPRESSURE], NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-permeability", "tissue permeability [...]", "ex45.c", tissue_permeability, &tissue_permeability, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-gamma", "time constanst permeability [...]", "ex45.c", gammaconst, &gammaconst, NULL);CHKERRQ(ierr);
  strcpy(options->imagefile, "./vessel.vtk");
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
       ierr = PetscPrintf(PETSC_COMM_WORLD, "ZBounds {%10.3e,%10.3e} \n",
                            options->bounds[4],options->bounds[5]);
       options->ImageData->GetSpacing(options->spacing);
       ierr = PetscPrintf(PETSC_COMM_WORLD, "spacing {%10.3e,%10.3e,%10.3e} \n",
                            options->spacing[0],options->spacing[1],options->spacing[2]);
     }
  else 
     {
       options->ImageData = 0;
     }

  char              *tmpstring;
  ierr = PetscStrcpy(options->filenosuffix,options->imagefile);CHKERRQ(ierr);
  ierr = PetscStrrstr(options->filenosuffix,".vtk",&tmpstring);CHKERRQ(ierr);
  if (tmpstring) tmpstring[0] = 0;
  ierr = PetscOptionsString("-o", "file output", "ex45.c", options->filenosuffix, options->filenosuffix, sizeof(options->filenosuffix), &flg);CHKERRQ(ierr);

  // FIXME error handle time steps. max time  should be  < epsilon^{-1}
  // FIXME epsilon^{-1} ~ 1/2 * voxel width
  options->parameters[PARAM_PHASEEPSILON] = options->spacing[0]; 
  PetscReal         max_time;               /* phase field max time allowed */
  max_time = 1./options->parameters[PARAM_PHASEEPSILON];
  ierr = PetscOptionsInt("-phasepresolve_ts_max_steps","Maximum number of time steps","TSSetMaxSteps",options->max_steps,&options->max_steps,NULL);CHKERRQ(ierr);
  options->time_step = max_time/options->max_steps;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "epsilon should be half voxel size %10.3e, max time step should be less than %10.3e, time step %10.3e \n",options->parameters[PARAM_PHASEEPSILON],max_time,options->time_step);

  

  ierr = PetscOptionsEnd();
  
  options->parameters[PARAM_MOBILITYOIL        ] = tissue_permeability/oil_viscosity  *atmosphericpressure; // [m^2/atm/s]
  options->parameters[PARAM_MOBILITYBLOOD      ] = tissue_permeability/water_viscosity*atmosphericpressure; // [m^2/atm/s]
  options->parameters[PARAM_ALPHA              ] = conduction/ options->parameters[PARAM_RHOBLOOD] / options->parameters[PARAM_SPECIFICHEATBLOOD] ;    // [W/m/K / (kg/m^3) / (J/kg/K)] =      m^s /s
  // convenience
  options->parameters[PARAM_SATURATION_SOURCE  ] = gammaconst*options->parameters[PARAM_RHODCACL]*options->parameters[PARAM_EPSILON]*options->parameters[PARAM_POROSITY]/options->parameters[PARAM_RHOBLOOD];
  options->parameters[PARAM_PRESSURE_SOURCE    ] = gammaconst*options->parameters[PARAM_RHODCACL]*options->parameters[PARAM_EPSILON]*options->parameters[PARAM_POROSITY]*(1./options->parameters[PARAM_RHOBLOOD] - 1./options->parameters[PARAM_RHOOIL]);
  options->parameters[PARAM_TEMPERATURE_SOURCE ] = gammaconst*options->parameters[PARAM_RHODCACL]*options->parameters[PARAM_EPSILON]*options->parameters[PARAM_POROSITY]/options->parameters[PARAM_RHOBLOOD] /options->parameters[PARAM_SPECIFICHEATBLOOD] /molecularmass * heatofreaction/ options->temperaturescaling ; // [(1/s) / (kg/mole) / (J/kg/K) * (J/mole)  (1hK/100K) ] = [hK/s]

  // echo parameters
  ierr = PetscPrintf(PETSC_COMM_WORLD, "MESH FILE                       = %s\n"    ,options->imagefile                              );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "OUTPUT FILE                     = %s\n"    ,options->filenosuffix                          );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "GAMMA                   [1/s]   = %12.5e\n",gammaconst                                     );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "temperaturescaling              = %12.5e\n",options->temperaturescaling                    );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_OMEGA                     = %12.5e\n",options->parameters[PARAM_OMEGA               ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_RHOBLOOD                  = %12.5e\n",options->parameters[PARAM_RHOBLOOD            ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_RHOOIL                    = %12.5e\n",options->parameters[PARAM_RHOOIL              ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_RHODCACL                  = %12.5e\n",options->parameters[PARAM_RHODCACL            ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_ALPHA             [m^2/s] = %12.5e\n",options->parameters[PARAM_ALPHA               ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_SPECIFICHEATBLOOD         = %12.5e\n",options->parameters[PARAM_SPECIFICHEATBLOOD   ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_SPECIFICHEATTISSUE        = %12.5e\n",options->parameters[PARAM_SPECIFICHEATTISSUE  ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_USALT               [hC]  = %12.5e\n",options->parameters[PARAM_USALT               ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_UARTERY             [hC]  = %12.5e\n",options->parameters[PARAM_UARTERY             ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_POROSITY                  = %12.5e\n",options->parameters[PARAM_POROSITY            ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_MOBILITYOIL   [m^2/atm/s] = %12.5e\n",options->parameters[PARAM_MOBILITYOIL         ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_MOBILITYBLOOD [m^2/atm/s] = %12.5e\n",options->parameters[PARAM_MOBILITYBLOOD       ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_INJECTIONVELOCITY         = %12.5e\n",options->parameters[PARAM_INJECTIONVELOCITY   ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_DISPLACEMENTPRESSURE[atm] = %12.5e\n",options->parameters[PARAM_DISPLACEMENTPRESSURE]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_BASELINEPRESSURE    [atm] = %12.5e\n",options->parameters[PARAM_BASELINEPRESSURE    ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_BOUNDARYPRESSURE    [atm] = %12.5e\n",options->parameters[PARAM_BOUNDARYPRESSURE    ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_EPSILON                   = %12.5e\n",options->parameters[PARAM_EPSILON             ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_SATURATION_SOURCE         = %12.5e\n",options->parameters[PARAM_SATURATION_SOURCE   ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_PRESSURE_SOURCE           = %12.5e\n",options->parameters[PARAM_PRESSURE_SOURCE     ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_TEMPERATURE_SOURCE [hK/s] = %12.5e\n",options->parameters[PARAM_TEMPERATURE_SOURCE  ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_ADVECTIONTERM             = %12.5e\n",options->parameters[PARAM_ADVECTIONTERM       ]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "PARAM_ARTIFICIALDIFFUSION       = %12.5e\n",options->parameters[PARAM_ARTIFICIALDIFFUSION ]);CHKERRQ(ierr);
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
  ierr = PetscPrintf(PETSC_COMM_WORLD, "lower = (%14.7e,%14.7e,%14.7e), upper = (%14.7e,%14.7e,%14.7e) \n",lower[0],lower[1],lower[2], upper[0],upper[1],upper[2]);CHKERRQ(ierr);
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
  ierr = DMViewFromOptions(*dm, NULL, "-dm_view");CHKERRQ(ierr);
  ierr = DMSetFromOptions(*dm);CHKERRQ(ierr);
  ierr = DMViewFromOptions(*dm, NULL, "-dm_view");CHKERRQ(ierr);

  // ierr = DMPlexSetRefinementUniform(*dm, PETSC_FALSE);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

static PetscErrorCode SetupProblem(PetscDS prob, AppCtx *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  // temperature equations
  // debug
  // ierr = PetscDSSetResidual(  prob, FIELD_TEMPERATURE, f0_damg, NULL);CHKERRQ(ierr);
  // ierr = PetscDSSetJacobian(  prob, FIELD_TEMPERATURE, FIELD_TEMPERATURE, g0_damg, NULL, NULL, NULL);CHKERRQ(ierr);
  // nonlinear equations for temperature, pressure, and saturation including change in advection velocity wrt (p,s)
  ierr = PetscDSSetResidual(prob, FIELD_TEMPERATURE, f0_temp, f1_temp);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(prob, FIELD_TEMPERATURE, FIELD_TEMPERATURE, g0_temp, g1_temp, NULL, g3_temp);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(prob, FIELD_TEMPERATURE, FIELD_PRESSURE   , NULL   , g1_temppres, NULL, g3_temppres);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(prob, FIELD_TEMPERATURE, FIELD_SATURATION , g0_tempsat, g1_tempsat, g2_tempsat, g3_tempsat);CHKERRQ(ierr);
  
  // //ierr = PetscDSSetBdJacobian(prob, 0, 0, g0_bd_uu_3d, NULL, NULL, NULL);CHKERRQ(ierr);

  // wetting phase pressure equations
  ierr = PetscDSSetResidual(  prob, FIELD_PRESSURE, f0_p, f1_p);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(  prob, FIELD_PRESSURE, FIELD_PRESSURE  , NULL, NULL, NULL,g3_pp);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(  prob, FIELD_PRESSURE, FIELD_SATURATION,g0_ps, NULL,g2_ps,g3_ps);CHKERRQ(ierr);
  ierr = PetscDSSetBdResidual(prob, FIELD_PRESSURE, f0_bd_p, NULL);CHKERRQ(ierr);
  // debug
  // ierr = PetscDSSetResidual(  prob, FIELD_PRESSURE, f0_damg, NULL);CHKERRQ(ierr);
  // ierr = PetscDSSetJacobian(  prob, FIELD_PRESSURE, FIELD_PRESSURE, g0_damg, NULL, NULL, NULL);CHKERRQ(ierr);

  // damage equations
  ierr = PetscDSSetResidual(  prob, FIELD_DAMAGE, f0_damg, NULL);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(  prob, FIELD_DAMAGE, FIELD_DAMAGE, g0_damg, NULL, NULL, NULL);CHKERRQ(ierr);

  // nonwetting phase saturation equations
  ierr = PetscDSSetResidual(  prob, FIELD_SATURATION, f0_conc, f1_conc);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(  prob, FIELD_SATURATION, FIELD_PRESSURE  ,    NULL, NULL, NULL, g3_conc);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(  prob, FIELD_SATURATION, FIELD_SATURATION, g0_conc, NULL, NULL,   NULL );CHKERRQ(ierr);
  // debug
  // ierr = PetscDSSetResidual(  prob, FIELD_SATURATION, f0_conc, NULL);CHKERRQ(ierr);
  // ierr = PetscDSSetJacobian(  prob, FIELD_SATURATION, FIELD_SATURATION, g0_conc, NULL, NULL,   NULL );CHKERRQ(ierr);
  //ierr = PetscDSSetBdResidual(prob, FIELD_SATURATION, f0_bd_conc, NULL);CHKERRQ(ierr);
  
  // phase field
  ierr = PetscDSSetResidual(prob, FIELD_PHASE, f0_phas, f1_phas);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(prob, FIELD_PHASE, FIELD_PHASE, g0_phas, NULL, NULL, g3_phas);CHKERRQ(ierr);

  // evaluate exact solution
  const PetscInt numberfields=5;
  ierr = PetscMalloc1(numberfields, &ctx->exactFuncs);CHKERRQ(ierr);
  ctx->exactFuncs[FIELD_TEMPERATURE] = analytic_temp;
  ctx->exactFuncs[FIELD_PRESSURE   ] = analytic_pres;
  ctx->exactFuncs[FIELD_DAMAGE     ] = analytic_damg;
  ctx->exactFuncs[FIELD_SATURATION ] = analytic_conc;
  ctx->exactFuncs[FIELD_PHASE ]      = analytic_phas;
  const PetscInt id = 1;
  const PetscInt nodeSetApplicatorValue = 2; // node set value assigned to exodus mesh
  const PetscInt nodeSetGroundValue = 3; // node set value assigned to exodus mesh
  const PetscInt nodeSetNeumannBoundaryValue = 4; // node set value assigned to exodus mesh
  //ierr = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, "wall", "marker", 0, 0, NULL, (void (*)(void)) ctx->exactFuncs[0], 1, &id, ctx);CHKERRQ(ierr);
  // Dirichlet BC
  //ierr = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, "applicator", "Face Sets"  , FIELD_TEMPERATURE,  0, NULL, (void(*)())salttemperature   , 1, &nodeSetApplicatorValue      , ctx);CHKERRQ(ierr);
  //ierr = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, "applicator", "Vertex Sets", FIELD_TEMPERATURE,  0, NULL, (void(*)())bodytemperature   , 1, &nodeSetNeumannBoundaryValue , ctx);CHKERRQ(ierr);
  //ierr = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, "applicator", "Face Sets"  , FIELD_PRESSURE   ,  0, NULL, (void(*)())vessel_pres       , 1, &nodeSetApplicatorValue      , ctx);CHKERRQ(ierr);
  ierr = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, "applicator", "marker", FIELD_PRESSURE   ,  0, NULL, (void(*)())baseline_pres     , 1, &id , ctx);CHKERRQ(ierr);
  //ierr = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, "applicator", "Face Sets"  , FIELD_DAMAGE     ,  0, NULL, (void(*)())tissuedamagefcn   , 1, &nodeSetApplicatorValue      , ctx);CHKERRQ(ierr);
  //ierr = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, "applicator", "Face Sets"  , FIELD_SATURATION ,  0, NULL, (void(*)())bolusinjection    , 1, &nodeSetApplicatorValue      , ctx);CHKERRQ(ierr);
  ierr = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, "applicator", "marker", FIELD_SATURATION ,  0, NULL, (void(*)())fieldzero         , 1, &id , ctx);CHKERRQ(ierr);
  // Neuman BC
  //ierr = PetscDSAddBoundary(prob, DM_BC_NATURAL  , "applicator", "Face Sets"  , FIELD_PRESSURE   ,  0, NULL, NULL                         , 1, &nodeSetApplicatorValue      , ctx);CHKERRQ(ierr);
  //ierr = PetscDSAddBoundary(prob, DM_BC_NATURAL  , "applicator", "Face Sets"  , FIELD_SATURATION ,  0, NULL, NULL                         , 1, &nodeSetApplicatorValue      , ctx);CHKERRQ(ierr);
  // DMAddBoundary(*dm, DM_BC_ESSENTIAL, "wall", "boundary", 0, 0, NULL, (void (*)(void)) user->bcFuncs[0], 1, &id, user);
  ierr = PetscDSSetConstants(prob, NUMPARAMETERS, ctx->parameters); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SetupMaterial(DM dm, DM dmAux, AppCtx *user)
{
  PetscErrorCode (*matFuncs[1])(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar u[], void *ctx) = {nu_2d};
  Vec            nu;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMCreateLocalVector(dmAux, &nu);CHKERRQ(ierr);
  ierr = DMProjectFunctionLocal(dmAux, 0.0, matFuncs, NULL, INSERT_ALL_VALUES, nu);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject) dm, "A", (PetscObject) nu);CHKERRQ(ierr);
  ierr = VecDestroy(&nu);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


static PetscErrorCode SetupDiscretization(DM dm, AppCtx* ctx)
{
  DM             cdm = dm;
  const PetscInt dim = ctx->dim;
  PetscFE        feAux = NULL;
  PetscDS        prob, probAux = NULL;
  PetscFE        fe[5];
  PetscQuadrature q;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* Create finite element */
  ierr = PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, ctx->simplex, "phas_", PETSC_DEFAULT, &fe[FIELD_PHASE]);CHKERRQ(ierr);
  ierr = PetscFEGetQuadrature(fe[FIELD_PHASE], &q);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fe[FIELD_PHASE], "phasefield");CHKERRQ(ierr);

  ierr = PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, ctx->simplex, "temp_", -1, &fe[FIELD_TEMPERATURE]);CHKERRQ(ierr);
  ierr = PetscFESetQuadrature(fe[FIELD_TEMPERATURE], q);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fe[FIELD_TEMPERATURE], "temperature");CHKERRQ(ierr);

  ierr = PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, ctx->simplex, "damg_", PETSC_DEFAULT, &fe[FIELD_DAMAGE]);CHKERRQ(ierr);
  ierr = PetscFESetQuadrature(fe[FIELD_DAMAGE], q);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fe[FIELD_DAMAGE], "damage");CHKERRQ(ierr);

  ierr = PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, ctx->simplex, "pres_", PETSC_DEFAULT, &fe[FIELD_PRESSURE]);CHKERRQ(ierr);
  ierr = PetscFESetQuadrature(fe[FIELD_PRESSURE], q);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fe[FIELD_PRESSURE], "pressure");CHKERRQ(ierr);

  ierr = PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, ctx->simplex, "conc_", PETSC_DEFAULT, &fe[FIELD_SATURATION]);CHKERRQ(ierr);
  ierr = PetscFESetQuadrature(fe[FIELD_SATURATION], q);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fe[FIELD_SATURATION], "concentration");CHKERRQ(ierr);

  /* Setup Variable Coeff */
  if (ctx->variableCoefficient == COEFF_FIELD) {

    ierr = PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, ctx->simplex, "mat_", -1, &feAux);CHKERRQ(ierr);
    ierr = PetscFESetQuadrature(feAux, q);CHKERRQ(ierr);
    ierr = PetscDSCreate(PetscObjectComm((PetscObject)dm),&probAux);CHKERRQ(ierr);
    ierr = PetscDSSetDiscretization(probAux, 0, (PetscObject) feAux);CHKERRQ(ierr);
  } 

  /* Set discretization and boundary conditions for each mesh */
  ierr = DMGetDS(dm, &prob);CHKERRQ(ierr);
  ierr = PetscDSSetDiscretization(prob, 0, (PetscObject) fe[0]);CHKERRQ(ierr);
  ierr = PetscDSSetDiscretization(prob, 1, (PetscObject) fe[1]);CHKERRQ(ierr);
  ierr = PetscDSSetDiscretization(prob, 2, (PetscObject) fe[2]);CHKERRQ(ierr);
  ierr = PetscDSSetDiscretization(prob, 3, (PetscObject) fe[3]);CHKERRQ(ierr);
  ierr = PetscDSSetDiscretization(prob, 4, (PetscObject) fe[4]);CHKERRQ(ierr);
  ierr = SetupProblem(prob, ctx);CHKERRQ(ierr);
  while (cdm) {
    PetscBool hasLabel;

    ierr = DMSetDS(cdm, prob);CHKERRQ(ierr);
    ierr = DMHasLabel(cdm, "marker", &hasLabel);CHKERRQ(ierr);
    if (!hasLabel) {ierr = CreateBCLabel(cdm, "marker");CHKERRQ(ierr);}
    ierr = DMGetCoarseDM(cdm, &cdm);CHKERRQ(ierr);
  }
  ierr = PetscFEDestroy(&fe[0]);CHKERRQ(ierr);
  ierr = PetscFEDestroy(&fe[1]);CHKERRQ(ierr);
  ierr = PetscFEDestroy(&fe[2]);CHKERRQ(ierr);
  ierr = PetscFEDestroy(&fe[3]);CHKERRQ(ierr);
  ierr = PetscFEDestroy(&fe[4]);CHKERRQ(ierr);
  ierr = PetscFEDestroy(&feAux);CHKERRQ(ierr);
  ierr = PetscDSDestroy(&probAux);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc, char **argv)
{
  AppCtx         ctx;
  DM             dm;
  TS             ts,prets;
  Vec            u, r;
  PetscReal      t       = 0.0;
  PetscReal      L2error = 0.0;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, help);if (ierr) return ierr;
  ierr = ProcessOptions(PETSC_COMM_WORLD, &ctx);CHKERRQ(ierr);
  ierr = CreateMesh(PETSC_COMM_WORLD, &dm, &ctx);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(dm, &ctx);CHKERRQ(ierr);
  ierr = SetupDiscretization(dm, &ctx);CHKERRQ(ierr);

  // get index subsets
  ierr = DMCreateFieldIS(dm, &ctx.numFields, &ctx.fieldNames, &ctx.fields);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm, &u);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) u, "solution");CHKERRQ(ierr);
  ierr = VecDuplicate(u, &r);CHKERRQ(ierr);

  // time stepper for phase field
  ierr = TSCreate(PETSC_COMM_WORLD, &prets);CHKERRQ(ierr);
  ierr = TSSetDM(prets, dm);CHKERRQ(ierr);
  ierr = DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, &ctx);CHKERRQ(ierr);
  ierr = DMTSSetIFunctionLocal(dm, DMPlexTSComputeIFunctionFEM, &ctx);CHKERRQ(ierr);
  ierr = DMTSSetIJacobianLocal(dm, DMPlexTSComputeIJacobianFEM, &ctx);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(prets, TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  ierr = TSSetOptionsPrefix(prets,"phasepresolve_");CHKERRQ(ierr);
  ierr = TSSetFromOptions(prets);CHKERRQ(ierr);
  ierr = TSSetTimeStep(prets,ctx.time_step);CHKERRQ(ierr);
  char              prevtkfilenametemplate[PETSC_MAX_PATH_LEN];
  ierr = PetscSNPrintf(prevtkfilenametemplate,sizeof(prevtkfilenametemplate),"%spre.%%04d.vtu",ctx.filenosuffix);CHKERRQ(ierr);
  ierr = TSMonitorSet(prets,TSMonitorSolutionVTK,&prevtkfilenametemplate,NULL);CHKERRQ(ierr);
   { // PCApply_FieldSplit
    SNES presnes;
    KSP  preksp;
    PC   prepc;
    ierr = TSGetSNES(prets,&presnes);
    ierr = SNESGetKSP(presnes,&preksp);CHKERRQ(ierr);
    ierr = KSPGetPC(preksp,&prepc);CHKERRQ(ierr);
    ierr = PCFieldSplitSetIS(prepc,"c",ctx.fields[FIELD_PHASE]);CHKERRQ(ierr);
   }

  // time stepper for coupled equations
  ierr = TSCreate(PETSC_COMM_WORLD, &ts);CHKERRQ(ierr);
  ierr = TSSetDM(ts, dm);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  // same application context for each field
  void          *ctxarray[ctx.numFields];
  ctxarray[0] = &ctx; ctxarray[1] = &ctx; ctxarray[2] = &ctx; ctxarray[3] = &ctx; ctxarray[4] = &ctx;
  ierr = DMProjectFunction(dm, t, ctx.exactFuncs, ctxarray, INSERT_ALL_VALUES, u);CHKERRQ(ierr);

  //ierr = TSMonitorSet(ts,TSMonitorSolutionVTK,&ctx,(void*)&TSMonitorSolutionVTKDestroy);CHKERRQ(ierr);
  // write vtk file at every time point
  char              vtkfilenametemplate[PETSC_MAX_PATH_LEN];
  ierr = PetscSNPrintf(vtkfilenametemplate,sizeof(vtkfilenametemplate),"%ssolution.%%04d.vtu",ctx.filenosuffix);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,TSMonitorSolutionVTK,&vtkfilenametemplate,NULL);CHKERRQ(ierr);
  ierr = TSSetPostStage(ts,TSUpdateArrhenius);CHKERRQ(ierr);

  IS   ispressuresaturation;
   { // PCApply_FieldSplit
    SNES mysnes;
    KSP  myksp;
    PC   mypc;
    SNESLineSearch mylinesearch;
    ierr = TSGetSNES(ts,&mysnes);
    ierr = SNESGetKSP(mysnes,&myksp);CHKERRQ(ierr);
    ierr = SNESGetLineSearch(mysnes,&mylinesearch);CHKERRQ(ierr);
    ierr = SNESSetObjective(mysnes,subspaceobjective,&ctx); CHKERRQ(ierr);
    ierr = SNESLineSearchSetPreCheck(mylinesearch,myprecheck,&ctx); CHKERRQ(ierr);

    ierr = KSPGetPC(myksp,&mypc);CHKERRQ(ierr);
    ierr = ISConcatenate(PETSC_COMM_WORLD,2,&ctx.fields[FIELD_PRESSURE],&ispressuresaturation); CHKERRQ(ierr);
    ierr = ISSort(ispressuresaturation);CHKERRQ(ierr);
    ierr = ISConcatenate(PETSC_COMM_WORLD,3,&ctx.fields[FIELD_PHASE],&ctx.isnotpressuresaturation); CHKERRQ(ierr);
    ierr = PCFieldSplitSetIS(mypc,"s",ispressuresaturation);CHKERRQ(ierr);
    ierr = PCFieldSplitSetIS(mypc,"u",ctx.fields[FIELD_TEMPERATURE]);CHKERRQ(ierr);
   }

  // first solve is for phase field
  char   phasefieldsolution[PETSC_MAX_PATH_LEN];
  ierr = PetscSNPrintf(phasefieldsolution,sizeof(phasefieldsolution),"./vector.%04d.dat",ctx.refine);CHKERRQ(ierr);

  PetscBool      flg;
  ierr = PetscTestFile(phasefieldsolution, 'r', &flg);CHKERRQ(ierr);

  if( flg == PETSC_TRUE ) 
    {/* Read in previously computed solution in binary format */
     PetscViewer    viewer;
     ierr = PetscPrintf(PETSC_COMM_WORLD,"reading vector in binary from vector.dat ...\n"); CHKERRQ(ierr);
     ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,phasefieldsolution,FILE_MODE_READ,&viewer); CHKERRQ(ierr);
     ierr = VecLoad(u,viewer); CHKERRQ(ierr);
     ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    } 
  else
    {/* solve and write phase field solution to disk in vector in binary format */
     PetscViewer    viewer;
     ierr = TSSolve(prets, u);CHKERRQ(ierr);
     ierr = PetscPrintf(PETSC_COMM_WORLD,"writing vector in binary to vector.dat ...\n"); CHKERRQ(ierr);
     ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,phasefieldsolution,FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
     ierr = VecView(u,viewer); CHKERRQ(ierr);
     ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    } 


  // solve full problem 
  ierr = TSSolve(ts, u);CHKERRQ(ierr);

  ierr = TSGetTime(ts, &t);CHKERRQ(ierr);
  // ierr = DMComputeL2Diff(dm, t, ctx.exactFuncs, NULL, u, &L2error);CHKERRQ(ierr);
  // if (L2error < 1.0e-11) {ierr = PetscPrintf(PETSC_COMM_WORLD, "L_2 Error: < 1.0e-11\n");CHKERRQ(ierr);}
  // else                   {ierr = PetscPrintf(PETSC_COMM_WORLD, "L_2 Error: %g\n", (double)L2error);CHKERRQ(ierr);}
  ierr = VecViewFromOptions(u, NULL, "-sol_vec_view");CHKERRQ(ierr);

  // clean up
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = ISDestroy(&ispressuresaturation);CHKERRQ(ierr);
  ierr = ISDestroy(&ctx.isnotpressuresaturation);CHKERRQ(ierr);

  for(PetscInt iii = 0 ; iii < ctx.numFields; iii++)
    {
     ierr = PetscFree(ctx.fieldNames[iii]);CHKERRQ(ierr);
     ierr = ISDestroy(&ctx.fields[iii]);CHKERRQ(ierr);
    }
  ierr = PetscFree(ctx.fieldNames);CHKERRQ(ierr);
  ierr = PetscFree(ctx.fields);CHKERRQ(ierr);
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
    args: -dm_refine 1 -temp_petscspace_order 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_p1_r3
    requires: triangle
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dm_refine 3 -temp_petscspace_order 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_p1_r5
    requires: triangle
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dm_refine 5 -temp_petscspace_order 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_p2_r1
    requires: triangle
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dm_refine 1 -temp_petscspace_order 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_p2_r3
    requires: triangle
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dm_refine 3 -temp_petscspace_order 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_p2_r5
    requires: triangle
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dm_refine 5 -temp_petscspace_order 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q1_r1
    requires: !single
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 1 -temp_petscspace_order 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q1_r3
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 3 -temp_petscspace_order 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q1_r5
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 5 -temp_petscspace_order 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q2_r1
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 1 -temp_petscspace_order 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q2_r3
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 3 -temp_petscspace_order 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q2_r5
    requires: !single
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 5 -temp_petscspace_order 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p1_r1
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 1 -temp_petscspace_order 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p1_r2
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 2 -temp_petscspace_order 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p1_r3
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 3 -temp_petscspace_order 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p2_r1
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 1 -temp_petscspace_order 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p2_r2
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 2 -temp_petscspace_order 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p2_r3
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 3 -temp_petscspace_order 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q1_r1
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 1 -temp_petscspace_order 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q1_r2
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 2 -temp_petscspace_order 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q1_r3
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 3 -temp_petscspace_order 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q2_r1
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 1 -temp_petscspace_order 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q2_r2
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 2 -temp_petscspace_order 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q2_r3
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 3 -temp_petscspace_order 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor

TEST*/
