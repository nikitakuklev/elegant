/* Copyright 2015 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* program: bratSubroutines.c
 * purpose: Bending magnet RAy Trace program for gridded data for dipole magnet.
 *           
 * The equations of motion are
 *        -
 *      d w    -1  -   - -
 *      --- = ---- w x B(q)
 *      d s    H
 *
 *        -
 *      d q    -
 *      --- =  w
 *      d s
 *
 * where B is the magnetic field, s is the distance,
 * q is the coordinate, w is the rate of change of q with
 * distance, and H is the rigidity
 * q(0) = z in the field input file
 * q(1) = x in the field input file
 * q(2) = y
 */
#include "mdb.h"
#include "SDDS.h"
#include "track.h"

void BRAT_B_field(double *B, double *Q);
double BRAT_setup_field_data(char *input, double xCenter, double zCenter);
double BRAT_setup_arc_field_data(char *input, char *sName, char *fieldName, double xCenter);
void BRAT_setup_single_scan_field_data(char *input);
double BRAT_exit_function(double *qp, double *q, double s);
void BRAT_deriv_function(double *qp, double *q, double s);
void BRAT_lorentz_integration(double *accelCoord, double *q, long doStoreData);
void BRAT_store_data(double *qp, double *q, double s, double exval);
void BRAT_optimize_magnet(unsigned long flags);
double refineAngle(double theta, double z0, double x0, double zv, double xv, 
                   double z1, double x1);
static long verbose_optimize = 0;

/* parameters of element needed for integration */
static double theta, fse=0;
static double dXOffset = 0, dZOffset = 0;
static double magnetYaw = 0;
static double fieldFactor = 1;
static long zDuplicate = 0;
static double rhoMax = 0;
static double fieldSign = 1;

/* parameters for field calculation: */
/* Bnorm is misnamed here, based on earlier versions of the program */
static double **Bnorm, **dBnormdz, **dBnormdx;
static double *BxNorm, *ByNorm, *BzNorm;
//static double Breference = 0;
static double xi, xf, dx;
static double yi, yf, dy;
static double zi, zf, dz, z_outer;
static long nx, ny, nz;
static long extend_data, particle_inside, single_scan=0, arc_scan=0;
static double extend_edge_angle;
static double gap = 0.04;
static double central_length;

static double *BArcNorm, *dBArcNormds, *sArc, rhoArc;
static long nArcPoints;

/* tolerance for integration, zero-finding */
static double integ_tol, zero_tol;
static double *X_stored, *Z_stored, *Y_stored;
static double *wX_stored, *wZ_stored, *wY_stored;
static double *aX_stored, *aZ_stored, *aY_stored;
static double *FX_stored, *FZ_stored, *FY_stored, *Fint_stored;
static double *s_stored;
static long n_stored, max_store;

/* particle trajectory parameters */
static double zStart, zEnd;            /* starting and ending point of integration */
static double xNomEntry, zNomEntry;    /* coordinates of hard-edge entry */
static double xNomExit, zNomExit;      /* coordinates of hard-edge exit */
static double xVertex, zVertex;
static double rigidity;
static double xCenter, zCenter;

#define TOLERANCE_FACTOR 1e-14
#define OPTIMIZE_ON        0x0001
#define OPTIMIZE_FSE       0x0002
#define OPTIMIZE_DX        0x0004
#define OPTIMIZE_VERBOSE   0x0008
#define OPTIMIZE_QUIET     0x0010
#define OPTIMIZE_DZ        0x0020
#define OPTIMIZE_YAW       0x0040

static long quiet=1;
static long useFTABLE = 0;

static short idealMode = 0, fieldMapDimension=2;
//static double idealB;
static double idealChord, idealEdgeAngle;

static double fseLimit[2] = {-1, 1};
static double dxLimit[2] = {-1, 1};
static double dzLimit[2] = {-1, 1};
static double yawLimit[2] = {-1, 1};

typedef struct {
  char *filename;
  long nx, ny, nz;
  double *Bx, *By, *Bz;
  double xi, xf, dx;
  double yi, yf, dy;
  double zi, zf, dz;
  double Bmin, Bmax;
} BRAT_3D_DATA;
static BRAT_3D_DATA *brat3dData = NULL;
static long nBrat3dData = 0;

long trackBRAT(double **part, long np, BRAT *brat, double pCentral, double **accepted)
{
  long ip, ic;
  
  if (!brat->initialized) {
    double *xd, *yd, *zd, *Bxd, *Byd, *Bzd;
    double Bmin, Bmax;
    long i, rows, idata;
    SDDS_DATASET SDDS_table;
    
    xd = yd = zd = Bxd = Byd = Bzd = NULL;
    
    /* See if we've read this file already---should use a hash table */
    for (i=0; i<nBrat3dData; i++) {
      if (strcmp(brat->filename, brat3dData[i].filename)==0)
        break;
    }

    if (i!=nBrat3dData) {
      brat->dataIndex = i;
    } else {
      /* read data from file */
      
      if (!SDDS_InitializeInput(&SDDS_table, brat->filename) || !SDDS_ReadPage(&SDDS_table)) {
        SDDS_SetError("Unable to read BRAT data file");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (!(rows=SDDS_CountRowsOfInterest(&SDDS_table)))
        bomb("no data in BRAT field file", NULL);
      if (SDDS_CheckColumn(&SDDS_table, "x", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY || 
          SDDS_CheckColumn(&SDDS_table, "y", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY ||
          SDDS_CheckColumn(&SDDS_table, "z", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY) {
        exit(1);
      }
      if (SDDS_CheckColumn(&SDDS_table, "Bx", "T", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY || 
          SDDS_CheckColumn(&SDDS_table, "By", "T", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY ||
          SDDS_CheckColumn(&SDDS_table, "Bz", "T", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY) {
        exit(1);
      }
      if (!(xd = SDDS_GetColumnInDoubles(&SDDS_table, "x")) ||
          !(yd = SDDS_GetColumnInDoubles(&SDDS_table, "y")) ||
          !(zd = SDDS_GetColumnInDoubles(&SDDS_table, "z")) ||
          !(Bxd = SDDS_GetColumnInDoubles(&SDDS_table, "Bx")) ||
          !(Byd = SDDS_GetColumnInDoubles(&SDDS_table, "By")) ||
          !(Bzd = SDDS_GetColumnInDoubles(&SDDS_table, "Bz"))) {
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      if (!SDDS_Terminate(&SDDS_table)) 
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

      /* It is assumed that the data is ordered so that x changes fastest.
       * This can be accomplished with sddssort -column=z,incr -column=y,incr -column=x,incr
       * The points are assumed to be equipspaced.
       */
      nx = 1;
      xi = xd[0];
      while (nx<rows) {
        if (xd[nx-1]>xd[nx])
          break;
        nx ++;
      }
      if (nx==rows) {
        fprintf(stderr, "BRAT file doesn't have correct structure or amount of data (x)\n");
        fprintf(stderr, "Use sddssort -column=z -column=y -column=x to sort the file\n");
        exit(1);
      }  
      xf = xd[nx-1];
      dx = (xf-xi)/(nx-1);
      
      ny = 1;
      yi = yd[0];
      while (ny<(rows/nx)) {
        if (yd[(ny-1)*nx]>yd[ny*nx])
          break;
        ny++;
      }
      if (ny==rows) {
        fprintf(stderr, "BRAT file doesn't have correct structure or amount of data (y)\n");
        fprintf(stderr, "Use sddssort -column=z -column=y -column=x to sort the file\n");
        exit(1);
      }
      yf = yd[(ny-1)*nx];
      dy = (yf-yi)/(ny-1);
      
      if (nx<=1 || ny<=1 || (nz = rows/(nx*ny))<=1) {
        fprintf(stderr, "BRAT file doesn't have correct structure or amount of data (z)\n");
        fprintf(stderr, "Use sddssort -column=z -column=y -column=x to sort the file\n");
        exit(1);
      }
      zi = zd[0];
      zf = zd[rows-1];
      dz = (zf-zi)/(nz-1);
      
      Bmin = -(Bmax = -DBL_MAX);
      for (idata=0; idata<rows; idata++) {
        if (fabs(yd[idata])<dy/2 && fabs(xd[idata])<dx/2) {
          if (Byd[idata]>Bmax)
            Bmax = Byd[idata];
          if (Byd[idata]<Bmin)
            Bmin = Byd[idata];
        }
      }
      free(xd);
      free(yd);
      free(zd);
      if (Bmax==-DBL_MAX) {
        fprintf(stderr, "BRAT file doesn't have valid Bmax value---verify that values with x=y=0 exist\n");
        exit(1);
      }
      if (fabs(Bmin)>fabs(Bmax))
        SWAP_DOUBLE(Bmin, Bmax);
      
      if (!quiet) 
        printf("3D BRAT field map data: nx=%ld, ny=%ld, nz=%ld\ndx=%e, dy=%e, dz=%e\nx:[%e, %e], y:[%e, %e], z:[%e, %e]\nBy:[%e, %e]\n",
                nx, ny, nz, dx, dy, dz,
                xi, xf, yi, yf, zi, zf,
                Bmin, Bmax
                );

      /* copy to storage area */
      if (!(brat3dData = SDDS_Realloc(brat3dData, sizeof(*brat3dData)*(nBrat3dData+1))))
        bombElegant("memory allocation failure storage BRAT data", NULL);

      brat3dData[nBrat3dData].Bmin = Bmin;
      brat3dData[nBrat3dData].Bmax = Bmax;
      
      brat3dData[nBrat3dData].Bx = Bxd;
      brat3dData[nBrat3dData].xi = xi;
      brat3dData[nBrat3dData].xf = xf;
      brat3dData[nBrat3dData].dx = dx;
      brat3dData[nBrat3dData].nx = nx;

      brat3dData[nBrat3dData].By = Byd;
      brat3dData[nBrat3dData].yi = yi;
      brat3dData[nBrat3dData].yf = yf;
      brat3dData[nBrat3dData].dy = dy;
      brat3dData[nBrat3dData].ny = ny;

      brat3dData[nBrat3dData].Bz = Bzd;
      brat3dData[nBrat3dData].zi = zi;
      brat3dData[nBrat3dData].zf = zf;
      brat3dData[nBrat3dData].dz = dz;
      brat3dData[nBrat3dData].nz = nz;
      cp_str(&brat3dData[nBrat3dData].filename, brat->filename);
      
      brat->dataIndex = nBrat3dData;
      nBrat3dData += 1;
    }
    brat->initialized = 1;
  }

  if (brat->dataIndex<0 || brat->dataIndex>=nBrat3dData)
    bombElegant("BRAT data indexing bug (1). Please report.", NULL);
  if (strcmp(brat->filename, brat3dData[brat->dataIndex].filename)!=0)
    bombElegant("BRAT data indexing bug (2). Please report.", NULL);
  
  BxNorm = brat3dData[brat->dataIndex].Bx;
  xi = brat3dData[brat->dataIndex].xi;
  xf = brat3dData[brat->dataIndex].xf;
  dx = brat3dData[brat->dataIndex].dx;
  nx = brat3dData[brat->dataIndex].nx;
  
  ByNorm = brat3dData[brat->dataIndex].By;
  yi = brat3dData[brat->dataIndex].yi;
  yf = brat3dData[brat->dataIndex].yf;
  dy = brat3dData[brat->dataIndex].dy;
  ny = brat3dData[brat->dataIndex].ny;
  
  BzNorm = brat3dData[brat->dataIndex].Bz;
  zi = brat3dData[brat->dataIndex].zi;
  zf = brat3dData[brat->dataIndex].zf;
  dz = brat3dData[brat->dataIndex].dz;
  nz = brat3dData[brat->dataIndex].nz;

  zStart = zi-dz;
  z_outer = MAX(fabs(zi), fabs(zf));

  fieldMapDimension = 3;

  fieldFactor = brat->fieldFactor;
  theta = brat->angle;
  fse = brat->fse;
  dXOffset = brat->dxMap;
  dZOffset = brat->dzMap;
  magnetYaw = brat->yawMap;
  
  xVertex = brat->xVertex;
  zVertex = brat->zVertex;
  xNomEntry = brat->xEntry;
  zNomEntry = brat->zEntry;
  xNomExit = brat->xExit;
  zNomExit = brat->zExit;
  zCenter = (zNomEntry+zNomExit)/2;
  xCenter = (xNomEntry+xNomExit)/2;
  central_length = 100*(zNomExit-zNomEntry);

  theta = refineAngle(theta, zNomEntry, xNomEntry, zVertex, xVertex,
                      zNomExit, xNomExit);
  
  rigidity = pCentral*particleMass*c_mks/particleCharge*particleRelSign;
  rhoMax = fabs(rigidity/brat3dData[brat->dataIndex].Bmax);
  
  integ_tol = brat->accuracy;
  zero_tol = integ_tol*10;
  useFTABLE = brat->useFTABLE;

  for (ip=0; ip<np; ip++) {
    double accelCoord[6], q[10];

    for (ic=0; ic<6; ic++)
      accelCoord[ic] = part[ip][ic];

    BRAT_lorentz_integration(accelCoord, q, 0);

    for (ic=0; ic<6; ic++)
      part[ip][ic] = accelCoord[ic];
  }

  return np;
}


double BRAT_optim_function(double *param, long *invalid)
{
  double q[9], result;
  double accelCoord[6] = {0,0,0,0,0,0};
  double *w;

  fse = param[0];
  dXOffset = param[1];
  dZOffset = param[2];
  magnetYaw = param[3];
  BRAT_lorentz_integration(accelCoord, q, 0);
  *invalid = 0;
  w = q+3;
  result = sqrt(sqr(accelCoord[0]) + sqr(accelCoord[1]));
  if (verbose_optimize && !quiet)
    printf("optim_function called: FSE = %e   dX = %e   w1/w0 = %e  sf = %e  penalty function = %e\n", 
           param[0], param[1], w[1]/w[0], fabs(accelCoord[4]), result); 
  return(result);
}

void BRAT_report_function(double result, double *param, long pass,
                     long n_evals, long n_dim)
{
  if (!quiet) {
    printf("report for pass %ld of optimization:\n", pass);
    printf("    %ld evaluations so far\n", n_evals);
    printf("    penalty function is %e\n", result);
    printf("    FSE      = %.15e\n", param[0]);
    printf("    dX       = %.15e\n", param[1]);
    printf("    dZ       = %.15e\n", param[2]);
    printf("    yaw      = %.15e\n", param[3]);
  }
}

#define N_EVAL_MAX 1500
#define N_PASS_MAX 3

void BRAT_optimize_magnet(unsigned long flags)
{
  double x[4], dx[4], xlo[4], xhi[4];
  short disable[4] = {0,0,0,0};
  double tolerance, result;
  long dummy;

  tolerance = 1e-10;
  x[0] = fse;
  x[1] = dXOffset;
  x[2] = dZOffset;
  x[3] = magnetYaw;
  dx[0] = dx[1] = dx[2] = dx[3] = 1e-4;
  xlo[0] = fseLimit[0]; 
  xhi[0] = fseLimit[1];
  xlo[1] = dxLimit[0]; 
  xhi[1] = dxLimit[1];
  xlo[2] = dzLimit[0]; 
  xhi[2] = dzLimit[1];
  xlo[3] = yawLimit[0]; 
  xhi[3] = yawLimit[1];
  if (!(flags&OPTIMIZE_FSE))
    disable[0] = 1;
  if (!(flags&OPTIMIZE_DX))
    disable[1] = 1;
  if (!(flags&OPTIMIZE_DZ))
    disable[2] = 1;
  if (!(flags&OPTIMIZE_YAW))
    disable[3] = 1;
  if (simplexMin(&result, x, dx, xlo, xhi, disable, 4, tolerance,
                 tolerance, BRAT_optim_function, BRAT_report_function,
                 N_EVAL_MAX, N_PASS_MAX, 12, 3.0, 1.0, 0)<0) {
    fprintf(stderr, "warning: optimization of magnet `failed'\n");
  }
  dx[0] = dx[1] = dx[2] = dx[3] = 1e-4;
  if (simplexMin(&result, x, dx, xlo, xhi, disable, 4, tolerance,
                 tolerance, BRAT_optim_function, BRAT_report_function,
                 N_EVAL_MAX, N_PASS_MAX, 12, 3.0, 1.0, 0)<0) {
    fprintf(stderr, "warning: optimization of magnet `failed'\n");
  }
  fse = x[0];
  dXOffset = x[1];
  dZOffset = x[2];
  magnetYaw = x[3];
  result = BRAT_optim_function(x, &dummy);
  if (!quiet) {
    printf("penalty function after optimization: %e\n", result);
    printf("FSE = %.15e\n", x[0]);
    printf("dX  = %.15e\n", x[1]);
    printf("dZ  = %.15e\n", x[2]);
    printf("Yaw = %.15e\n", x[3]);
  }
}

#define MIN_N_STEPS 100

static double global_delta;

void BRAT_lorentz_integration(
                            double *accelCoord, 
                            double *q,   /* z, x, y, dz/dS, dx/dS, dy/dS, dF/dS */
                            long doStoreData
                            )
{
  long n_eq = 10;
  double exvalue, qptest[10];
  long misses[10], accmode[10];
  double accuracy[10], tiny[10];
  double hrec, hmax, dSds;
  double s_start=0, s_end, exit_toler, ds, dz;
  long int_return;
  double *w, *IF;
  //double xStart, dx;
  double zStart, slope, phi;
  
  global_delta = accelCoord[5];
  particle_inside = 0;
  n_stored = 0;
  /* drift to parallel plane for start of integration */
  slope = (xVertex-xNomEntry)/(zVertex-zNomEntry);
  phi = atan(slope);
  if (idealMode) {
    dz = idealChord;
    //dx = slope*dz;
  } else {
    dz = zf-zi;
    //dx = slope*dz;
  }
  //xStart = xNomEntry - dx;
  zStart = zNomEntry - dz;

  /* Compute drift-back distance and apply to path length */
  /* ds = (zStart-(zNomEntry+accelCoord[0]*sin(phi)))/(cos(phi+atan(accelCoord[1]))*cos(atan(accelCoord[3]))); */
  /* New expression from R. Lindberg */
  ds = (zStart-(zNomEntry-accelCoord[0]*sin(phi)))*sqrt(1+sqr(accelCoord[1])+sqr(accelCoord[3]))/(cos(phi)-accelCoord[1]*sin(phi));
  accelCoord[4] += ds;

#ifdef DEBUG
  fprintf(stderr, "drift-back: %e\n", ds);
#endif

  /* transform from accel coords to cartesian coords */

  /* note that w0^2+w1^2+w2^2 = 1 */
  /* dS/ds is the rate of change of distance traveled w.r.t. distance of the fiducial particle */
  dSds = sqrt(1+sqr(accelCoord[1])+sqr(accelCoord[3]));
  w = q+3;
  w[0] = (-accelCoord[1]*sin(phi)+cos(phi))/dSds;
  w[1] = ( accelCoord[1]*cos(phi)+sin(phi))/dSds;
  w[2] = accelCoord[3]/dSds;

  q[0] = zNomEntry - accelCoord[0]*sin(phi) + ds*w[0];
  q[1] = xNomEntry + accelCoord[0]*cos(phi) + ds*w[1];
  q[2] = accelCoord[2] + ds*w[2];

#ifdef DEBUG
  fprintf(stderr, "initial: q[0] = %21.15e, q[1] = %21.15e, q[2] = %21.15e\n",
          q[0], q[1], q[2]);
#endif


  IF = q+6;
  IF[0] = IF[1] = IF[2] = IF[3] = 0;
  BRAT_deriv_function(qptest, q, 0.0L);
  if (qptest[3]!=0 || qptest[4]!=0 || qptest[5]!=0)
    bomb("particle started inside the magnet!", NULL);

  if (!useFTABLE) {
    fill_long_array(accmode, n_eq, 2);
    fill_double_array(tiny, n_eq, TOLERANCE_FACTOR);
    fill_double_array(accuracy, n_eq, integ_tol);
    fill_long_array(misses, n_eq, 0);
    
    /* use adaptive integration */
    s_start = 0;
    s_end   = central_length*2;
    hmax    = fabs(rhoMax*theta)/MIN_N_STEPS;
    hrec    = hmax/10;
    if (idealMode)
      zEnd = 2*idealChord;
    else
      zEnd = zf + (zf-zi);
    if ((exit_toler = sqr(integ_tol)*s_end)<central_length*TOLERANCE_FACTOR)
      exit_toler = central_length*TOLERANCE_FACTOR;
    switch (int_return = bs_odeint(q, BRAT_deriv_function, n_eq, accuracy,
                                   accmode, tiny, misses,
                                   &s_start, s_end, exit_toler, hrec, hmax, &hrec, BRAT_exit_function,
                                   exit_toler, 0, doStoreData?BRAT_store_data:NULL)) {
    case DIFFEQ_ZERO_STEPSIZE:
    case DIFFEQ_CANT_TAKE_STEP:
    case DIFFEQ_OUTSIDE_INTERVAL:
    case DIFFEQ_XI_GT_XF:
      printf("Integration failure---may be program bug\n");
      exit(1);
      break;
    case DIFFEQ_END_OF_INTERVAL:
      break;
    default:
      if ((exvalue = BRAT_exit_function(NULL, q, 0.0L))>exit_toler)
        bomb("warning: exit value of exceeds tolerance", NULL);
      break;
    }
#ifdef DEBUG
    fprintf(stderr, "final: q[0] = %21.15e, q[1] = %21.15e, q[2] = %21.15e\n",
            q[0], q[1], q[2]);
#endif
  } else {
    /* FTABLE mode */
#ifdef FTABLE_WORKING
    double eomc;
    double **A, step;
    long  nKicks;
    nKicks = useFTABLE;
    eomc = -particleCharge/particleMass/c_mks;
    step = (zf-zi)/nKicks;
    A = (double**)czarray_2d(sizeof(double), 3, 3);

    for (ik=0; ik<nKicks; ik++) {
      /* 1. get particle's coordinates */
      coord = particle[ip];
      factor = sqrt(1+sqr(coord[1])+sqr(coord[3]));
      p0 = (1.+coord[5])*Po;
      p[2] = p0/factor;
      p[0] = coord[1]*p[2];
      p[1] = coord[3]*p[2];
      
      /* 2. get field at the middle point */
      xyz[0] = coord[0] + p[0]/p[2]*step/2.0;
      xyz[1] = coord[2] + p[1]/p[2]*step/2.0;
      xyz[2] = s_location; 
      interpolateFTable(B, xyz, ftable);
      if (fpdebug) 
        fprintf(fpdebug, "%ld 0 %ld %21.15e %21.15e %21.15e %21.15e %21.15e %21.15e\n", ik, ip, xyz[0], xyz[1], xyz[2], B[0], B[1], B[2]);
      BA = sqrt(sqr(B[0]) + sqr(B[1]) + sqr(B[2]));
      /* 3. calculate the rotation matrix */
      A[0][0] = -(p[1]*B[2] - p[2]*B[1]);
      A[0][1] = -(p[2]*B[0] - p[0]*B[2]);
      A[0][2] = -(p[0]*B[1] - p[1]*B[0]);
      pA = sqrt(sqr(A[0][0]) + sqr(A[0][1]) + sqr(A[0][2]));
      /* When field not equal to zero or not parallel to the particles motion */
      if (BA>ftable->threshold && pA) {
        A[0][0] /= pA;
        A[0][1] /= pA;
        A[0][2] /= pA;
        A[1][0] = B[0]/BA;
        A[1][1] = B[1]/BA;
        A[1][2] = B[2]/BA;
        A[2][0] = A[0][1]*A[1][2]-A[0][2]*A[1][1];
        A[2][1] = A[0][2]*A[1][0]-A[0][0]*A[1][2];
        A[2][2] = A[0][0]*A[1][1]-A[0][1]*A[1][0];
        
        /* 4. rotate coordinates from (x,y,z) to (u,v,w) with u point to BxP, v point to B */
        pz0 = p[2];
        rotate_coordinate(A, p, 0);
        if (p[2] < 0)
          bombElegant("Table function doesn't support particle going backward", NULL);
        rotate_coordinate(A, B, 0);
        
        /* 5. apply kick */
        rho = p[2]/(eomc*B[1]);
        theta0=theta1=theta2=0.;
        if (A[2][2]) {
          tm_a =  3.0*A[0][2]/A[2][2];
          tm_b = -6.0*A[1][2]*p[1]/p[2]/A[2][2]-6.0;
          tm_c =  6.0*step/rho/A[2][2];
#ifdef USE_GSL
          gsl_poly_solve_cubic (tm_a, tm_b, tm_c, &theta0, &theta1, &theta2);
#else
          bombElegant("gsl_poly_solve_cubic function is not available becuase this version of elegant was not built against the gsl library", NULL);
#endif
        } else if (A[0][2]) {
          tm_a = A[1][2]*p[1]/p[2]+A[2][2];
          theta0 = (tm_a-sqrt(sqr(tm_a)-2.*A[0][2]*step/rho))/A[0][2];
          theta1 = (tm_a+sqrt(sqr(tm_a)-2.*A[0][2]*step/rho))/A[0][2];
        } else {
          tm_a = A[1][2]*p[1]/p[2]+A[2][2];          
          theta0 = step/rho/tm_a;
        }
        theta=choose_theta(rho, theta0, theta1, theta2);
        if (fpdebug) 
          fprintf(fpdebug, "%ld 1 %ld %21.15e %21.15e %21.15e %21.15e %21.15e %21.15e\n", ik, ip, xyz[0], xyz[1], xyz[2], B[0], B[1], B[2]);
        
        p[0] = -p[2]*sin(theta);
        p[2] *= cos(theta);
        xyz[0] = rho*(cos(theta)-1);
        xyz[1] = (p[1]/p[2])*rho*theta;
        xyz[2] = rho*sin(theta);
        
        /* 6. rotate back to (x,y,z) */
        rotate_coordinate(A, xyz, 1);
        rotate_coordinate(A, p, 1);
        coord[0] += xyz[0];
        coord[2] += xyz[1];
        coord[4] += sqrt(sqr(rho*theta)+sqr(xyz[1]));
        coord[1] = p[0]/p[2];
        coord[3] = p[1]/p[2];
      } else {
        coord[0] += coord[1]*step;
        coord[2] += coord[3]*step;
        coord[4] += step*factor;         
      }
      s_location += step;
    }
#endif
    /* convert back to (Z, X, Y, WZ, WX, WY) */
  }
  
  /* drift back to reference plane */
  slope = (xVertex-xNomExit)/(zVertex-zNomExit);
  phi = -atan(slope);
  ds = ((zNomExit-q[0])*cos(phi) - (xNomExit-q[1])*sin(phi))/(w[0]*cos(phi) - w[1]*sin(phi));
#ifdef DEBUG
  printf("exit: phi = %le, ds = %le\n", phi, ds);
#endif
  q[0] += ds*w[0];
  q[1] += ds*w[1];
  q[2] += ds*w[2];
  
  /* convert back to accelerator coordinates */
  dSds = 1./(-w[1]*sin(phi)+w[0]*cos(phi));
  
  accelCoord[1] = (w[1]*cos(phi)+w[0]*sin(phi))*dSds;
  accelCoord[3] = w[2]*dSds;
  accelCoord[0] = (q[0]-zNomExit)/sin(phi);
  accelCoord[2] = q[2];
  accelCoord[4] += s_start+ds;
#ifdef DEBUG
  fprintf(stderr, "final: x=%e, xp=%e, y=%e, yp=%e\n",
          accelCoord[0], accelCoord[1], accelCoord[2], accelCoord[3]);
#endif

#ifdef DEBUG
  fprintf(stderr, "path length: %e\n", s_start);
  fprintf(stderr, "drift back exit: %e\n", ds);
  fprintf(stderr, "ref plane: x=%e, xp=%e, y=%e, yp=%e, s=%e\n",
          accelCoord[0], accelCoord[1], accelCoord[2], accelCoord[3], accelCoord[4]);
#endif
}

#undef DEBUG

double BRAT_exit_function(double *qp, double *q, double s)
{
  return q[0] - zEnd;
}

void BRAT_deriv_function(double *qp, double *q, double s)
{
  double *F;
  static double *w, *wp, factor;

  F = qp+6;
  BRAT_B_field(F, q);
  F[3] = (1-F[2])*F[2];
  w  = q+3;
  wp = qp+3;
  qp[0] = w[0];
  qp[1] = w[1];
  qp[2] = w[2];
  factor = fieldFactor/(rigidity*(1+global_delta));
  wp[0] = (w[2]*F[1]-w[1]*F[2])*factor;
  wp[1] = (w[0]*F[2]-w[2]*F[0])*factor;
  wp[2] = (w[1]*F[0]-w[0]*F[1])*factor;
}

void BRAT_store_data(double *qp, double *q, double s, double exval)
{
  if (fabs(q[0])>2*z_outer)
    return;
  if (n_stored>=max_store) {
    max_store += 100;
    X_stored = trealloc(X_stored, sizeof(*X_stored)*max_store);
    Y_stored = trealloc(Y_stored, sizeof(*Y_stored)*max_store);
    Z_stored = trealloc(Z_stored, sizeof(*Z_stored)*max_store);
    wX_stored = trealloc(wX_stored, sizeof(*X_stored)*max_store);
    wY_stored = trealloc(wY_stored, sizeof(*Y_stored)*max_store);
    wZ_stored = trealloc(wZ_stored, sizeof(*Z_stored)*max_store);
    aX_stored = trealloc(aX_stored, sizeof(*X_stored)*max_store);
    aY_stored = trealloc(aY_stored, sizeof(*Y_stored)*max_store);
    aZ_stored = trealloc(aZ_stored, sizeof(*Z_stored)*max_store);
    FX_stored = trealloc(FX_stored, sizeof(*FX_stored)*max_store);
    FY_stored = trealloc(FY_stored, sizeof(*FY_stored)*max_store);
    FZ_stored = trealloc(FZ_stored, sizeof(*FZ_stored)*max_store);
    Fint_stored = trealloc(Fint_stored, sizeof(*Fint_stored)*max_store);
    s_stored = trealloc(s_stored, sizeof(*s_stored)*max_store);
  }
  Z_stored[n_stored] = q[0];
  X_stored[n_stored] = q[1];
  Y_stored[n_stored] = q[2];
  wZ_stored[n_stored] = q[3];
  wX_stored[n_stored] = q[4];
  wY_stored[n_stored] = q[5];
  aZ_stored[n_stored] = qp[3];
  aX_stored[n_stored] = qp[4];
  aY_stored[n_stored] = qp[5];
  FZ_stored[n_stored] = qp[6];
  FX_stored[n_stored] = qp[7];
  FY_stored[n_stored] = qp[8];
  Fint_stored[n_stored] = q[9]/gap/2;
  s_stored[n_stored] = s;
  n_stored++;
}

double BRAT_setup_arc_field_data(char *input, char *sName, char *fieldName, double xCenter)
{
  SDDS_DATASET SDDS_table;
  long iData, increasing;
  double BRef, BMin, BMax, sMin, sMax;
  
  if (!SDDS_InitializeInput(&SDDS_table, input) || !SDDS_ReadPage(&SDDS_table)) {
    SDDS_SetError("Unable to read data file");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!(nArcPoints=SDDS_CountRowsOfInterest(&SDDS_table)))
    bomb("no data in field file", NULL);
  if (nArcPoints<2)
    bomb("too little data in field file", NULL);
  
  if (SDDS_CheckColumn(&SDDS_table, sName, "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY) 
    exit(1);
  if (SDDS_CheckColumn(&SDDS_table, fieldName, "T", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY)
    exit(1);

  if (!(sArc = SDDS_GetColumnInDoubles(&SDDS_table, sName)) ||
      !(BArcNorm = SDDS_GetColumnInDoubles(&SDDS_table, fieldName)))  {
    SDDS_SetError("Unable to retrieve data columns");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!(dBArcNormds = malloc(sizeof(*dBArcNormds)*nArcPoints))) {
    fprintf(stderr, "brat: memory allocation failure\n");
    exit(1);
  }
  
  increasing = 0;
  if (sArc[0]<sArc[1])
    increasing = 1;
  BMin = -(BMax = -DBL_MAX);
  sMin = -(sMax = -DBL_MAX);
  for (iData=0; iData<nArcPoints; iData++) {
    if (BArcNorm[iData]>BMax)
      BMax = BArcNorm[iData];
    if (BArcNorm[iData]<BMin)
      BMin = BArcNorm[iData];
    if (sArc[iData]>sMax)
      sMax = sArc[iData];
    if (sArc[iData]<sMin)
      sMin = sArc[iData];
    if (iData==0)
      continue;
    if ((increasing  && sArc[iData-1]>=sArc[iData]) ||
        (!increasing && sArc[iData-1]<=sArc[iData]))
      bomb("arc s values not monotonic", NULL);
  }
  
  if (fabs(BMin)>fabs(BMax))
    BRef = BMin;
  else
    BRef = BMax;
  if (sMin!=0) {
    fprintf(stderr, "minimum value of s is %le, but must be zero for arc data\n",
            sMin);
    exit(1);
  }

  if (BRef==0)
    bomb("field appears to be zero in arc data", NULL);
  /*
  for (iData=0; iData<nArcPoints; iData++)
    BArcNorm[iData] /= BRef;
  */

  dBArcNormds[0] = (BArcNorm[1]-BArcNorm[0])/(sArc[1]-sArc[0]);
  dBArcNormds[nArcPoints-1] 
    = (BArcNorm[nArcPoints-1]-BArcNorm[nArcPoints-2])
      /(sArc[nArcPoints-1]-sArc[nArcPoints-2]);

  for (iData=1; iData<(nArcPoints-1); iData++) {
    dBArcNormds[iData] 
      = (BArcNorm[iData+1]-BArcNorm[iData-1])
        /(sArc[iData+1]-sArc[iData-1]);
  }
  
  return BRef;
}

double BRAT_setup_field_data(char *input, double xCenter, double zCenter)
{
  SDDS_TABLE SDDS_table;
  long idata, ix, iz, rows;
  double x, z, B, *xd=NULL, *zd=NULL, *Bd=NULL, xCheck, zCheck;
  double *yd = NULL, *Bxd = NULL, *Byd = NULL, *Bzd = NULL;
  double Bmin, Bmax;
  long ix0, ix1, iz0, iz1;
  double dz0=0, dx0=0;
#ifdef DEBUG
  double xc;
#endif

  if (!SDDS_InitializeInput(&SDDS_table, input) || !SDDS_ReadPage(&SDDS_table)) {
    SDDS_SetError("Unable to read data file");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!(rows=SDDS_CountRowsOfInterest(&SDDS_table)))
    bomb("no data in field file", NULL);

  if (fieldMapDimension==2) {
    if (!single_scan) {
      if (SDDS_CheckColumn(&SDDS_table, "x", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY || 
          SDDS_CheckColumn(&SDDS_table, "z", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY ||
          SDDS_CheckColumn(&SDDS_table, "B", "T", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY) {
        exit(1);
      }
      if (!(xd=(double*)SDDS_GetColumnInDoubles(&SDDS_table, "x")) ||
          !(zd=(double*)SDDS_GetColumnInDoubles(&SDDS_table, "z")) ||
          !((Bd=(double*)SDDS_GetColumnInDoubles(&SDDS_table, "B")))) {
        SDDS_SetError("Unable to retrieve data column(s)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      
      if (SDDS_CheckParameter(&SDDS_table, "xMinimum", "m", SDDS_DOUBLE, stderr)!=SDDS_CHECK_OKAY || 
          SDDS_CheckParameter(&SDDS_table, "xInterval", "m", SDDS_DOUBLE, stderr)!=SDDS_CHECK_OKAY ||
          SDDS_CheckParameter(&SDDS_table, "xDimension", NULL, SDDS_LONG, stderr)!=SDDS_CHECK_OKAY ||
          SDDS_CheckParameter(&SDDS_table, "zMinimum", "m", SDDS_DOUBLE, stderr)!=SDDS_CHECK_OKAY || 
          SDDS_CheckParameter(&SDDS_table, "zInterval", "m", SDDS_DOUBLE, stderr)!=SDDS_CHECK_OKAY ||
          SDDS_CheckParameter(&SDDS_table, "zDimension", NULL, SDDS_LONG, stderr)!=SDDS_CHECK_OKAY)
        exit(1);
      if (!SDDS_GetParameter(&SDDS_table, "xMinimum", &xi) || !SDDS_GetParameter(&SDDS_table, "xDimension", &nx) ||
          !SDDS_GetParameter(&SDDS_table, "xInterval", &dx) ||
          !SDDS_GetParameter(&SDDS_table, "zMinimum", &zi) || !SDDS_GetParameter(&SDDS_table, "zDimension", &nz) ||
          !SDDS_GetParameter(&SDDS_table, "zInterval", &dz)) {
        SDDS_SetError("Unable to retrieve parameter value(s)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (nx<=0 || nz<=0) 
        bomb("non-positive value for dx and/or dz", NULL);
      xf = (nx-1)*dx+xi;
      zf = (nz-1)*dz+zi;
      if (nx==1) {
        single_scan = 1;
        xi = xf = 0;
        if (!quiet)
          printf("data treated as single scan\n");
      }
      else if (!quiet)
        printf("data grid:\n    x: [%.8f, %.8f]m   dx = %.8fm   nx = %ld\n",
               xi, xf, dx, nx);
      if (!quiet)
      printf("    z: [%.8f, %.8f]m   dz = %.8fm   nz = %ld\n",
             zi, zf, dz, nz);
    }
    else {
      if (SDDS_CheckColumn(&SDDS_table, "z", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY) 
        exit(1);
      if (SDDS_CheckColumn(&SDDS_table, "B", "T", SDDS_ANY_FLOATING_TYPE, stderr)==SDDS_CHECK_OKAY)
        exit(1);
      if (!(zd=(double*)SDDS_GetColumnInDoubles(&SDDS_table, "z"))) {
        SDDS_SetError("Unable to retrieve data z column");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (!(Bd=(double*)SDDS_GetColumnInDoubles(&SDDS_table, "B"))) {
        SDDS_SetError("Unable to retrieve data field column");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      nx = 1;
      xi = xf = 0;
      zi = zd[0];
      zf = zd[rows-1];
      nz = rows;
      if (nz<1)
        bomb("nz < 1 in single scan data", NULL);
      dz = (zf-zi)/(nz-1);
      if (!quiet) {
        printf("data grid:\n    ");
        printf("    z: [%.8f, %.8f]m   dz = %.8fm   nz = %ld\n",
               zi, zf, dz, nz);
      }
    }

#ifdef DEBUG
    fprintf(stderr, "Allocating %ld x %ld field and gradient arrays\n", 
            nx, nz);
#endif
    Bnorm = (double**)zarray_2d(sizeof(**Bnorm), nx, nz);
    dBnormdz = (double**)zarray_2d(sizeof(**dBnormdz), nx, nz);
    dBnormdx = (double**)zarray_2d(sizeof(**dBnormdx), nx, nz);
#ifdef DEBUG
    fprintf(stderr, "Initializing Bnorm array\n");
#endif
    for (ix=0; ix<nx; ix++)
      for (iz=0; iz<nz; iz++)
        Bnorm[ix][iz] = DBL_MAX;
    idata = 0;
#ifdef DEBUG
    xc = xi;
#endif
    Bmin = -(Bmax = -DBL_MAX);
    ix = 0;
    if (single_scan)
      dx = 1;
    for (idata=0; idata<rows; idata++) {
#ifdef DEBUG
      fprintf(stderr, "Checking idata=%ld\n", idata);
#endif
      x = single_scan?0:xd[idata];
      z = zd[idata];
      B = Bd[idata];
#ifdef DEBUG
      fprintf(stderr, "x=%e, z=%e, B=%e\n", x, z, B);
#endif
      if (!single_scan)
        ix = (x-xi)/dx+0.5;
      iz = (z-zi)/dz+0.5;
      xCheck = ix*dx+xi;
      zCheck = iz*dz+zi;
#ifdef DEBUG
      fprintf(stderr, "ix = %ld, iz = %ld, xc=%e, zc=%e\n", ix, iz,
              xCheck, zCheck);
#endif
      if (fabs(x-xCheck)/dx>1e-2 || fabs(z-zCheck)/dz>1e-2) {
        fprintf(stderr, "Grid data not uniform: (%e, %e) -> (%e, %e)\n",
                x, z, xCheck, zCheck);
        exit(1);
      }
      if (ix<0 || ix>nx || iz<0 || iz>nz) {
        fprintf(stderr, "Point outside grid: (%e, %e) -> (%ld, %ld)\n",
                x, z, ix, iz);
        exit(1);
      }
      if (Bnorm[ix][iz]!=DBL_MAX)
        bomb("two data points map to one grid point", NULL);
      Bnorm[ix][iz] = B;
      if (B<Bmin)
        Bmin = B;
      if (B>Bmax)
        Bmax = B;
    }
    if (!single_scan) 
      free(xd);
    free(zd);
    free(Bd);
    if (!SDDS_Terminate(&SDDS_table))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    
    idata = 0;
    for (ix=0; ix<nx; ix++)
      for (iz=0; iz<nz; iz++) {
        if (Bnorm[ix][iz]==DBL_MAX) {
          idata ++;
          Bnorm[ix][iz] = 0;
        }
      }
    if (idata && !quiet)
      fprintf(stderr, "Warning: no data for %ld of %ld points\n",
              idata, nx*nz);
    
    if (fabs(Bmin)>fabs(Bmax))
      Bmax = Bmin;
    if (Bmax==0)
      bomb("extremal B value is zero", NULL);
    if (!quiet) {
      printf("Maximum value of B is %.8f T\n", Bmax);
    }
    /* 
      if (!quiet)
      printf("reference B value is %.8f T\n", Bmax);
      for (ix=0; ix<nx; ix++)
      for (iz=0; iz<nz; iz++)
      Bnorm[ix][iz] /= Bmax;
      */

    /* compute derivatives */
    for (ix=0; ix<nx; ix++) {
      if (nx==1)
        ix0 = ix1 = 0;
      else {
        if (ix==0) {
          dx0 = dx;
          ix1 = ix+1;
          ix0 = ix;
        }
        else if (ix==(nx-1)) {
          dx0 = dx;
          ix1 = ix;
          ix0 = ix-1;
        }
        else {
          dx0 = 2*dx;
          ix1 = ix+1;
          ix0 = ix-1;
        }
      }
      for (iz=0; iz<nz; iz++) {
        if (iz==0) {
          dz0 = dz;
          iz1 = iz+1;
          iz0 = iz;
        }
        else if (iz==(nz-1)) {
          dz0 = dz;
          iz1 = iz;
          iz0 = iz-1;
        }
        else {
          dz0 = 2*dz;
          iz1 = iz+1;
          iz0 = iz-1;
        }
        dBnormdz[ix][iz] = (Bnorm[ix][iz1]-Bnorm[ix][iz0])/dz0;
        dBnormdx[ix][iz] = (Bnorm[ix1][iz]-Bnorm[ix0][iz])/dx0;
      }
    }
  } else {
    /* 3d field map */
    if (SDDS_CheckColumn(&SDDS_table, "x", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY || 
        SDDS_CheckColumn(&SDDS_table, "y", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY ||
        SDDS_CheckColumn(&SDDS_table, "z", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY) {
      exit(1);
    }
    if (SDDS_CheckColumn(&SDDS_table, "Bx", "T", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY || 
        SDDS_CheckColumn(&SDDS_table, "By", "T", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY ||
        SDDS_CheckColumn(&SDDS_table, "Bz", "T", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY) {
      exit(1);
    }
    if (!(xd = SDDS_GetColumnInDoubles(&SDDS_table, "x")) ||
        !(yd = SDDS_GetColumnInDoubles(&SDDS_table, "y")) ||
        !(zd = SDDS_GetColumnInDoubles(&SDDS_table, "z")) ||
        !(Bxd = SDDS_GetColumnInDoubles(&SDDS_table, "Bx")) ||
        !(Byd = SDDS_GetColumnInDoubles(&SDDS_table, "By")) ||
        !(Bzd = SDDS_GetColumnInDoubles(&SDDS_table, "Bz"))) {
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }

    /* It is assumed that the data is ordered so that x changes fastest.
     * This can be accomplished with sddssort -column=z,incr -column=y,incr -column=x,incr
     * The points are assumed to be equipspaced.
     */
    nx = 1;
    xi = xd[0];
    while (nx<rows) {
      if (xd[nx-1]>xd[nx])
        break;
      nx ++;
    }
    if (nx==rows) {
      fprintf(stderr, "file doesn't have correct structure or amount of data (x)\n");
      fprintf(stderr, "Use sddssort -column=z -column=y -column=x to sort the file\n");
      exit(1);
    }  
    xf = xd[nx-1];
    dx = (xf-xi)/(nx-1);

    ny = 1;
    yi = yd[0];
    while (ny<(rows/nx)) {
      if (yd[(ny-1)*nx]>yd[ny*nx])
        break;
      ny++;
    }
    if (ny==rows) {
      fprintf(stderr, "file doesn't have correct structure or amount of data (y)\n");
      fprintf(stderr, "Use sddssort -column=z -column=y -column=x to sort the file\n");
      exit(1);
    }
    yf = yd[(ny-1)*nx];
    dy = (yf-yi)/(ny-1);

    if (nx<=1 || ny<=1 || (nz = rows/(nx*ny))<=1) {
      fprintf(stderr, "file doesn't have correct structure or amount of data (z)\n");
      fprintf(stderr, "Use sddssort -column=z -column=y -column=x to sort the file\n");
      exit(1);
    }
    zi = zd[0];
    zf = zd[rows-1];
    dz = (zf-zi)/(nz-1);

    Bmin = -(Bmax = -DBL_MAX);
    for (idata=0; idata<rows; idata++) {
      if (fabs(yd[idata])<dy/2 && fabs(xd[idata])<dx/2) {
        if (Byd[idata]>Bmax)
          Bmax = Byd[idata];
        if (Byd[idata]<Bmin)
          Bmin = Byd[idata];
      }
    }
    free(xd);
    free(yd);
    free(zd);
    if (Bmax==-DBL_MAX) {
      fprintf(stderr, "file doesn't have valid Bmax value---verify that values with x=y=0 exist\n");
      exit(1);
    }
    if (fabs(Bmin)>fabs(Bmax))
      SWAP_DOUBLE(Bmin, Bmax);
    BxNorm = Bxd;
    ByNorm = Byd;
    BzNorm = Bzd;
    
    if (!quiet) 
      printf("3D field map data: nx=%ld, ny=%ld, nz=%ld\ndx=%e, dy=%e, dz=%e\nx:[%e, %e], y:[%e, %e], z:[%e, %e]\nBy:[%e, %e]\n",
              nx, ny, nz, dx, dy, dz,
              xi, xf, yi, yf, zi, zf,
              Bmin, Bmax
              );
  }

  return Bmax;
}

void BRAT_B_field(double *F, double *Qg)
{
  long ix, iy, iz, j;
  double fx, fy, fz, val_z1, val_z2, f[3];
  double x, y, z, derivSign=1;
  double Q[3];

  memcpy(Q, Qg, 3*sizeof(*Q));
  if (idealMode) {
    double Lx;
    Lx = idealChord + 2*Qg[1]*tan(idealEdgeAngle*PI/180);
    if (fabs(Qg[0])>Lx/2) {
      F[0] = F[1] = F[2] = 0;
    } else {
      F[0] = F[1] = 0;
      F[2] = 1;
    }
    for (j=0; j<3; j++)
      F[j] *= fieldSign*(1+fse);
    return;
  }
  if (arc_scan) {
    double phi, s, dByds;
    //double R;
    long inStraight;
    OUTRANGE_CONTROL belowRange = {0.0,  OUTRANGE_SATURATE } ;
    OUTRANGE_CONTROL aboveRange = {0.0,  OUTRANGE_VALUE } ;
    unsigned long code;
    phi = fabs(atan2(Q[0], Q[1]));
#ifdef DEBUG
    fprintf(stderr, "computing arc field for Z=%e, X=%e, phi=%e: ",
            Q[0], Q[1], phi);
#endif
    if (phi<=theta/2) {
      s = phi*rhoArc;
      inStraight = 0;
    } else {
      //R = sqrt(sqr(Q[0])+sqr(Q[1]));
      s = rhoArc*((theta/2) + tan(phi-theta/2));
      inStraight = 1;
    }
    F[2] = interpolate(BArcNorm, sArc, nArcPoints, s, &belowRange, &aboveRange, 1, &code, 1.0);
    dByds = interpolate(dBArcNormds, sArc, nArcPoints, s, &belowRange, &aboveRange, 1, &code, 1.0);
    if (inStraight) {
      phi = theta/2;
    }
    F[0] = SIGN(Q[0])*Q[2]*dByds*cos(phi);
    F[1] = -Q[2]*dByds*sin(phi);
#ifdef DEBUG
    fprintf(stderr, "s=%e, B=%e\n", s, F[2]);
#endif
    for (j=0; j<3; j++)
      F[j] *= fieldSign*(1+fse);
    return ;
  }

  if (magnetYaw) {
    double Qy[2];
    Qy[0] = Q[0];
    Qy[1] = Q[1];
    Q[0] =  Qy[0]*cos(magnetYaw) + Qy[1]*sin(magnetYaw);
    Q[1] = -Qy[0]*sin(magnetYaw) + Qy[1]*cos(magnetYaw);
  }
  x = Q[1] - dXOffset;
  y = Q[2];

  z = Q[0];
  if (zDuplicate && z>0) {
    z = -fabs(z);
    derivSign = -1;
  }
  z -= dZOffset;
  
  F[0] = F[1] = F[2] = 0;

#ifdef DEBUG
  fprintf(stderr, "Finding field at x=%e (%e) z=%e\n", x, Q[1], z);
#endif

  if (!single_scan) {
    if (!particle_inside) {
      if ((xi-x)>dx*1e-6 || (x-xf)>dx*1e-6)
        return;
      if ((zi-z)>dz*1e-6 || (z-zf)>dz*1e-6)
        return;
      if (fieldMapDimension==3 &&
          ((yi-y)>dy*1e-6 || (y-yf)>dy*1e-6))
        return;
    }

    particle_inside = 1;
    if (!extend_data) {
      if ((xi-x)>dx*1e-6 || (x-xf)>dx*1e-6)
        return;
      if ((zi-z)>dz*1e-6 || (z-zf)>dz*1e-6)
        return;
      if (fieldMapDimension==3 &&
          ((yi-y)>dy*1e-6 || (y-yf)>dy*1e-6))
        return;
    }
  }
  else {
    if (extend_data && extend_edge_angle)
      z = z-x*tan(extend_edge_angle);
    if ((zi-z)>dz*1e-6 || (z-zf)>dz*1e-6)
      return;
  }

  iz = (z-zi)/dz;
  if (iz>=nz-1) {
    iz = nz-2;
    z = zf;
  }
  if (iz<0) {
    iz = 0;
    z = zi;
  }
  fz = (z-zi)/dz-iz;

  if (single_scan) {
    F[2] = Bnorm[0][iz] + fz*(Bnorm[0][iz+1]-Bnorm[0][iz]);
    if (y) {
      F[0] = derivSign*y*(dBnormdz[0][iz] + fz*(dBnormdz[0][iz+1]-dBnormdz[0][iz]));
      F[1] = y*(dBnormdx[0][iz] + fz*(dBnormdx[0][iz+1]-dBnormdx[0][iz]));
    }
    for (j=0; j<3; j++)
      F[j] *= fieldSign*(1+fse);
    return;
  }

  ix = (x-xi)/dx;
  if (ix>=nx-1) {
    ix = nx-2;
    x  = xf;
  }
  if (ix<0) {
    ix = 0;
    x = xi;
  }
  fx = (x-xi)/dx-ix;

  if (fieldMapDimension==3) {
    iy = (y-yi)/dy;
    if (iy>=ny-1) {
      iy = ny-2;
      y = yf;
    }
    if (iy<0) {
      iy = 0;
      y = yi;
    }
    fy = (y-yi)/dy-iy;
  } else {
    iy = 0;
    fy = 0;
  }
  
  if (fieldMapDimension==2) {
    val_z1 = Bnorm[ix][iz] + fx*(Bnorm[ix+1][iz]-Bnorm[ix][iz]);
    val_z2 = Bnorm[ix][iz+1] + fx*(Bnorm[ix+1][iz+1]-Bnorm[ix][iz+1]);
    F[2] = val_z1 + fz*(val_z2-val_z1);
    val_z1 = derivSign*(dBnormdz[ix][iz] + fx*(dBnormdz[ix+1][iz]-dBnormdz[ix][iz]));
    val_z2 = derivSign*(dBnormdz[ix][iz+1] + fx*(dBnormdz[ix+1][iz+1]-dBnormdz[ix][iz+1]));
    f[0] = y*(val_z1 + fz*(val_z2-val_z1));
    val_z1 = dBnormdx[ix][iz] + fx*(dBnormdx[ix+1][iz]-dBnormdx[ix][iz]);
    val_z2 = dBnormdx[ix][iz+1] + fx*(dBnormdx[ix+1][iz+1]-dBnormdx[ix][iz+1]);
    f[1] = y*(val_z1 + fz*(val_z2-val_z1));
    F[0] =  f[0];
    F[1] =  f[1];
  } else {
    long iq;
    double Finterp1[2][2], Finterp2[2];
    double *Fq[3], Freturn[3];
    Fq[0] = BzNorm;
    Fq[1] = BxNorm;
    Fq[2] = ByNorm;

#ifdef DEBUG
    fprintf(stderr, "Doing 3d interpolation: x=(%le, %le, %le), i=(%ld, %ld, %ld), f=(%le, %le, %le)\n", 
            x, y, z, ix, iy, iz, fx, fy, fz);
#endif

    for (iq=0; iq<3; iq++) {
      /* interpolate vs z to get four points in a x-y grid */
      /* (ix, iy) */
      Finterp1[0][0] = (1-fz)*(*(Fq[iq]+(ix+0)+iy*nx + iz*nx*ny)) +
        fz*(*(Fq[iq]+(ix+0)+iy*nx + (iz+1)*nx*ny));
      /* (ix+1, iy) */        
      Finterp1[1][0] = (1-fz)*(*(Fq[iq]+(ix+1)+(iy+0)*nx + iz*nx*ny)) +
        fz*(*(Fq[iq]+(ix+1)+(iy+0)*nx + (iz+1)*nx*ny));
      /* (ix, iy+1) */
      Finterp1[0][1] = (1-fz)*(*(Fq[iq]+(ix+0)+(iy+1)*nx + iz*nx*ny)) +
        fz*(*(Fq[iq]+(ix+0)+(iy+1)*nx + (iz+1)*nx*ny));
      /* (ix+1, iy+1) */
      Finterp1[1][1] = (1-fz)*(*(Fq[iq]+(ix+1)+(iy+1)*nx + iz*nx*ny)) +
        fz*(*(Fq[iq]+(ix+1)+(iy+1)*nx + (iz+1)*nx*ny));
      
      /* interpolate vs y to get two points spaced by dx */
      /* ix */
      Finterp2[0] = (1-fy)*Finterp1[0][0] + fy*Finterp1[0][1];
      /* ix+1 */
      Finterp2[1] = (1-fy)*Finterp1[1][0] + fy*Finterp1[1][1];
    
      /* interpolate vs x to get 1 value */
      Freturn[iq] = (1-fx)*Finterp2[0] + fx*Finterp2[1];
    }

    F[0] = Freturn[0];
    F[1] = Freturn[1];
    F[2] = Freturn[2];
  }
  
  for (j=0; j<3; j++)
    F[j] *= fieldSign*(1+fse);

#ifdef DEBUG
  fprintf(stderr, "F[0] = %e, F[1] = %e, F[2] = %e, FSE = %e\n", F[0], F[1], F[2], fse);
#endif

  return;
}

double refineAngle(double theta, 
                   double z0, double x0,
                   double zv, double xv, 
                   double z1, double x1)
{
  double dot, mag0, mag1;
  double cos_theta, theta1;
  
  dot = (xv-x0)*(x1-xv) + (zv-z0)*(z1-zv);
  mag0 = sqrt(sqr(xv-x0) + sqr(zv-z0));
  mag1 = sqrt(sqr(xv-x1) + sqr(zv-z1));
  cos_theta = dot/(mag0*mag1);

  theta1 = acos(cos_theta);
  if (fabs(theta1-theta)<1e-2)
    return theta1;

  theta1 = PIx2-theta1;
  if (fabs(theta1-theta)<1e-2)
    return theta1;
  
  theta1 = acos(cos_theta);
  printf("Error: failed to figure out refined angle: theta = %le, theta1 = %le\n",
         theta, theta1);
  exit(1);
}

