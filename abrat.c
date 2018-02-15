/* Copyright 2015 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* program: brat.c
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
#include "scan.h"
#include "table.h"
#include "SDDS.h"

#define SCAN_X 0
#define SCAN_XP 1
#define SCAN_Y 2
#define SCAN_YP 3
#define SCAN_DELTA 4
static char *variable_name[5] = {
    "x", "xp", "y", "yp", "delta"
    };
static char *initial_variable_name[5] = {
    "xi", "xpi", "yi", "ypi" , "deltai"
    };
static char *initial_variable_symbol[5] = {
    "x$bi$n", "x$a'$n$bi$n",
    "y$bi$n", "y$a'$n$bi$n",
    "$gd$r$bi$n"
    }  ;
static char *variable_units[5] = {
    "m", "", "m", "", ""
    } ;
static char *variable_description[5] = {
    "Initial x", "Initial x'",
    "Initial y", "Initial y'", "Initial $gd$r"
    } ;

#define SET_TOLERANCE 0
#define SET_THETA 1
#define SET_OUTPUT 2
#define SET_SCAN 3
#define SET_SINGLE_SCAN 4
#define SET_OPTIMIZE 5
#define SET_FIELDMAP_OUTPUT 6
#define SET_EXTEND_DATA 7
#define SET_FIELDSIGN 8
#define SET_VERTEX 9
#define SET_GAP 10
#define SET_BEAMFILES 11
#define SET_QUIET 12
#define SET_FSC 13
#define SET_DX_OFFSET 14
#define SET_Z_DUPLICATE 15
#define SET_DZ_OFFSET 16
#define SET_ARC_SCAN 17
#define SET_ENTRY 18
#define SET_EXIT 19
#define SET_RIGIDITY 20
#define SET_IDEAL 21
#define SET_FSE_LIMIT 22
#define SET_DX_LIMIT 23
#define SET_DZ_LIMIT 24
#define SET_YAW_LIMIT 25
#define SET_YAW 26
#define SET_3DFIELDMAP 27
#define N_OPTIONS 28

static char *option[N_OPTIONS] = {
    "tolerance", "theta", "output",
    "scan", "singlescan", "optimize",
    "fieldmapoutput", "extenddata", "fieldsign", "vertex",
    "gap", "beamfiles", "quiet", "fsc", 
    "dxdipole", "zduplicate", "dzdipole", "arcscan",
    "nominalEntrance", "nominalExit",
    "rigidity", "ideal", "fseLimit", "dxLimit", "dzLimit", "yawLimit", "yawDipole",
    "3dfieldfile",
    } ;

static char *USAGE = "abrat {<field-file>|-ideal=<fieldInTesla>,<chordInMeters>,<edgeAngleInDeg>} \n"
" [-3dFieldFile] [-zDuplicate] [-extendData[=edge-angle]] [-fieldSign={+|-}]\n"
" [{-scan={x | xp | y | yp | delta},lower,upper,number | \n"
" -beamFiles=<input>,<output> }]\n"
" -vertex=<x-in-meters>,<z-in-meters> -nominalEntrance=<x>,<y> -nominalExit=<x>,<y>\n"
" -theta=<targetInDegrees> -rigidity=<Tesla-meters>\n"
" [-output=filename] [-singleScan | -arcScan=<sName>,<fieldName>,<rhoIdeal>]\n"
" [-fsc=<value>] [-dxDipole=<m>] [-dzDipole=<m>] [-yawDipole=<value>]\n"
" {[-optimize[=verbose][{fse,dx,dz,yaw}]]\n"
"  -fseLimit=<min>,<max> -dxLimit=<min>,<max> -dzLimit=<min>,<max> -yawLimit=<min>,<max>}\n"
" [-fieldmapOutput=filename,zmin,zmax,nz,xmin,xmax,nx]\n"
" [-tolerance=integration-tolerance]\n"
" [-gap=<meters>] [-quiet]\n"
" Integrates particle trajectories through a symmetric or asymmetric bending magnet.\n\n"
" The data is in a SDDS-format file of (x, z, B), where x is parallel to\n"
" the line through the center of the magnet and the center of curvature,\n"
" and z is perpedicular to x.  The integration takes place in the (x, z) coordinate\n"
" system, starting at z<0.  x>0 is to the left as the particle enters the magnet,\n"
" with y pointing upwards.\n\n"
" The central bending radius can be optimized to give the desired bending angle.\n"
" -ideal            specifies parameters of an ideal hard-edge field for testing.\n"
" -3dFieldFile      specifies that field-file has (x, y, z, Bx, By, Bz).\n"
"                   The grid must be equi-spaced and sorted in (z, y, x) order.\n"
"                   E.g., sddssort <filename> -col=z,incr -col=y,incr -col=x,incr\n"
" -scan             specifies which accelerator coordinate to scan and over what range\n"
" -beamFiles        specifies elegant-style SDDS beam file, <input>, to track and\n"
"                   file, <output>, in which to place resulting beam \n"
" -vertex           Specifies (x, z) coordinates of the vertex of the magnet.\n"
" -nominalEntry     \n"
" -nominalExit      Give nominal entry and exit points for the trajectory,\n"
"                   to define the hard-edge entrance and exit points relative to which\n"
"                   accelerator coordinates and transfer matrices are defined.\n"
" -theta            bending angle in degrees (positive value) desired for this magnet.\n"
" -rigidity         Rigidity for the beam.\n"
" -fsc              fractional strength change, which may be the FSE value\n"
"                   obtained from a previous optimization with this field map\n"
" -dxDipole         Set the offset in the transverse dipole position along the \n"
"                   magnet vertical midplane.  This value can be supplied from a\n"
"                   previous run.  A positive value moves the magnet away from the\n"
"                   center of curvature.\n"
" -dzDipole         Set the offset in the dipole position perpendicular the \n"
"                   magnet vertical midplane.  A positive value moves the magnet\n"
"                   away from the incoming beam.\n"
" -yawDipole        Set the yaw angle about (0,0). A positive value rotates the\n"
"                   magnet clockwise as see from above.\n"
" -output           give file for particle trajectory output\n"
" -singleScan       input file has a single scan vs z. Untested in this version.\n"
" -arcScan          input file has a single arc scan vs s, distance along the ideal\n"
"                   beam path.  This is only correct for a sector magnet!\n"
"                   Untested in this version.\n"
" -optimize         optimize bending radius, offset, or rotation to get desired bending angle and\n"
"                   exit trajectory\n"
" -fseLimit         Set range for FSE in optimization\n"
" -dxLimit          Set range for dx of dipole in optimization\n"
" -dzLimit          Set range for dz of dipole in optimization\n"
" -yawLimit         Set range for yaw of dipole in optimization\n"
" -fieldMapOutput   outputs the field data to a file for checking\n"
" -extendData       extend data in x dimension as needed so particles see field\n"
" -zDuplicate       reflect copy of -z data into +z\n"
" -tolerance        tolerance for integration\n"
" -fieldSign        specify the sign of the bending radius (normally +)\n"
" -gap              specify the full gap of the magnet in meters.  Used to\n"
"                   compute the edge-field integral.\n"
" -quiet            Suppress informational printouts.\n\n"
"Program by Michael Borland  (Version 6, May 2016).\n";

unsigned long optimizeFlags;
#define OPTIMIZE_ON        0x0001
#define OPTIMIZE_FSE       0x0002
#define OPTIMIZE_DX        0x0004
#define OPTIMIZE_VERBOSE   0x0008
#define OPTIMIZE_QUIET     0x0010
#define OPTIMIZE_DZ        0x0020
#define OPTIMIZE_YAW       0x0040

#define TOLERANCE_FACTOR 1e-14

double particleMass = me_mks;
double particleCharge = e_mks;
double particleMassMV = me_mev;
double particleRadius = re_mks;
double particleRelSign = 1;   /* relative to electron */
long particleIsElectron = 1;

void make_fieldmap_file(char *filename, char *data_file, double Zi, double Zf, long nZ, double Xi, double Xf, long nX);
void setup_integration_output(SDDS_TABLE *SDDS_output, char *filename, char *inputfile, char *field_map_file,
                              double theta, double fse, long variable);

#include "bratSubroutines.c"

int main(int argc, char **argv)
{
  int i_arg;
  long iv;
  SCANNED_ARG *scanned;
  char *output, *input;
  char *fmap_output;
  long nx_fmap, nz_fmap;
  double zmin_fmap, zmax_fmap, xmin_fmap, xmax_fmap;
  SDDS_DATASET SDDS_output;
  double initial, final, value, delta;
  long number, variable;
  double accelCoord[6], q[10];
  long vertexGiven, entryGiven, exitGiven;
  char *inputBeamFile, *outputBeamFile;
  //double *xB, *xpB, *yB, *ypB, *sB, *deltaB;
  double xCenter, zCenter;
  char *arcSName, *arcFieldName;
  double Breference = 0;
  double idealB;

  if ((argc = scanargsg(&scanned, argc, argv))<2)
    bomb(NULL, USAGE);

  nx_fmap = nz_fmap = 0;
  integ_tol = 1e-12;
  theta = 0;
  input = fmap_output = output = NULL;
  extend_data = 0;
  single_scan = 0;
  variable = -1;
  xNomEntry = zNomEntry = zStart = zEnd = rigidity = 0;
  xVertex =  0;
  vertexGiven = entryGiven = exitGiven = 0;
  inputBeamFile = outputBeamFile = NULL;
  //xB = xpB = yB = ypB = sB = deltaB = NULL;
  arcSName = arcFieldName = NULL;
  idealMode = 0;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (scanned[i_arg].arg_type==OPTION) {
      /* process options here */
      switch (match_string(scanned[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_QUIET:
        quiet = 1;
        break;
      case SET_3DFIELDMAP:
        fieldMapDimension = 3;
        break;
      case SET_IDEAL:
        if (scanned[i_arg].n_items!=4 ||
            sscanf(scanned[i_arg].list[1], "%lf", &idealB)!=1 ||
            sscanf(scanned[i_arg].list[2], "%lf", &idealChord)!=1 ||
            sscanf(scanned[i_arg].list[3], "%lf", &idealEdgeAngle)!=1)
          bomb("invalid -ideal syntax/values", NULL);
        idealMode = 1;
        break;
      case SET_BEAMFILES:
        if (scanned[i_arg].n_items!=3)
          bomb("invalid -beamFiles syntax", NULL);
        inputBeamFile = scanned[i_arg].list[1];
        outputBeamFile = scanned[i_arg].list[2];
        break;
      case SET_TOLERANCE:
        if (scanned[i_arg].n_items!=2 ||
            sscanf(scanned[i_arg].list[1], "%lf", &integ_tol)!=1 ||
            integ_tol<=TOLERANCE_FACTOR)
          bomb("invalid -tolerance syntax/value", NULL);
        break;
      case SET_THETA:
        if (scanned[i_arg].n_items!=2 ||
            sscanf(scanned[i_arg].list[1], "%lf", &theta)!=1 ||
            theta<=0)
          bomb("invalid -theta syntax/value", NULL);
        break;
      case SET_RIGIDITY:
        if (scanned[i_arg].n_items!=2 ||
            sscanf(scanned[i_arg].list[1], "%lf", &rigidity)!=1 ||
            rigidity<=0)
          bomb("invalid -rigidity syntax/value", NULL);
        break;
      case SET_FSC:
        if (scanned[i_arg].n_items!=2 ||
            sscanf(scanned[i_arg].list[1], "%lf", &fse)!=1)
          bomb("invalid fsc syntax/value", NULL);
        break;
      case SET_DX_OFFSET:
        if (scanned[i_arg].n_items!=2 ||
            sscanf(scanned[i_arg].list[1], "%lf", &dXOffset)!=1)
          bomb("invalid dxDipole syntax/value", NULL);
        break;
      case SET_DZ_OFFSET:
        if (scanned[i_arg].n_items!=2 ||
            sscanf(scanned[i_arg].list[1], "%lf", &dZOffset)!=1)
          bomb("invalid dzDipole syntax/value", NULL);
        break;
      case SET_YAW:
        if (scanned[i_arg].n_items!=2 ||
            sscanf(scanned[i_arg].list[1], "%lf", &magnetYaw)!=1)
          bomb("invalid -yaw syntax/value", NULL);
        break;
      case SET_OUTPUT:
        if (scanned[i_arg].n_items!=2)
          bomb("invalid -output syntax", NULL);
        output = scanned[i_arg].list[1];
        break;
      case SET_SCAN:
        if (scanned[i_arg].n_items!=5 ||
            (variable=match_string(scanned[i_arg].list[1], variable_name, 5, EXACT_MATCH))<0 ||
            sscanf(scanned[i_arg].list[2], "%lf", &initial)!=1 ||
            sscanf(scanned[i_arg].list[3], "%lf", &final)!=1 ||
            sscanf(scanned[i_arg].list[4], "%ld", &number)!=1 ||
            initial>=final || number<1)
          bomb("invalid -scan syntax", NULL);
        break;
      case SET_FIELDMAP_OUTPUT:
        if (scanned[i_arg].n_items!=8)
          bomb("invalid -fieldmap_output syntax", NULL);
        if (sscanf(scanned[i_arg].list[2], "%lf", &zmin_fmap)!=1 ||
            sscanf(scanned[i_arg].list[3], "%lf", &zmax_fmap)!=1 ||
            sscanf(scanned[i_arg].list[4], "%ld", &nz_fmap)!=1 ||
            sscanf(scanned[i_arg].list[5], "%lf", &xmin_fmap)!=1 ||
            sscanf(scanned[i_arg].list[6], "%lf", &xmax_fmap)!=1 ||
            sscanf(scanned[i_arg].list[7], "%ld", &nx_fmap)!=1 ||
            zmin_fmap>=zmax_fmap || xmin_fmap>=xmax_fmap || 
            nz_fmap<2 || nx_fmap<2)
            bomb("invalid -fieldmap_output syntax", NULL);
        fmap_output = scanned[i_arg].list[1];
        break;
      case SET_OPTIMIZE:
        optimizeFlags = OPTIMIZE_ON;
        scanned[i_arg].n_items -= 1;
        if (!scanItemList(&optimizeFlags, scanned[i_arg].list+1, &scanned[i_arg].n_items, 0, 
                          "verbose", -1, NULL, 0, OPTIMIZE_VERBOSE,
                          "fse",  -1, NULL, 0, OPTIMIZE_FSE,
                          "dx",   -1, NULL, 0, OPTIMIZE_DX,
                          "dz",   -1, NULL, 0, OPTIMIZE_DZ,
                          "yaw",  -1, NULL, 0, OPTIMIZE_YAW,
                          NULL))
          bomb("invalid -optimize syntax", NULL);
        if (optimizeFlags&OPTIMIZE_VERBOSE)
          verbose_optimize = 1;
        break;
      case SET_EXTEND_DATA:
        extend_data = 1;
        if (scanned[i_arg].n_items!=2 ||
            sscanf(scanned[i_arg].list[1], "%lf", &extend_edge_angle)!=1)
          bomb("invalid -extend syntax", USAGE);
        break;
      case SET_SINGLE_SCAN:
        single_scan = 1;
        break;
      case SET_ARC_SCAN:
        if (scanned[i_arg].n_items!=4 ||
            !strlen(arcSName = scanned[i_arg].list[1]) ||
            !strlen(arcFieldName = scanned[i_arg].list[2]) ||
            sscanf(scanned[i_arg].list[3], "%lf", &rhoArc)!=1 ||
            rhoArc<=0)
          bomb("invalid -arcScan syntax/values", USAGE);
        arc_scan = 1;
        break;
      case SET_VERTEX:
        if (scanned[i_arg].n_items!=3 ||
            sscanf(scanned[i_arg].list[1], "%lf", &xVertex)!=1 ||
            sscanf(scanned[i_arg].list[2], "%lf", &zVertex)!=1)
          bomb("invalid -vertex syntax/values", USAGE);
        vertexGiven = 1;
        break;
      case SET_ENTRY:
        if (scanned[i_arg].n_items!=3 ||
            sscanf(scanned[i_arg].list[1], "%lf", &xNomEntry)!=1 ||
            sscanf(scanned[i_arg].list[2], "%lf", &zNomEntry)!=1)
          bomb("invalid -entry syntax/values", USAGE);
        entryGiven = 1;
        break;
      case SET_EXIT:
        if (scanned[i_arg].n_items!=3 ||
            sscanf(scanned[i_arg].list[1], "%lf", &xNomExit)!=1 ||
            sscanf(scanned[i_arg].list[2], "%lf", &zNomExit)!=1)
          bomb("invalid -exit syntax/values", USAGE);
        exitGiven = 1;
        break;
      case SET_FIELDSIGN:
        if (scanned[i_arg].n_items!=2)
          bomb("invalid -fieldSign syntax/value", USAGE);
        switch (scanned[i_arg].list[1][0]) {
        case '-':
          fieldSign = -1;
          break;
        case '+':
          fieldSign = 1;
          break;
        default:
          bomb("invalid -fieldSign syntax/value", USAGE);
          break;
        }
        break;
      case SET_GAP:
        if (scanned[i_arg].n_items!=2 ||
            sscanf(scanned[i_arg].list[1], "%lf", &gap)!=1 ||
            gap<=0)
          bomb("invalid -gap syntax/values", USAGE);
        break;
      case SET_FSE_LIMIT:
        if (scanned[i_arg].n_items!=3 ||
            sscanf(scanned[i_arg].list[1], "%lf", &fseLimit[0])!=1 ||
            sscanf(scanned[i_arg].list[2], "%lf", &fseLimit[1])!=1 ||
            fseLimit[0]>=fseLimit[1])
          bomb("invalid -fseLimit syntax/values", USAGE);
        break;
      case SET_DX_LIMIT:
        if (scanned[i_arg].n_items!=3 ||
            sscanf(scanned[i_arg].list[1], "%lf", &dxLimit[0])!=1 ||
            sscanf(scanned[i_arg].list[2], "%lf", &dxLimit[1])!=1 ||
            dxLimit[0]>=dxLimit[1])
          bomb("invalid -dxLimit syntax/values", USAGE);
        break;
      case SET_DZ_LIMIT:
        if (scanned[i_arg].n_items!=3 ||
            sscanf(scanned[i_arg].list[1], "%lf", &dzLimit[0])!=1 ||
            sscanf(scanned[i_arg].list[2], "%lf", &dzLimit[1])!=1 ||
            dzLimit[0]>=dzLimit[1])
          bomb("invalid -dzLimit syntax/values", USAGE);
        break;
      case SET_YAW_LIMIT:
        if (scanned[i_arg].n_items!=3 ||
            sscanf(scanned[i_arg].list[1], "%lf", &yawLimit[0])!=1 ||
            sscanf(scanned[i_arg].list[2], "%lf", &yawLimit[1])!=1 ||
            yawLimit[0]>=yawLimit[1])
          bomb("invalid -yawLimit syntax/values", USAGE);
        break;
      case SET_Z_DUPLICATE:
        zDuplicate = 1;
        break;
      default:
        fprintf(stderr, "option %s not recognized\n", scanned[i_arg].list[0]);
        bomb("unknown option given", USAGE);
        break;
      }
    }
    else {
      /* filename argument */
      if (!input)
        input = scanned[i_arg].list[0];
      else
        bomb("too many filenames", NULL);
    }
  }
  if (!idealMode && !input)
    bomb("no input file given", USAGE);

  zero_tol = integ_tol*10;
  delta = (final-initial)/(number-1);

  if (theta<=0)
    bomb("you must specify theta > 0", USAGE);
  if (!entryGiven)
    bomb("you must give -entry", USAGE);
  if (!vertexGiven)
    bomb("you must give -vertex", USAGE);
  if (!exitGiven)
    bomb("you must give -exit", USAGE);
  if (rigidity<=0) 
    bomb("you must give -rigidity", USAGE);
    
  theta *= PI/180;        /* degrees to radians */

  zCenter = (zNomEntry+zNomExit)/2;
  xCenter = (xNomEntry+xNomExit)/2;
  central_length = 100*(zNomExit-zNomEntry);
  if (!idealMode) {
    if (arc_scan) {
      Breference = BRAT_setup_arc_field_data(input, arcSName, arcFieldName, xCenter);
    } else {
      Breference = BRAT_setup_field_data(input, xCenter, zCenter);
      zStart = zi-dz;
      if (zStart>zNomEntry) {
        fprintf(stderr, "zStart = %e, zNomEntry = %e\n",
                zStart, zNomEntry);
        exit(1);
      }
    }
  } else {
    Breference = idealB;
    zStart = -idealChord;
  }
  rhoMax = rigidity/Breference;
  
  /* set variables for saving the path */
  max_store = 0;
  X_stored = Z_stored = Y_stored = s_stored = NULL;
  wX_stored = wZ_stored = wY_stored = NULL;
  aX_stored = aZ_stored = aY_stored = NULL;
  FX_stored = FZ_stored = FY_stored = Fint_stored = NULL;
  if (arc_scan)
    z_outer = 3*(zNomExit-zNomEntry);
  else if (idealMode)
    z_outer = 3*idealChord;
  else
    z_outer = MAX(fabs(zi), fabs(zf));
  
  if (optimizeFlags)
    BRAT_optimize_magnet(optimizeFlags);

  accelCoord[0] = accelCoord[1] = accelCoord[2] = accelCoord[3] = accelCoord[4] = accelCoord[5] = 0;
  BRAT_lorentz_integration(accelCoord, q, 1, &Breference);
  rhoMax = rigidity/Breference;
  if (!quiet) fprintf(stderr, "Breference = %le\n", Breference);

  if (fmap_output) {
    if (single_scan || arc_scan)
      puts("warning: contour output not supported in -singleScan or -arcScan mode");
    else {
      if (!nx_fmap || !nz_fmap) {
        nx_fmap = 2*nx;
        nz_fmap = 2*nz;
      }
      if (!quiet) fprintf(stderr, "Making %ld by %ld field map\n", nx_fmap, nz_fmap);
      make_fieldmap_file(fmap_output, input, zmin_fmap, zmax_fmap, nz_fmap, xmin_fmap, xmax_fmap, nx_fmap);
    }
 }

  if (output) {
    if (inputBeamFile)
      bomb("Use of -beamFiles option is incompatible with -output option", NULL);
    setup_integration_output(&SDDS_output, output, input, fmap_output, theta*180/PI, fse,
                             variable);
  }
  
  if (!inputBeamFile) {
    if (variable==-1) {
      number = 1;
      initial = value = delta = 0;
    }
    
    /* do integration for series of particles */
    for (iv=0, value=initial; iv<number; iv++, value+=delta) {
      accelCoord[0] = accelCoord[1] = accelCoord[2] = accelCoord[3] = accelCoord[4] = accelCoord[5] = 0;
      switch (variable) {
      case SCAN_X:
        if (!quiet) printf("working on xa=%e m:  ", value);
        accelCoord[0] = value;
        break;
      case SCAN_XP:
        if (!quiet) printf("working on xpa=%e m:  ", value);
        accelCoord[1] = value;
        break;
      case SCAN_Y:
        if (!quiet) printf("working on ya=%e m:  ", value);
        accelCoord[2] = value;
        break;
      case SCAN_YP:
        if (!quiet) printf("working on ypa=%e m:  ", value);
        accelCoord[3] = value;
        break;
      case SCAN_DELTA:
        if (!quiet) printf("working on delta=%e m:  ", value);
        accelCoord[5] = value;
        break;
      default:
        break;
      }
      BRAT_lorentz_integration(accelCoord, q, 1, NULL);
      if (!quiet) printf("\n");
      fflush(stdout);
      if (n_stored>0 && Z_stored[n_stored-1]<0)
        n_stored--;
      if (output) {
        if (n_stored<1) {
          fprintf(stderr, "z_outer = %le\n", z_outer);
          bomb("no data stored for particle path output", NULL);
        }
        if (!SDDS_StartTable(&SDDS_output, n_stored) ||
            !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, Z_stored, n_stored, "Z") ||
            !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, X_stored, n_stored, "X") ||
            !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, Y_stored, n_stored, "Y") ||
            !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, wZ_stored, n_stored, "wZ") ||
            !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, wX_stored, n_stored, "wX") ||
            !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, wY_stored, n_stored, "wY") ||
            !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, aZ_stored, n_stored, "aZ") ||
            !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, aX_stored, n_stored, "aX") ||
            !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, aY_stored, n_stored, "aY") ||
            !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, s_stored, n_stored, "s") ||
            !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, FZ_stored, n_stored, "BZ") ||
            !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, FX_stored, n_stored, "BX") ||
            !SDDS_SetColumn(&SDDS_output, SDDS_SET_BY_NAME, FY_stored, n_stored, "BY") ||
            !SDDS_SetParameters(&SDDS_output, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                                "xf", accelCoord[0], "xpf", accelCoord[1],
                                "yf", accelCoord[2], "ypf", accelCoord[3],
                                "sf", accelCoord[4], "FSE", fse,
                                "dXDipole", dXOffset, "dZDipole", dZOffset,
                                "YawDipole", magnetYaw,
                                "BReference", Breference,
                                (variable==-1?NULL:initial_variable_name[variable]), value, 
                                NULL) ||
            !SDDS_WriteTable(&SDDS_output)) {
          SDDS_SetError("Problem writing SDDS table for integration output");
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
      }
    }
  } else {
    double *x=NULL, *xp=NULL, *y=NULL, *yp=NULL, *t=NULL, *p=NULL;
    double pCentral, beta;
    SDDS_DATASET inputBeam, outputBeam;
    if (!SDDS_InitializeInput(&inputBeam, inputBeamFile) ||
        !SDDS_ReadPage(&inputBeam)) {
      SDDS_SetError("Unable to read beam file");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (!(number=SDDS_RowCount(&inputBeam)))
      bomb("no data in beam input file", NULL);
    if (SDDS_CheckColumn(&inputBeam, "x", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY ||
        SDDS_CheckColumn(&inputBeam, "y", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY ||
        SDDS_CheckColumn(&inputBeam, "xp", NULL, SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY ||
        SDDS_CheckColumn(&inputBeam, "yp", NULL, SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY ||
        SDDS_CheckColumn(&inputBeam, "t", "s", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OKAY)
      bomb("missing data or wrong units in input beam file", NULL);
    if (SDDS_CheckColumn(&inputBeam, "p", "m$be$nc", SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OKAY &&
        SDDS_CheckColumn(&inputBeam, "p", "", SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OKAY) {
      bomb("column p is missing or has wrong units", NULL);
    }
    if (SDDS_CheckParameter(&inputBeam, "pCentral", "m$be$nc", SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OKAY &&
        SDDS_CheckParameter(&inputBeam, "pCentral", "", SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OKAY) {
      bomb("Parameter pCentral is missing or has wrong units", NULL);
    }
    if (!(x = SDDS_GetColumnInDoubles(&inputBeam, "x")) ||
        !(xp = SDDS_GetColumnInDoubles(&inputBeam, "xp")) ||
        !(y = SDDS_GetColumnInDoubles(&inputBeam, "y")) ||
        !(yp = SDDS_GetColumnInDoubles(&inputBeam, "yp")) ||
        !(t = SDDS_GetColumnInDoubles(&inputBeam, "t")) ||
        !(p = SDDS_GetColumnInDoubles(&inputBeam, "p")) ||
        !SDDS_GetParameter(&inputBeam, "pCentral", &pCentral))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (!SDDS_InitializeCopy(&outputBeam, &inputBeam, outputBeamFile, "w") ||
        !SDDS_SaveLayout(&outputBeam) || !SDDS_WriteLayout(&outputBeam) || 
        !SDDS_StartPage(&outputBeam, number) ||
        !SDDS_CopyParameters(&outputBeam, &inputBeam))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    for (iv=0; iv<number; iv++) {
      accelCoord[0] = x[iv];
      accelCoord[1] = xp[iv];
      accelCoord[2] = y[iv];
      accelCoord[3] = yp[iv];
      accelCoord[4] = 0;
      accelCoord[5] = (p[iv]-pCentral)/pCentral;
      BRAT_lorentz_integration(accelCoord, q, 0, NULL);
      x[iv] = accelCoord[0];
      xp[iv] = accelCoord[1];
      y[iv] = accelCoord[2];
      yp[iv] = accelCoord[3];
      p[iv] = (1+accelCoord[5])*pCentral;
      beta = p[iv]/sqrt(p[iv]*p[iv]+1);
      t[iv] += accelCoord[4]/(beta*c_mks);
    }
    if (!SDDS_SetColumnFromDoubles(&outputBeam, SDDS_SET_BY_NAME, x, number, "x") ||
        !SDDS_SetColumnFromDoubles(&outputBeam, SDDS_SET_BY_NAME, xp, number, "xp") ||
        !SDDS_SetColumnFromDoubles(&outputBeam, SDDS_SET_BY_NAME, y, number, "y") ||
        !SDDS_SetColumnFromDoubles(&outputBeam, SDDS_SET_BY_NAME, yp, number, "yp") ||
        !SDDS_SetColumnFromDoubles(&outputBeam, SDDS_SET_BY_NAME, t, number, "t") ||
        !SDDS_SetColumnFromDoubles(&outputBeam, SDDS_SET_BY_NAME, p, number, "p") ||
        !SDDS_WritePage(&outputBeam) || !SDDS_Terminate(&outputBeam) ||
        !SDDS_Terminate(&inputBeam) )
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  
  return(0);
}

void bombElegant(const char *error, const char *usage)
{
  if (error)
    fprintf(stderr, "error: %s\n", error);
  if (usage)
    fprintf(stderr, "usage: %s\n", usage);
  exit(1);
}

void make_fieldmap_file(char *filename, char *data_file, double Zi, double Zf, long nZ, double Xi, double Xf, long nX)
{
  char s[200];
  double dX, dZ, X, Z, B[3], Q[3];
  long iZ, iX, row;
  SDDS_TABLE SDDS_table;

  if (Xi>Xf)
    SWAP_DOUBLE(Xi, Xf);
  if (Zi>Zf)
    SWAP_DOUBLE(Zi, Zf);
  dX = (Xf - Xi)/(nX-1);
  dZ = (Zf - Zi)/(nZ-1);

  sprintf(s, "Contours of constant field for file %s", data_file);
  if (!SDDS_InitializeOutput(&SDDS_table, SDDS_BINARY, 1, s, "fieldmap", filename) ||
      SDDS_DefineColumn(&SDDS_table, "By", "B$by$n", "T", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDS_table, "Bx", "B$bx$n", "T", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDS_table, "Bz", "B$bz$n", "T", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDS_table, "x", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDS_table, "z", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      !SDDS_WriteLayout(&SDDS_table) || !SDDS_StartTable(&SDDS_table, nX*nZ)) {
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
     
  Q[2] = 0;
  for (iZ=row=0; iZ<nZ; iZ++) {
    particle_inside = 0;
    if ((Z = Zi + iZ*dZ)>Zf)
      Z = Zf;
    for (iX=0; iX<nX; iX++) {
      if ((X = Xi + iX*dX)>Xf)
        X = Xf;
      Q[0] = Z;
      Q[1] = X;
      BRAT_B_field(B, Q);
      if (!SDDS_SetRowValues(&SDDS_table, SDDS_PASS_BY_VALUE|SDDS_SET_BY_INDEX, row++, 0, B[2],
                             1, B[1], 2, B[0], 3, X, 4, Z, -1)) {
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
    }
  }

  if (!SDDS_WriteTable(&SDDS_table) || !SDDS_Terminate(&SDDS_table))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
}

void setup_integration_output(SDDS_TABLE *SDDS_output, char *filename, char *inputfile, char *field_map_file,
                              double theta, double fse, long variable)
{
  char s[200];
  sprintf(s, "brat output for file %s with $gQ$r=%g$ao$n, fse=%.4f",
          inputfile, theta, fse);
  if (!SDDS_InitializeOutput(SDDS_output, SDDS_BINARY, 1,s, "brat output", filename) ||
      SDDS_DefineParameter(SDDS_output, "theta", "$gQ$r", "degrees", "Total bending angle of magnet", NULL, 
                           SDDS_DOUBLE, (sprintf(s, "%.15e", theta), s))<0 ||
      SDDS_DefineParameter(SDDS_output, "rigidity", "H", "T*m", "Rigidity of electron beam", NULL, 
                           SDDS_DOUBLE, (sprintf(s, "%.15e", rigidity), s))<0 ||
      SDDS_DefineParameter(SDDS_output, "FSE", "$gr$r", "", "Optimal FSE", NULL, 
                           SDDS_DOUBLE, (sprintf(s, "%.15e", fse), s))<0 ||
      SDDS_DefineParameter(SDDS_output, "xf", "x$bf$n", "m", "Final x", NULL, SDDS_DOUBLE, NULL)<0 ||
      SDDS_DefineParameter(SDDS_output, "xpf", "x'$bf$n", NULL, "Final x'", NULL, SDDS_DOUBLE, NULL)<0 ||
      SDDS_DefineParameter(SDDS_output, "yf", "y$bf$n", "m", "Final y", NULL, SDDS_DOUBLE, NULL)<0 ||
      SDDS_DefineParameter(SDDS_output, "ypf", "y'$bf$n", NULL, "Final y'", NULL, SDDS_DOUBLE, NULL)<0 ||
      SDDS_DefineParameter(SDDS_output, "sf", "s$bf$n", NULL, "Final path length", NULL, SDDS_DOUBLE, NULL)<0 ||
      SDDS_DefineParameter(SDDS_output, "BReference", "B$bref$n", "T", "Reference field value", NULL, SDDS_DOUBLE, NULL)<0 ||
      SDDS_DefineParameter(SDDS_output, "dXDipole", "", "m", "Displacement of dipole in X",
                           NULL, SDDS_DOUBLE, NULL)<0 ||
      SDDS_DefineParameter(SDDS_output, "dZDipole", "", "m", "Displacement of dipole in Z",
                           NULL, SDDS_DOUBLE, NULL)<0 ||
      SDDS_DefineParameter(SDDS_output, "YawDipole", "", "rad", "Yaw of dipole",
                           NULL, SDDS_DOUBLE, NULL)<0 ||
      SDDS_DefineColumn(SDDS_output, "s", "s", "m", "Path length", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_output, "Z", "Z", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_output, "X", "X", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_output, "Y", "Y", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_output, "wZ", "wZ", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_output, "wX", "wX", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_output, "wY", "wY", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_output, "aZ", "aZ", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_output, "aX", "aX", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_output, "aY", "aY", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_output, "BZ", "BZ", "T", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_output, "BX", "BX", "T", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDS_output, "BY", "BY", "T", NULL, NULL, SDDS_DOUBLE, 0)<0) {
      SDDS_SetError("Problem setting up SDDS table for integration output");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    
  }
  if (variable!=-1 &&
      SDDS_DefineParameter(SDDS_output, initial_variable_name[variable], 
                           initial_variable_symbol[variable], variable_units[variable],
                           variable_description[variable], NULL, SDDS_DOUBLE, NULL)<0) {
    SDDS_SetError("Problem setting up SDDS table for integration output");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (!SDDS_WriteLayout(SDDS_output)) {
    SDDS_SetError("Problem setting up SDDS table for integration output");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
}

