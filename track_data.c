/************************************************************************* \
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: track_data4.c
 * contents: data arrays for elegant 
 *
 * Michael Borland, 1989-1993
 */
#include "mdb.h"
#include "track.h"

/* flag word for trace mode */
long trace_mode = 0;

/* default values for electron tracking */
double particleMass = me_mks;
double particleCharge = e_mks;
double particleMassMV = me_mev;
double particleRadius = re_mks;
double particleRelSign = 1;   /* relative to electron */
long particleIsElectron = 1;

/* A hash table for loading parameters effectively */
htab *load_hash;

/* various user-controlled global flags and settings (global_settings namelist) */
long inhibitFileSync = 0;
long echoNamelists = 1;
long mpiRandomizationMode = 3;
long exactNormalizedEmittance = 0;
long shareTrackingBasedMatrices = 1;
long parallelTrackingBasedMatrices = 1;
double trackingMatrixStepFactor = 1;
long trackingMatrixPoints = 9;
double trackingMatrixStepSize[6] = {5e-5, 5e-5, 5e-5, 5e-5, 5e-5, 5e-5};
#if USE_MPI
short mpiAbort = 0;
#endif

long trajectoryTracking = 0;

long particleID = 1;

/* number of sigmas for gaussian random numbers in radiation emission simulation in CSBEND, KQUAD, etc. */
double srGaussianLimit = 3.0;

char *entity_name[N_TYPES] = {
    "LINE", "QUAD", "SBEN", "RBEN", "DRIF", "SEXT", "OCTU", "MULT", "SOLE", 
    "HKICK", "VKICK", "RFCA", "ELSE", "HMON", "VMON", "MONI", "RCOL", "ECOL", 
    "MARK", "MATR", "ALPH", "RFDF", "RFTMEZ0", "RMDF", "TMCF", "CEPL", "WATCH",
    "TWPL", "MALIGN", "TWLA", "PEPPOT", "ENERGY", "MAXAMP", "ROTATE",
    "TRCOUNT", "RECIRC", "QUFRINGE", "SCRAPER", "CENTER", "BUMPER",
    "KSEXT", "KSBEND", "KQUAD", "MAGNIFY", "SAMPLE", "KICKER", "SCATTER",
    "NIBEND", "KPOLY", "NISEPT", "RAMPRF", "RAMPP", "STRAY", "CSBEND",
    "TWMTA", "MATTER", "RFMODE", "TRFMODE", "ZLONGIT", "SREFFECTS",
    "MODRF", "BMAPXY", "ZTRANSVERSE", "IBSCATTER", "FMULT",
    "WAKE", "TRWAKE", "TUBEND", "CHARGE", "PFILTER", "HISTOGRAM",
    "CSRCSBEND", "CSRDRIFT", "RFCW", "REMCOR", "MAPSOLENOID",
    "REFLECT", "CLEAN", "TWISS", "WIGGLER", "SCRIPT", "FLOOR",
    "LTHINLENS", "LMIRROR", "EMATRIX", "FRFMODE", "FTRFMODE",
    "TFBPICKUP", "TFBDRIVER", "LSCDRIFT", "DSCATTER", "LSRMDLTR",
    "TAYLORSERIES", "RFTM110", "CWIGGLER", "EDRIFT", "SCMULT", "ILMATRIX",
    "TSCATTER", "KQUSE", "UKICKMAP", "MBUMPER", "EMITTANCE", "MHISTOGRAM", 
    "FTABLE", "KOCT", "RIMULT", "GFWIGGLER", "MRFDF", "CORGPIPE", "LRWAKE",
    "EHKICK", "EVKICK", "EKICKER", "BMXYZ", "BRAT", "BGGEXP", "BRANCH",
    "IONEFFECTS", "SLICE", "SPEEDBUMP", "CCBEND", "HKPOLY", "BOFFAXE",
    "APCONTOUR", "TAPERAPC", "TAPERAPE", "TAPERAPR", "SHRFDF",
    };

char *madcom_name[N_MADCOMS] = {
    "", "USE", "TITLE", "RETURN", "RENAME"
        } ;

char *entity_text[N_TYPES] = {
    NULL, 
    "A quadrupole implemented as a matrix, up to 3rd order. Use KQUAD for symplectic tracking.",
    "A sector dipole implemented as a matrix, up to 2nd order. Use CSBEND for symplectic tracking.",
    "A rectangular dipole, implemented as a SBEND with edge angles, up to 2nd order. Use CSBEND for symplectic tracking.",
    "A drift space implemented as a matrix, up to 2nd order. Use EDRIFT for symplectic tracking.",
    "A sextupole implemented as a matrix, up to 3rd order. Use KSEXT for symplectic tracking.",
    "An octupole implemented as a third-order matrix. Use KOCT for symplectic tracking.",
    "A canonical kick multipole.",
    "A solenoid implemented as a matrix, up to 2nd order.",
    "A horizontal steering dipole implemented as a matrix, up to 2nd order. Use EHKICK for symplectic tracking.",
    "A vertical steering dipole implemented as a matrix, up to 2nd order. Use EVKICK for symplectic tracking.",
    "A first-order matrix RF cavity with exact phase dependence.",
    "Not implemented.",
    "A horizontal position monitor, accepting a rpn equation for the readout as a function of the actual position (x).",
    "A vertical position monitor, accepting a rpn equation for the readout as a function of the actual position (y).",
    "A two-plane position monitor, accepting two rpn equations for the readouts as a function of the actual positions (x and y).",
    "A rectangular collimator.",
    "An elliptical collimator.",
    "A marker, equivalent to a zero-length drift space.",
    "Explicit matrix input from a text file, in the format written by the print_matrix command.",
    "An alpha magnet implemented as a matrix, up to 3rd order.",
    "A simple traveling or standing wave deflecting RF cavity.",
    "A TM-mode RF cavity specified by the on-axis Ez field.",
    "A linearly-ramped electric field deflector, using an approximate analytical solution FOR LOW ENERGY PARTICLES.",
    "A numerically-integrated accelerating TM RF cavity with spatially-constant fields.",
    "A numerically-integrated linearly-ramped electric field deflector.",
    "A beam property/motion monitor--allowed modes are centroid, parameter, coordinate, and fft.",
    "A numerically-integrated traveling-wave stripline deflector.",
    "A misalignment of the beam, implemented as a zero-order matrix.",
    "A numerically-integrated first-space-harmonic traveling-wave linear accelerator.",
    "A pepper-pot plate.",
    "An element that matches the central momentum to the beam momentum, or changes the central momentum or energy to a specified value.",
    "A collimating element that sets the maximum transmitted particle amplitudes for all following elements, until the next MAXAMP.",
    "An element that rotates the beam about the longitudinal axis.",
    "An element that defines the point from which transmission calculations are made.",
    "An element that defines the point to which particles recirculate in multi-pass tracking",
    "An element consisting of a linearly increasing or decreasing quadrupole field.",
    "A collimating element that sticks into the beam from one side only.  The directions 0, 1, 2, and 3 are from +x, +y, -x, and -y, respectively.",
    "An element that centers the beam transversely on the ideal trajectory.",
    "A time-dependent kicker magnet with optional spatial dependence of the kick and no fringe effects. The waveform is in SDDS format, with time in seconds and amplitude normalized to 1. The optional spatial dependence is also specified as an SDDS file.",
    "A canonical kick sextupole, which differs from the MULT element with ORDER=2 in that it can be used for chromaticity correction.",
    "A kick bending magnet which is NOT canonical, but is better than a 2nd order matrix implementation.  Use CSBEND instead.",
    "A canonical kick quadrupole.",
    "An element that allows multiplication of phase-space coordinates of all particles by constants.",
    "An element that reduces the number of particles in the beam by interval-based or random sampling.",
    "A combined horizontal-vertical steering magnet implemented as a matrix, up to 2nd order. For time-dependent kickers, see BUMPER.",
    "A scattering element to add gaussian random numbers to particle coordinates.",
    "A numerically-integrated dipole magnet with various extended-fringe-field models.",
    "A thin kick element with polynomial dependence on the coordinates in one plane.",
    "A numerically-integrated dipole magnet with a Cartesian gradient.",
    "A voltage-, phase-, and/or frequency-ramped RF cavity, implemented like RFCA.",
    "A momentum-ramping element that changes the central momentum according to an SDDS-format file of the momentum factor vs time in seconds.",
    "A stray field element with local and global components.  Global components are defined relative to the initial beamline direction.",
    "A canonical kick sector dipole magnet.",
    "A numerically-integrated traveling-wave muffin-tin accelerator.",
    "A Coulomb-scattering and energy-absorbing element simulating material in the beam path.",
    "A simulation of a beam-driven TM monopole mode of an RF cavity.",
    "A simulation of a beam-driven TM dipole mode of an RF cavity.",
    "A simulation of a single-pass broad-band or functionally specified longitudinal impedance.",
    "Lumped simulation of synchrotron radiation effects (damping and quantum excitation) for rings.",
    "A first-order matrix RF cavity with exact phase dependence, plus optional amplitude and phase modulation.",
    "A map of Bx and By vs x and y.",
    "A simulation of a single-pass broad-band or functionally-specified transverse impedance.",
    "A simulation of intra-beam scattering.",
    "Multipole kick element with coefficient input from an SDDS file.",
    "Longitudinal wake specified as a function of time lag behind the particle.",
    "Transverse wake specified as a function of time lag behind the particle.",
    "A special rectangular bend element for top-up backtracking.",
    "An element to establish the total charge of a beam.  Active on first pass only.  If given, overrides all charge specifications on other elements.",
    "An element for energy and momentum filtration.",
    "Request for histograms of particle coordinates to be output to SDDS file.",
    "Like CSBEND, but incorporates a simulation of Coherent Synchrotron radiation.",
    "A follow-on element for CSRCSBEND that applies the CSR wake over a drift.",
    "A combination of RFCA, WAKE, TRWAKE, and LSCDRIFT.",
    "An element to remove correlations from the tracked beam to simulate certain types of correction.",
    "A numerically-integrated solenoid specified as a map of (Bz, Br) vs (z, r).",
    "Reflects the beam back on itself, which is useful for multiple beamline matching.",
    "Cleans the beam by removing outlier particles.",
    "Sets Twiss parameter values.",
    "A wiggler or undulator for damping or excitation of the beam.",
    "An element that allows transforming the beam using an external script.",
    "Sets floor coordinates",
    "A thin lens for light optics",
    "A mirror for light optics",
    "Explicit matrix input with data in the element definition, rather than in a file.",
    "One or more beam-driven TM monopole modes of an RF cavity, with data from a file.",
    "One or more beam-driven TM dipole modes of an RF cavity, with data from a file.",
    "Pickup for a turn-by-turn feedback loop",
    "Driver for a turn-by-turn feedback loop",
    "Longitudinal space charge impedance",
    "A scattering element to add random changes to particle coordinates according to a user-supplied distribution function",
    "A non-symplectic numerically integrated planar undulator including optional co-propagating laser beam for laser modulation of the electron beam.",
    "Tracks through a Taylor series map specified by a file containing coefficients.",
    "Tracks through a TM110-mode (deflecting) rf cavity with all magnetic and electric field components.  NOT RECOMMENDED---See below.",
    "Tracks through a wiggler using canonical integration routines of Y. Wu (Duke University).",
    "Tracks through a drift with no approximations (Exact DRIFT).",
    "Tracks through a zero length multipole to simulate space charge effects",  
    "An Individualized Linear Matrix for each particle for fast symplectic tracking with chromatic and amplitude-dependent effects",
    "An element to simulate Touschek scattering.",
    "A canonical kick element combining quadrupole and sextupole fields.",
    "An undulator kick map (e.g., using data from RADIA).",
    "A time-dependent multipole kicker magnet. The waveform is in SDDS format, with time in seconds and amplitude normalized to 1.",
    "Applies a linear transformation to the beam to force the emittance to given values.",
    "Request for multiple dimensions (1, 2, 4 or 6) histogram output of particle coordinates.",
    "Tracks through a magnetic field which is expressed by a SDDS table.",
    "A canonical kick octupole.",
    "Multiplies radiation integrals by a given factor.  Use to compute emittance for collection of various types of cells.",
    "Tracks through a wiggler using generate function method of J. Bahrdt and G. Wuestefeld (BESSY, Berlin, Germany).",
    "Zero-length Multipole RF DeFlector from dipole to decapole",
    "A corrugated round pipe, commonly used as a dechirper in linacs.",
    "Long-range (inter-bunch and inter-turn) longitudinal and transverse wake",
    "A horizontal steering dipole implemented using an exact hard-edge model",
    "A vertical steering dipole implemented using an exact hard-edge model",
    "A combined horizontal/vertical steering dipole implemented using an exact hard-edge model",
    "A map of (Bx, By, Bz) vs (x, y, z), for straight elements only",
    "Bending magnet RAy Tracing using (Bx, By, Bz) vs (x, y, z).",
    "A magnetic field element using generalized gradient expansion.",
    "Conditional branch instruction to jump to another part of the beamline",
    "Simulates ionization of residual gas and interaction with the beam.",
    "Performs slice-by-slice analysis of the beam for output to a file.",
    "Simulates a semi-circular protuberance from one or both walls of the chamber.",
    "A canonically-integrated straight dipole magnet, assumed to have multipoles defined in Cartesian coordinates.",
    "Applies kick according to a Hamiltonian that's a polynomial function of x and y together with a generalized drift also given as a polynomial of qx and qy",
    "A straight magnetic field element using off-axis expansion from an on-axis derivative.",
    "An aperture (or its inverse) defined by (x, y) points in an SDDS file.",
    "A tapered aperture that is a section of a circular cylinder.",
    "A tapered elliptical aperture.",
    "A tapered rectangular aperture.",
    "Simulation through space harmonics of zero length deflecting cavity.",
    } ;

QUAD quad_example;
/* quadrupole physical parameters */
PARAMETER quad_param[N_QUAD_PARAMS]={
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&quad_example.length), NULL, 0.0, 0, "length"},
    {"K1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.k1), NULL, 0.0, 0, "geometric strength"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"HKICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&quad_example.xkick), NULL, 0.0, 0, "horizontal correction kick"},
    {"VKICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&quad_example.ykick), NULL, 0.0, 0, "vertical correction kick"},
    {"HCALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.xKickCalibration), NULL, 1.0, 0, "calibration factor for horizontal correction kick"},
    {"VCALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.yKickCalibration), NULL, 1.0, 0, "calibration factor for vertical correction kick"},
    {"HSTEERING", "", IS_SHORT, 0, (long)((char *)&quad_example.xSteering), NULL, 0.0, 0, "use for horizontal steering?"},
    {"VSTEERING", "", IS_SHORT, 0, (long)((char *)&quad_example.ySteering), NULL, 0.0, 0, "use for vertical steering?"},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.order), NULL, 0.0, 0, "matrix order"},
    {"EDGE1_EFFECTS", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.edge1_effects), NULL, 0.0, 1, "include entrance edge effects?"},
    {"EDGE2_EFFECTS", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.edge2_effects), NULL, 0.0, 1, "include exit edge effects?"},
    {"FRINGE_TYPE", "", IS_STRING, 0, (long)((char *)&quad_example.fringeType), "fixed-strength", 0.0, 0, "type of fringe: \"inset\", \"fixed-strength\", or \"integrals\""},
    {"FFRINGE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.ffringe), NULL, 0.0, 0, "For non-integrals mode, fraction of length occupied by linear fringe region."},
    {"LEFFECTIVE", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&quad_example.lEffective), NULL, -1.0, 0, "Effective length. Ignored if non-positive. Cannot be used with non-zero FFRINGE."},
    {"I0P", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.fringeIntP[0]), NULL, 0.0, 0, "i0+ fringe integral"},
    {"I1P", "M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.fringeIntP[1]), NULL, 0.0, 0, "i1+ fringe integral"},
    {"I2P", "M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.fringeIntP[2]), NULL, 0.0, 0, "i2+ fringe integral"},
    {"I3P", "M$a4$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.fringeIntP[3]), NULL, 0.0, 0, "i3+ fringe integral"},
    {"LAMBDA2P", "M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.fringeIntP[4]), NULL, 0.0, 0, "lambda2+ fringe integral"},
    {"I0M", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.fringeIntM[0]), NULL, 0.0, 0, "i0- fringe integral"},
    {"I1M", "M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.fringeIntM[1]), NULL, 0.0, 0, "i1- fringe integral"},
    {"I2M", "M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.fringeIntM[2]), NULL, 0.0, 0, "i2- fringe integral"},
    {"I3M", "M$a4$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.fringeIntM[3]), NULL, 0.0, 0, "i3- fringe integral"},
    {"LAMBDA2M", "M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.fringeIntM[4]), NULL, 0.0, 0, "lambda2- fringe integral"},
    {"RADIAL", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&quad_example.radial), NULL, 0.0, 0, "If non-zero, converts the quadrupole into a radially-focusing lens"},

    };

BEND sbend_example;
/* bending magnet physical parameters */
PARAMETER sbend_param[N_BEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&sbend_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&sbend_example.angle), NULL, 0.0, 0, "bend angle"},
    {"K1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.k1), NULL, 0.0, 0, "geometric focusing strength"},
    {"E1", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.e[0]), NULL, 0.0, 0, "entrance edge angle"},
    {"E2", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.e[1]), NULL, 0.0, 0, "exit edge angle"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"K2", "1/M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.k2), NULL, 0.0, 0, "geometric sextupole strength"},
    {"H1", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.h[0]), NULL, 0.0, 0, "entrance pole-face curvature"},
    {"H2", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.h[1]), NULL, 0.0, 0, "exit pole-face curvature"},
    {"HGAP", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.hgap), NULL, 0.0, 0, "half-gap between poles"},
    {"FINT", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.fint), NULL, DEFAULT_FINT, 0, "edge-field integral"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.dx), NULL, 0.0, 0, "misaligment of entrance"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.dy), NULL, 0.0, 0, "misalignment of entrance"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.dz), NULL, 0.0, 0, "misalignment of entrance"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.fse), NULL, 0.0, 0, "fractional strength error of all components"},
    {"FSE_DIPOLE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.fseDipole), NULL, 0.0, 0, "fractional strength error of dipole component"},
    {"FSE_QUADRUPOLE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.fseQuadrupole), NULL, 0.0, 0, "fractional strength error of quadrupole component"},
    {"ETILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.etilt), NULL, 0.0, 0, "error rotation about incoming longitudinal axis"},
    {"EDGE1_EFFECTS", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.edge_effects[0]), NULL, 0.0, 1, "include entrance edge effects?"},
    {"EDGE2_EFFECTS", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.edge_effects[1]), NULL, 0.0, 1, "include exit edge effects?"},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.order), NULL, 0.0, 0, "matrix order"},
    {"EDGE_ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.edge_order), NULL, 0.0, 0, "edge matrix order"},
    {"TRANSPORT", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.TRANSPORT), NULL, 0.0, 0, "use (incorrect) TRANSPORT equations for T436 of edge?"},
    {"USE_BN", "", IS_SHORT, 0, (long)((char *)&sbend_example.use_bn), NULL, 0.0, 0, "use B1 and B2 instead of K1 and K2 values?"},
    {"B1", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.b1), NULL, 0.0, 0, "K1 = B1/rho, where rho is bend radius"},
    {"B2", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sbend_example.b2), NULL, 0.0, 0, "K2 = B2/rho"},
    };

BEND rbend_example;
/* rectangular bending magnet physical parameters */
/* The only difference between this and S-bend is in the description of L */
PARAMETER rbend_param[N_BEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&rbend_example.length), NULL, 0.0, 0, "magnet (straight) length"},
    {"ANGLE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&rbend_example.angle), NULL, 0.0, 0, "bend angle"},
    {"K1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.k1), NULL, 0.0, 0, "geometric focusing strength"},
    {"E1", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.e[0]), NULL, 0.0, 0, "entrance edge angle"},
    {"E2", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.e[1]), NULL, 0.0, 0, "exit edge angle"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"K2", "1/M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.k2), NULL, 0.0, 0, "geometric sextupole strength"},
    {"H1", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.h[0]), NULL, 0.0, 0, "entrance pole-face curvature"},
    {"H2", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.h[1]), NULL, 0.0, 0, "exit pole-face curvature"},
    {"HGAP", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.hgap), NULL, 0.0, 0, "half-gap between poles"},
    {"FINT", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.fint), NULL, DEFAULT_FINT, 0, "edge-field integral"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.dx), NULL, 0.0, 0, "misaligment of entrance"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.dy), NULL, 0.0, 0, "misalignment of entrance"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.dz), NULL, 0.0, 0, "misalignment of entrance"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.fse), NULL, 0.0, 0, "fractional strength error of all components"},
    {"FSE_DIPOLE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.fseDipole), NULL, 0.0, 0, "fractional strength error of dipole component"},
    {"FSE_QUADRUPOLE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.fseQuadrupole), NULL, 0.0, 0, "fractional strength error of quadrupole component"},
    {"ETILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.etilt), NULL, 0.0, 0, "error rotation about incoming longitudinal axis"},
    {"EDGE1_EFFECTS", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.edge_effects[0]), NULL, 0.0, 1, "include entrance edge effects?"},
    {"EDGE2_EFFECTS", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.edge_effects[1]), NULL, 0.0, 1, "include exit edge effects?"},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.order), NULL, 0.0, 0, "matrix order"},
    {"EDGE_ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.edge_order), NULL, 0.0, 0, "edge matrix order"},
    {"TRANSPORT", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.TRANSPORT), NULL, 0.0, 0, "use (incorrect) TRANSPORT equations for T436 of edge?"},
    {"USE_BN", "", IS_SHORT, 0, (long)((char *)&rbend_example.use_bn), NULL, 0.0, 0, "use B1 and B2 instead of K1 and K2 values?"},
    {"B1", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.b1), NULL, 0.0, 0, "K1 = B1/rho, where rho is bend radius"},
    {"B2", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rbend_example.b2), NULL, 0.0, 0, "K2 = B2/rho"},
    };

DRIFT drift_example;
/* drift length physical parameters */
PARAMETER drift_param[N_DRIFT_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&drift_example.length), NULL, 0.0, 0, "length"},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&drift_example.order), NULL, 0.0, 0, "matrix order"}
    };

SEXT sext_example;
/* sextupole physical parameters */
PARAMETER sext_param[N_SEXT_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&sext_example.length), NULL, 0.0, 0, "length"},
    {"K2", "1/M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sext_example.k2), NULL, 0.0, 0, "geometric strength"},
    {"K1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sext_example.k1), NULL, 0.0, 0, "geometric quadrupole strength error. See notes below!"},
    {"J1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sext_example.j1), NULL, 0.0, 0, "geometric skew quadrupole strength error. See notes below!"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sext_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sext_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sext_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sext_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sext_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"FFRINGE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sext_example.ffringe), NULL, 0.0, 0, "Length occupied by linear fringe regions as fraction hard-edge length L."},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&sext_example.order), NULL, 0.0, 0, "matrix order"},
    };
   
OCTU oct_example;
/* octupole physical parameters */
PARAMETER octu_param[N_OCTU_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&oct_example.length), NULL, 0.0, 0, "length"},
    {"K3", "1/M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&oct_example.k3), NULL, 0.0, 0, "geometric strength"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&oct_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&oct_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&oct_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&oct_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&oct_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&oct_example.order), NULL, 0.0, 0, "matrix order"}
    };
   
MULT mult_example;
/* multipole physical parameters */
PARAMETER mult_param[N_MULT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.length), NULL, 0.0, 0, "length"},
    {"KNL", "M$a-ORDER$n", IS_DOUBLE, 0, (long)((char *)&mult_example.KnL), NULL, 0.0, 0, "integrated geometric strength"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&mult_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"BORE", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.bore), NULL, 0.0, 0, "bore radius"},
    {"BTIPL", "T M", IS_DOUBLE, 0, (long)((char *)&mult_example.BTipL), NULL, 0.0, 0, "integrated field at pole tip, used if BORE nonzero"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FACTOR", "", IS_DOUBLE, 0, (long)((char *)&mult_example.factor), NULL, 1.0, 0, "factor by which to multiply strength"},
    {"ORDER", "", IS_SHORT, 0, (long)((char *)&mult_example.order), NULL, 0.0, 1, "multipole order"},
    {"N_KICKS", "", IS_SHORT, 0, (long)((char *)&mult_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks"},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&mult_example.synch_rad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"EXPAND_HAMILTONIAN", "", IS_SHORT, 0, (long)((char *)&mult_example.expandHamiltonian), NULL, 0.0, 0, "If 1, Hamiltonian is expanded to leading order."},
    };

FMULT fmult_example;
/* multipole physical parameters */
PARAMETER fmult_param[N_FMULT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&fmult_example.length), NULL, 0.0, 0, "length"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&fmult_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&fmult_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&fmult_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 0, (long)((char *)&fmult_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, 0, (long)((char *)&fmult_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&fmult_example.n_kicks), NULL, 0.0, 1, "number of kicks"},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&fmult_example.synch_rad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&fmult_example.filename), NULL, 0.0, 0, "name of file containing multipole data"},
    {"SQRT_ORDER", "", IS_SHORT, 0, (long)((char *)&fmult_example.sqrtOrder), NULL, 0.0, 0, "Ignored, kept for backward compatibility only."},
    };

TAYLORSERIES taylorSeries_example;
/* taylor series physical parameters */
PARAMETER taylorSeries_param[N_TAYLORSERIES_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&taylorSeries_example.length), NULL, 0.0, 0, "length"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&taylorSeries_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&taylorSeries_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&taylorSeries_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 0, (long)((char *)&taylorSeries_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&taylorSeries_example.filename), NULL, 0.0, 0, "name of file containing talor series map data"},
    };

SOLE sole_example;
/* solenoid physical parameters */
PARAMETER sole_param[N_SOLE_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&sole_example.length), NULL, 0.0, 0, "length"},
    {"KS", "RAD/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sole_example.ks), NULL, 0.0, 0, "geometric strength, -Bs/(B*Rho)"},
    {"B", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sole_example.B), NULL, 0.0, 0, "field strength (used if KS is zero)"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sole_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sole_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&sole_example.dz), NULL, 0.0, 0, "misalignment"},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&sole_example.order), NULL, 0.0, 0, "matrix order"},
    };
   
HCOR hcor_example;
/* horizontal corrector physical parameters */
PARAMETER hcor_param[N_HCOR_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hcor_example.length), NULL, 0.0, 0, "length"},
    {"KICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hcor_example.kick), NULL, 0.0, 0, "kick strength"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hcor_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"B2", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hcor_example.b2), NULL, 0.0, 0, "normalized sextupole strength (kick = KICK*(1+B2*x^2) when y=0)"},
    {"CALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hcor_example.calibration), NULL, 1.0, 0, "strength multiplier"},
    {"EDGE_EFFECTS", "", IS_SHORT, 0, (long)((char *)&hcor_example.edge_effects), NULL, 0.0, 0, "include edge effects?"},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&hcor_example.order), NULL, 0.0, 0, "matrix order"},
    {"STEERING", "", IS_SHORT, 0, (long)((char *)&hcor_example.steering), NULL, 0.0, 1, "use for steering?"},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&hcor_example.synchRad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&hcor_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"LERAD", "", IS_DOUBLE, 0, (long)((char *)&hcor_example.lEffRad), NULL, 0.0, 0, "if L=0, use this length for radiation computations"},
    };

VCOR vcor_example;
/* vertical corrector physical parameters */
PARAMETER vcor_param[N_VCOR_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&vcor_example.length), NULL, 0.0, 0, "length"},
    {"KICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&vcor_example.kick), NULL, 0.0, 0, "kick strength"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&vcor_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"B2", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&vcor_example.b2), NULL, 0.0, 0, "normalized sextupole strength (kick = KICK*(1+B2*y^2))"},
    {"CALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&vcor_example.calibration), NULL, 1.0, 0, "strength multiplier"},
    {"EDGE_EFFECTS", "", IS_SHORT, 0, (long)((char *)&vcor_example.edge_effects), NULL, 0.0, 0, "include edge effects?"},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&vcor_example.order), NULL, 0.0, 0, "matrix order"},
    {"STEERING", "", IS_SHORT, 0, (long)((char *)&vcor_example.steering), NULL, 0.0, 1, "use for steering?"},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&vcor_example.synchRad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&vcor_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"LERAD", "", IS_DOUBLE, 0, (long)((char *)&vcor_example.lEffRad), NULL, 0.0, 0, "if L=0, use this length for radiation computations"},
    };

RFCA rfca_example;
/* rf cavity physical parameters */
PARAMETER rfca_param[N_RFCA_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfca_example.length), NULL, 0.0, 0, "length"},
    {"VOLT", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfca_example.volt), NULL, 0.0, 0, "peak voltage"},
    {"PHASE", "DEG", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfca_example.phase), NULL, 0.0, 0, "phase"},
    {"FREQ", "Hz", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfca_example.freq), NULL, 500.0e6, 0, "frequency"},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&rfca_example.Q), NULL, 0.0, 0, "cavity Q (for cavity that charges up to given voltage from 0)"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&rfca_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"CHANGE_P0", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&rfca_example.change_p0), NULL, 0.0, 0, "does cavity change central momentum?"}, 
    {"CHANGE_T", "", IS_SHORT, 0, (long)((char *)&rfca_example.change_t), NULL, 0.0, 0, "set to 1 for long runs to avoid rounding error in phase"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&rfca_example.fiducial), NULL, 0.0, 0, "mode for determining fiducial arrival time (light, tmean, first, pmaximum)"},
    {"END1_FOCUS", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&rfca_example.end1Focus), NULL, 0.0, 0, "include focusing at entrance?"},
    {"END2_FOCUS", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&rfca_example.end2Focus), NULL, 0.0, 0, "include focusing at exit?"},
    {"BODY_FOCUS_MODEL", "", IS_STRING, PARAM_CHANGES_MATRIX, (long)((char *)&rfca_example.bodyFocusModel), NULL, 0.0, 0, "None (default) or SRS (simplified Rosenzweig/Serafini for standing wave)"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&rfca_example.nKicks), NULL, 0.0, 0, "Number of kicks to use for kick method.  Set to zero for matrix method."},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfca_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfca_example.dy), NULL, 0.0, 0, "misalignment"},
    {"T_REFERENCE", "S", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfca_example.tReference), NULL, -1.0, 0, "arrival time of reference particle"},
    {"LINEARIZE", "", IS_SHORT, 0, (long)((char *)&rfca_example.linearize), NULL, 0.0, 0, "Linearize phase dependence?"},
    {"LOCK_PHASE", "", IS_SHORT, 0, (long)((char *)&rfca_example.lockPhase), NULL, 0.0, 0, "Lock phase to given value regardless of bunch centroid motion?"},
    };
   
HMON hmon_example;
/* name for horizontal monitor physical parameters */
PARAMETER hmon_param[N_HMON_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hmon_example.length), NULL, 0.0, 0, "length"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hmon_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hmon_example.dy), NULL, 0.0, 0, "misalignment"},
    {"WEIGHT", "", IS_DOUBLE, 0, (long)((char *)&hmon_example.weight), NULL, 1.0, 0, "weight in correction"},
    {"TILT", "", IS_DOUBLE, 0, (long)((char *)&hmon_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"CALIBRATION", "", IS_DOUBLE, 0, (long)((char *)&hmon_example.calibration), NULL, 1.0, 0, "calibration factor for readout"},
    {"SETPOINT", "M", IS_DOUBLE, 0, (long)((char *)&hmon_example.setpoint), NULL, 0.0, 0, "steering setpoint"},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&hmon_example.order), NULL, 1.0, 0, "matrix order"},
    {"READOUT", "", IS_STRING, 0, (long)((char *)&hmon_example.readout), NULL, 0.0, 1, "rpn expression for readout (actual position supplied in variable x)"},
    {"CO_FITPOINT", "", IS_SHORT, 0, (long)((char *)&hmon_example.coFitpoint), NULL, 0.0, 0, "If nonzero, then closed orbit value is placed in variable <name>\\#<occurence>.xco"},
    } ;
   
VMON vmon_example;
/* name for vertical monitor physical parameters */
PARAMETER vmon_param[N_VMON_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&vmon_example.length), NULL, 0.0, 0, "length"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&vmon_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&vmon_example.dy), NULL, 0.0, 0, "misalignment"},
    {"WEIGHT", "", IS_DOUBLE, 0, (long)((char *)&vmon_example.weight), NULL, 1.0, 0, "weight in correction"},
    {"TILT", "", IS_DOUBLE, 0, (long)((char *)&vmon_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"CALIBRATION", "", IS_DOUBLE, 0, (long)((char *)&vmon_example.calibration), NULL, 1.0, 0, "calibration factor for readout"},
    {"SETPOINT", "M", IS_DOUBLE, 0, (long)((char *)&vmon_example.setpoint), NULL, 0.0, 0, "steering setpoint"},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&vmon_example.order), NULL, 1.0, 0, "matrix order"},
    {"READOUT", "", IS_STRING, 0, (long)((char *)&vmon_example.readout), NULL, 0.0, 1, "rpn expression for readout (actual position supplied in variable y)"},
    {"CO_FITPOINT", "", IS_SHORT, 0, (long)((char *)&vmon_example.coFitpoint), NULL, 0.0, 0, "If nonzero, then closed orbit value is placed in variable <name>\\#<occurence>.yco"},
    } ;

MONI moni_example;   
/* name for two-plane monitor physical parameters */
PARAMETER moni_param[N_MONI_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&moni_example.length), NULL, 0.0, 0, "length"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&moni_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&moni_example.dy), NULL, 0.0, 0, "misalignment"},
    {"WEIGHT", "", IS_DOUBLE, 0, (long)((char *)&moni_example.weight), NULL, 1.0, 0, "weight in correction"},
    {"TILT", "", IS_DOUBLE, 0, (long)((char *)&moni_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"XCALIBRATION", "", IS_DOUBLE, 0, (long)((char *)&moni_example.xcalibration), NULL, 1.0, 0, "calibration factor for x readout"},
    {"YCALIBRATION", "", IS_DOUBLE, 0, (long)((char *)&moni_example.ycalibration), NULL, 1.0, 0, "calibration factor for y readout"},
    {"XSETPOINT", "M", IS_DOUBLE, 0, (long)((char *)&moni_example.xsetpoint), NULL, 0.0, 0, "x steering setpoint"},
    {"YSETPOINT", "M", IS_DOUBLE, 0, (long)((char *)&moni_example.ysetpoint), NULL, 0.0, 0, "y steering setpoint"},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&moni_example.order), NULL, 1.0, 0, "matrix order"},
    {"XREADOUT", "", IS_STRING, 0, (long)((char *)&moni_example.x_readout), NULL, 0.0, 1, "rpn expression for x readout (actual position supplied in variables x, y"},
    {"YREADOUT", "", IS_STRING, 0, (long)((char *)&moni_example.y_readout), NULL, 0.0, 1, "rpn expression for y readout (actual position supplied in variables x, y"},
    {"CO_FITPOINT", "", IS_SHORT, 0, (long)((char *)&moni_example.coFitpoint), NULL, 0.0, 0, "If nonzero, then closed orbit values are placed in variables <name>\\#<occurence>.xco and <name>\\#<occurence>.yco"},
    } ;

RCOL rcol_example;   
/* name for rectangular collimator physical parameters */
PARAMETER rcol_param[N_RCOL_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rcol_example.length), NULL, 0.0, 0, "length"},
    {"X_MAX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rcol_example.x_max), NULL, 0.0, 0, "half-width in x"},
    {"Y_MAX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rcol_example.y_max), NULL, 0.0, 0, "half-width in y"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rcol_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rcol_example.dy), NULL, 0.0, 0, "misalignment"},
    {"OPEN_SIDE", "", IS_STRING, 0, (long)((char *)&rcol_example.openSide), NULL, 0.0, 0, "which side, if any, is open (+x, -x, +y, -y)"},
    {"INVERT", "", IS_SHORT, 0, (long)((char *)&rcol_example.invert), NULL, 0.0, 0, "If non-zero, particles inside the aperture are lost while those outside are transmitted."},
    } ;
   
ECOL ecol_example;
/* name for elliptical collimator physical parameters */
PARAMETER ecol_param[N_ECOL_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ecol_example.length), NULL, 0.0, 0, "length"},
    {"X_MAX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ecol_example.x_max), NULL, 0.0, 0, "half-axis in x"},
    {"Y_MAX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ecol_example.y_max), NULL, 0.0, 0, "half-axis in y"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ecol_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ecol_example.dy), NULL, 0.0, 0, "misalignment"},
    {"OPEN_SIDE", "", IS_STRING, 0, (long)((char *)&ecol_example.openSide), NULL, 0.0, 0, "which side, if any, is open (+x, -x, +y, -y)"},
    {"EXPONENT", "", IS_SHORT, 0, (long)((char *)&ecol_example.exponent), NULL, 0.0, 2, "Exponent for boundary equation.  2 is ellipse."},
    {"YEXPONENT", "", IS_SHORT, 0, (long)((char *)&ecol_example.yExponent), NULL, 0.0, 0, "y exponent for boundary equation.  2 is ellipse.  If 0, defaults to EXPONENT"},
    {"INVERT", "", IS_SHORT, 0, (long)((char *)&ecol_example.invert), NULL, 0.0, 0, "If non-zero, particles inside the aperture are lost while those outside are transmitted."},
    } ;

CLEAN clean_example;
/* name for beam cleaner physical parameters */
PARAMETER clean_param[N_CLEAN_PARAMS] = {
    {"MODE", "", IS_STRING, 0, (long)((char *)&clean_example.mode), "stdeviation", 0.0, 0, "stdeviation, absdeviation, or absvalue"},
    {"XLIMIT", "", IS_DOUBLE, 0, (long)((char *)&clean_example.xLimit), NULL, 0.0, 0, "Limit for x"},
    {"XPLIMIT", "", IS_DOUBLE, 0, (long)((char *)&clean_example.xpLimit), NULL, 0.0, 0, "Limit for x'"},
    {"YLIMIT", "", IS_DOUBLE, 0, (long)((char *)&clean_example.yLimit), NULL, 0.0, 0, "Limit for y"},
    {"YPLIMIT", "", IS_DOUBLE, 0, (long)((char *)&clean_example.ypLimit), NULL, 0.0, 0, "Limit for y'"},
    {"TLIMIT", "", IS_DOUBLE, 0, (long)((char *)&clean_example.tLimit), NULL, 0.0, 0, "Limit for t"},
    {"DELTALIMIT", "", IS_DOUBLE, 0, (long)((char *)&clean_example.deltaLimit), NULL, 0.0, 0, "Limit for (p-p0)/p0"},
    } ;

MARK mark_example;
/* name for marker parameters */
PARAMETER mark_param[N_MARK_PARAMS] = {
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&mark_example.dx), NULL, 0.0, 0, "non-functional misalignment (e.g., for girder)"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&mark_example.dy), NULL, 0.0, 0, "non-functional misalignment (e.g., for girder)"},
    {"FITPOINT", "", IS_SHORT, 0, (long)((char *)&mark_example.fitpoint), NULL, 0.0, 0, "Supply local values of Twiss parameters, moments, floor coordinates, matrices, etc. for optimization?"},
    } ;

MATR matr_example;
/* name for matrix parameters */
PARAMETER matr_param[N_MATR_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&matr_example.length), NULL, 0.0, 0, "length"},
    {"FRACTION", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&matr_example.fraction), NULL, 1.0, 0, "The provided matrix M is interpolated with the identity matrix I according to f*M+(1-f)*I."},
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&matr_example.filename), "", 0.0, 0, "input file"},
    {"ORDER", "", IS_SHORT, 0, (long)((char *)&matr_example.order), NULL, 0.0, 1, "matrix order"},
    } ;

ALPH alph_example;
/* names for alpha magnet parameters */
PARAMETER alph_param[N_ALPH_PARAMS] = {
    {"XMAX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&alph_example.xmax), NULL, 0.0, 0, "size of alpha"},
    {"XS1", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&alph_example.xs1), NULL, 0.0, 0, "inner scraper position relative to XMAX"},
    {"XS2", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&alph_example.xs2), NULL, 0.0, 0, "outer scraper position relative to XMAX"},
    {"DP1", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&alph_example.dp1), NULL, -1, 0, "inner scraper fractional momentum deviation"},
    {"DP2", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&alph_example.dp2), NULL, 1, 0, "outer scraper fractional momentum deviation"},
    {"XPUCK", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&alph_example.xPuck), NULL, -1, 0, "position of scraper puck"},
    {"WIDTHPUCK", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&alph_example.widthPuck), NULL, 0, 0, "size of scraper puck"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&alph_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&alph_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&alph_example.dz), NULL, 0.0, 0, "misalignment"},
    {"TILT", "", IS_DOUBLE, 0, (long)((char *)&alph_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"PART", "", IS_SHORT, 0, (long)((char *)&alph_example.part), NULL, 0.0, 0, "0=full, 1=first half, 2=second half"},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&alph_example.order), NULL, 0.0, 0, "matrix order [1,3]"}
    } ;

RFDF rfdf_example;
/* names for rf deflector parameters */
PARAMETER rfdf_param[N_RFDF_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.length), NULL, 0.0, 0, "length"},
    {"PHASE", "DEG", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.phase), NULL, 0.0, 0, "phase"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"FREQUENCY", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.frequency), NULL, DEFAULT_FREQUENCY, 0, "frequency"},
    {"VOLTAGE", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.voltage), NULL, 0.0, 0, "voltage"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.fse), NULL, 0.0, 0, "Fractional Strength Error"},
    {"B2", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.b2), NULL, 0.0, 0, "Normalized sextupole strength, kick=(1+b2*(x^2-y^2)/2)..."},
    {"TIME_OFFSET", "S", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.time_offset), NULL, 0.0, 0, "time offset (adds to phase)"},
    {"N_KICKS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.n_kicks), NULL, 0.0, 0, "number of kicks (0=autoscale)"},
    {"PHASE_REFERENCE", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"STANDING_WAVE", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.standingWave), NULL, 0.0, 0, "If nonzero, then cavity is standing wave."},
    {"VOLTAGE_WAVEFORM", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&rfdf_example.voltageWaveform), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving voltage waveform factor vs time"},
    {"VOLTAGE_PERIODIC", "", IS_SHORT, 0, (long)((char *)&rfdf_example.voltageIsPeriodic), NULL, 0.0, 0, "If non-zero, voltage waveform is periodic with period given by time span."},
    {"ALIGN_WAVEFORMS", "", IS_SHORT, 0,  (long)((char *)&rfdf_example.alignWaveforms), NULL, 0.0, 0, "If non-zero, waveforms' t=0 is aligned with first bunch arrival time."},
    {"VOLTAGE_NOISE", "", IS_DOUBLE, 0, (long)((char *)&rfdf_example.voltageNoise), NULL, 0.0, 0, "Rms fractional noise level for voltage."},
    {"PHASE_NOISE", "DEG", IS_DOUBLE, 0, (long)((char *)&rfdf_example.phaseNoise), NULL, 0.0, 0, "Rms noise level for phase."},
    {"GROUP_VOLTAGE_NOISE", "", IS_DOUBLE, 0, (long)((char *)&rfdf_example.groupVoltageNoise), NULL, 0.0, 0, "Rms fractional noise level for voltage linked to group."},
    {"GROUP_PHASE_NOISE", "DEG", IS_DOUBLE, 0, (long)((char *)&rfdf_example.groupPhaseNoise), NULL, 0.0, 0, "Rms noise level for phase linked to group."},
    {"VOLTAGE_NOISE_GROUP", "", IS_LONG, 0, (long)((char *)&rfdf_example.voltageNoiseGroup), NULL, 0.0, 0, "Group number for voltage noise."},
    {"PHASE_NOISE_GROUP", "", IS_LONG, 0, (long)((char *)&rfdf_example.phaseNoiseGroup), NULL, 0.0, 0, "Group number for phase noise."},
    {"START_PASS", "", IS_LONG, 0, (long)((char *)&rfdf_example.startPass), NULL, 0.0, -1, "If non-negative, pass on which to start modeling cavity."},    
    {"END_PASS", "", IS_LONG, 0, (long)((char *)&rfdf_example.endPass), NULL, 0.0, -1, "If non-negative, pass on which to end modeling cavity."},
    {"DRIFT_MATRIX", "", IS_SHORT, 0, (long)((char *)&rfdf_example.driftMatrix), NULL, 0.0, 0, "If non-zero, calculations involving matrices assume this element is a drift space."},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfdf_example.dz), NULL, 0.0, 0, "misalignment"},
    {"MAGNETIC_DEFLECTION", "", IS_SHORT, 0, (long)((char *)&rfdf_example.magneticDeflection), NULL, 0.0, 0, "If non-zero, deflection is assumed to be performed by a magnetic field, rather than electric field (default)."},
    } ;

RFTM110 rftm110_example;
/* names for rf tm110 deflecting cavity parameters */
PARAMETER rftm110_param[N_RFTM110_PARAMS] = {
    {"PHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&rftm110_example.phase), NULL, 0.0, 0, "phase"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&rftm110_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"FREQUENCY", "HZ", IS_DOUBLE, 0, (long)((char *)&rftm110_example.frequency), NULL, DEFAULT_FREQUENCY, 0, "frequency"},
    {"VOLTAGE", "V", IS_DOUBLE, 0, (long)((char *)&rftm110_example.voltage), NULL, 0.0, 0, "peak deflecting voltage"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&rftm110_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"VOLTAGE_WAVEFORM", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&rftm110_example.voltageWaveform), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving voltage waveform factor vs time"},
    {"VOLTAGE_PERIODIC", "", IS_SHORT, 0, (long)((char *)&rftm110_example.voltageIsPeriodic), NULL, 0.0, 0, "If non-zero, voltage waveform is periodic with period given by time span."},
    {"ALIGN_WAVEFORMS", "", IS_SHORT, 0,  (long)((char *)&rftm110_example.alignWaveforms), NULL, 0.0, 0, "If non-zero, waveforms' t=0 is aligned with first bunch arrival time."},
    {"VOLTAGE_NOISE", "", IS_DOUBLE, 0, (long)((char *)&rftm110_example.voltageNoise), NULL, 0.0, 0, "Rms fractional noise level for voltage."},
    {"PHASE_NOISE", "DEG", IS_DOUBLE, 0, (long)((char *)&rftm110_example.phaseNoise), NULL, 0.0, 0, "Rms noise level for phase."},
    {"GROUP_VOLTAGE_NOISE", "", IS_DOUBLE, 0, (long)((char *)&rftm110_example.groupVoltageNoise), NULL, 0.0, 0, "Rms fractional noise level for voltage linked to group."},
    {"GROUP_PHASE_NOISE", "DEG", IS_DOUBLE, 0, (long)((char *)&rftm110_example.groupPhaseNoise), NULL, 0.0, 0, "Rms noise level for phase linked to group."},
    {"VOLTAGE_NOISE_GROUP", "", IS_LONG, 0, (long)((char *)&rftm110_example.voltageNoiseGroup), NULL, 0.0, 0, "Group number for voltage noise."},
    {"PHASE_NOISE_GROUP", "", IS_LONG, 0, (long)((char *)&rftm110_example.phaseNoiseGroup), NULL, 0.0, 0, "Group number for phase noise."},
    {"START_PASS", "", IS_LONG, 0, (long)((char *)&rftm110_example.startPass), NULL, 0.0, -1, "If non-negative, pass on which to start modeling cavity."},    
    {"END_PASS", "", IS_LONG, 0, (long)((char *)&rftm110_example.endPass), NULL, 0.0, -1, "If non-negative, pass on which to end modeling cavity."},    
    } ;

RFTMEZ0 rftmez0_example;
/* names for tm mode cavity from Ez(z,r=0)
 */
PARAMETER rftmez0_param[N_RFTMEZ0_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.length), NULL, 0.0, 0, "length"},
    {"FREQUENCY", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.frequency), NULL, DEFAULT_FREQUENCY, 0, "frequency"},
    {"PHASE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.phase), NULL, 0.0, 0, "phase"},
    {"EZ_PEAK", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.Ez_peak), NULL, 0.0, 0, "Peak on-axis longitudinal electric field"},
    {"TIME_OFFSET", "S", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.time_offset), NULL, 0.0, 0, "time offset (adds to phase)"},
    {"PHASE_REFERENCE", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link to other time-dependent elements)"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.dzMA), NULL, 0.0, 0, "misalignment"},
    {"ETILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.eTilt), NULL, 0.0, 0, "misalignment"},
    {"EYAW", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.eYaw), NULL, 0.0, 0, "misalignment"},
    {"EPITCH", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.ePitch), NULL, 0.0, 0, "misalignment"},
    {"N_STEPS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.n_steps), NULL, 0.0, 100, "number of steps (for nonadaptive integration)"},
    {"RADIAL_ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.radial_order), NULL, 0.0, 1, "highest order in off-axis expansion"},
    {"CHANGE_P0", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.change_p0), NULL, 0.0, 0, "does element change central momentum?"},
    {"INPUTFILE", "", IS_STRING, 0, (long)((char *)&rftmez0_example.inputFile), NULL, 0.0, 0, "file containing Ez vs z at r=0"},
    {"ZCOLUMN", "", IS_STRING, 0, (long)((char *)&rftmez0_example.zColumn), NULL, 0.0, 0, "column containing z values"},
    {"EZCOLUMN", "", IS_STRING, 0, (long)((char *)&rftmez0_example.EzColumn), NULL, 0.0, 0, "column containing Ez values"},
    {"SOLENOID_FILE", "", IS_STRING, 0, (long)((char *)&rftmez0_example.solenoidFile), NULL, 0.0, 0, "file containing map of Bz and Br vs z and r.  Each page contains values for a single r."},
    {"SOLENOID_ZCOLUMN", "", IS_STRING, 0, (long)((char *)&rftmez0_example.solenoid_zColumn), NULL, 0.0, 0, "column containing z values for solenoid map."},
    {"SOLENOID_RCOLUMN", "", IS_STRING, 0, (long)((char *)&rftmez0_example.solenoid_rColumn), NULL, 0.0, 0, "column containing r values for solenoid map.  If omitted, data is assumed to be for r=0 and an on-axis expansion is performed."},
    {"SOLENOID_BZCOLUMN", "", IS_STRING, 0, (long)((char *)&rftmez0_example.solenoidBzColumn), NULL, 0.0, 0, "column containing Bz values for solenoid map."},
    {"SOLENOID_BRCOLUMN", "", IS_STRING, 0, (long)((char *)&rftmez0_example.solenoidBrColumn), NULL, 0.0, 0, "column containing Br values for solenoid map. If omitted, data is assumed to be for r=0 and an on-axis expansion is performed."},
    {"SOLENOID_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.solenoidFactor), NULL, 1.0, 0, "factor by which to multiply solenoid fields."},        
    {"SOLENOID_DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.dxSol), NULL, 0.0, 0, "misalignment"},
    {"SOLENOID_DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.dySol), NULL, 0.0, 0, "misalignment"},
    {"SOLENOID_DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.dzSolMA), NULL, 0.0, 0, "misalignment"},
    {"SOLENOID_ETILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.eTiltSol), NULL, 0.0, 0, "misalignment"},
    {"SOLENOID_EYAW", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.eYawSol), NULL, 0.0, 0, "misalignment"},
    {"SOLENOID_EPITCH", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.ePitchSol), NULL, 0.0, 0, "misalignment"},
    {"BX_STRAY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.BxStray), NULL, 0.0, 0, "Uniform stray horizontal field"},
    {"BY_STRAY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.ByStray), NULL, 0.0, 0, "Uniform stray vertical field"},
    {"ACCURACY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rftmez0_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&rftmez0_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, non-adaptive runge-kutta, modified midpoint)"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&rftmez0_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0, "{t|p},{median|min|max|ave|first|light} (e.g., \"t,median\")"},
    {"FIELD_TEST_FILE", "", IS_STRING, 0, (long)((char *)&rftmez0_example.fieldTestFile), NULL, 0.0, 0, "filename for output of test fields (r=0)"},
    } ;

RMDF rmdf_example;
/* names for rf deflector parameters */
PARAMETER rmdf_param[N_RMDF_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&rmdf_example.length), NULL, 0.0, 0, "length"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&rmdf_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"RAMP_TIME", "S", IS_DOUBLE, 0, (long)((char *)&rmdf_example.ramp_time), NULL, DEFAULT_RAMP_TIME, 0, "length of ramp"},
    {"VOLTAGE", "V", IS_DOUBLE, 0, (long)((char *)&rmdf_example.voltage), NULL, 0.0, 0, "full voltage"},
    {"GAP", "M", IS_DOUBLE, 0, (long)((char *)&rmdf_example.gap), NULL, DEFAULT_GAP, 0, "gap between plates"},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&rmdf_example.time_offset), NULL, 0.0, 0, "time offset of ramp start"},
    {"N_SECTIONS", "", IS_LONG, 0, (long)((char *)&rmdf_example.n_sections), NULL, 0.0, 10, "number of sections"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&rmdf_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&rmdf_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&rmdf_example.dy), NULL, 0.0, 0, "misalignment"},
    } ;

TMCF_MODE tmcf_example;
/* names for tm mode cavity with spatially constant fields, using nag 
 * integrator parameters 
 */
PARAMETER tmcf_param[N_TMCF_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&tmcf_example.length), NULL, 0.0, 0, "length"},
    {"FREQUENCY", "HZ", IS_DOUBLE, 0, (long)((char *)&tmcf_example.frequency), NULL, DEFAULT_FREQUENCY, 0, "frequency"},
    {"PHASE", "S", IS_DOUBLE, 0, (long)((char *)&tmcf_example.phase), NULL, 0.0, 0, "phase"},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&tmcf_example.time_offset), NULL, 0.0, 0, "time offset (adds to phase)"},
    {"RADIAL_OFFSET", "M", IS_DOUBLE, 0, (long)((char *)&tmcf_example.radial_offset), NULL, DEFAULT_RADIAL_OFFSET, 0, "not recommended"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&tmcf_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"ER", "V", IS_DOUBLE, 0, (long)((char *)&tmcf_example.Er), NULL, 0.0, 0, "radial electric field"},
    {"BPHI", "T", IS_DOUBLE, 0, (long)((char *)&tmcf_example.Bphi), NULL, 0.0, 0, "azimuthal magnetic field"},
    {"EZ", "V", IS_DOUBLE, 0, (long)((char *)&tmcf_example.Ez), NULL, 0.0, 0, "longitudinal electric field"},
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&tmcf_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&tmcf_example.x_max), NULL, 0.0, 0, "x half-aperture"},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&tmcf_example.y_max), NULL, 0.0, 0, "y half-aperture"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&tmcf_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&tmcf_example.dy), NULL, 0.0, 0, "misalignment"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&tmcf_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"N_STEPS", "", IS_LONG, 0, (long)((char *)&tmcf_example.n_steps), NULL, 0.0, 100, "number of steps (for nonadaptive integration)"},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&tmcf_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, non-adaptive runge-kutta, modified midpoint)"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&tmcf_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0, "{t|p},{median|min|max|ave|first|light} (e.g., \"t,median\")"}
    } ;

CE_PLATES cepl_example;
/* names for ramped deflector using nag integrator parameters 
 */
PARAMETER cepl_param[N_CEPL_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&cepl_example.length), NULL, 0.0, 0, "length"},
    {"RAMP_TIME", "S", IS_DOUBLE, 0, (long)((char *)&cepl_example.ramp_time), NULL, DEFAULT_RAMP_TIME, 0, "time to ramp to full strenth"},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&cepl_example.time_offset), NULL, 0.0, 0, "offset of ramp-start time"},
    {"VOLTAGE", "V", IS_DOUBLE, 0, (long)((char *)&cepl_example.voltage), NULL, 0.0, 0, "maximum voltage between plates due to ramp"},
    {"GAP", "M", IS_DOUBLE, 0, (long)((char *)&cepl_example.gap), NULL, DEFAULT_GAP, 0, "gap between plates"},
    {"STATIC_VOLTAGE", "V", IS_DOUBLE, 0, (long)((char *)&cepl_example.static_voltage), NULL, 0.0, 0, "static component of voltage"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&cepl_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&cepl_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&cepl_example.x_max), NULL, 0.0, 0, "x half-aperture"},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&cepl_example.y_max), NULL, 0.0, 0, "y half-aperture"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&cepl_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&cepl_example.dy), NULL, 0.0, 0, "misalignment"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&cepl_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"N_STEPS", "", IS_LONG, 0, (long)((char *)&cepl_example.n_steps), NULL, 0.0, 100, "number of steps (for nonadaptive integration)"},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&cepl_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, non-adaptive runge-kutta, modified midpoint)"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&cepl_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0, "{t|p},{median|min|max|ave|first|light} (e.g., \"t,median\")"}
    } ;

WATCH watch_example;
/* names for watch point */
PARAMETER watch_param[N_WATCH_PARAMS] = {
    {"FRACTION", "", IS_DOUBLE, 0, (long)((char *)&watch_example.fraction), NULL, 1.0, 0, "fraction of particles to dump (coordinate mode)"},
    {"START_PID", "", IS_LONG, 0, (long)((char *)&watch_example.startPID), NULL, 0.0, -1, "starting particleID for particles to dump"},
    {"END_PID", "", IS_LONG, 0, (long)((char *)&watch_example.endPID), NULL, 0.0, -1, "ending particleID for particles to dump"},
    {"INTERVAL", "", IS_LONG, 0, (long)((char *)&watch_example.interval), NULL, 0.0, 1, "interval for data output (in turns)"},
    {"START_PASS", "", IS_LONG, 0, (long)((char*)&watch_example.start_pass), NULL, 0.0, 0, "pass on which to start"},
    {"END_PASS", "", IS_LONG, 0, (long)((char*)&watch_example.end_pass), NULL, 0.0, -1, "pass on which to end (inclusive).  Ignored if negative."},
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&watch_example.filename), "", 0.0, 0, "output filename, possibly incomplete (see below)"},
    {"LABEL", "", IS_STRING, 0, (long)((char *)&watch_example.label), "", 0.0, 0, "output label"},
    {"MODE", "", IS_STRING, 0, (long)((char *)&watch_example.mode), "coordinates", 0.0, 0, "coordinate, parameter, centroid, or fft.  For fft mode, you may add a space and a qualifer giving the window type: hanning (default), parzen, welch, or uniform."},
    {"X_DATA", "", IS_SHORT, 0, (long)((char*)&watch_example.xData), NULL, 0.0, 1, "include x data in coordinate mode?"},
    {"Y_DATA", "", IS_SHORT, 0, (long)((char*)&watch_example.yData), NULL, 0.0, 1, "include y data in coordinate mode?"},
    {"LONGIT_DATA", "", IS_SHORT, 0, (long)((char*)&watch_example.longitData), NULL, 0.0, 1, "include longitudinal data in coordinate mode?"},
    {"EXCLUDE_SLOPES", "", IS_SHORT, 0, (long)((char*)&watch_example.excludeSlopes), NULL, 0.0, 0, "exclude slopes in coordinate mode?"},
    {"FLUSH_INTERVAL", "", IS_LONG, 0, (long)((char *)&watch_example.flushInterval), NULL, 0.0, 100, "file flushing interval (parameter or centroid mode)"},    
    {"SPARSE_INTERVAL", "", IS_LONG, 0, (long)((char *)&watch_example.sparseInterval), NULL, 0.0, 1, "interval for particle output (coordinate mode)"},
    {"DISABLE", "", IS_SHORT, 0, (long)((char *)&watch_example.disable), NULL, 0.0, 0, "If nonzero, no output will be generated."},    
    {"USE_DISCONNECT", "", IS_SHORT, 0, (long)((char *)&watch_example.useDisconnect), NULL, 0.0, 0, "If nonzero, files are disconnected between each write operation. May be useful for parallel operation.  Ignored otherwise."},
    {"INDEX_OFFSET", "", IS_LONG, 0, (long)((char *)&watch_example.indexOffset), NULL, 0.0, 0, "Offset for file indices for sequential file naming."},
    {"REFERENCE_FREQUENCY", "", IS_DOUBLE, 0, (long)((char *)&watch_example.referenceFrequency), NULL, -1.0, -1, "If non-zero, the indicated frequency is used to define the bucket center for purposes of computing time offsets."},
    } ;

TW_PLATES twpl_example;
/* names for tw (tem) ramped deflector using nag integrator parameters 
 */
PARAMETER twpl_param[N_TWPL_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&twpl_example.length), NULL, 0.0, 0, "length"},
    {"RAMP_TIME", "S", IS_DOUBLE, 0, (long)((char *)&twpl_example.ramp_time), NULL, DEFAULT_RAMP_TIME, 0, "time to ramp to full strenth"},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&twpl_example.time_offset), NULL, 0.0, 0, "offset of ramp-start time"},
    {"VOLTAGE", "V", IS_DOUBLE, 0, (long)((char *)&twpl_example.voltage), NULL, 0.0, 0, "maximum voltage between plates due to ramp"},
    {"GAP", "M", IS_DOUBLE, 0, (long)((char *)&twpl_example.gap), NULL, DEFAULT_GAP, 0, "gap between plates"},
    {"STATIC_VOLTAGE", "V", IS_DOUBLE, 0, (long)((char *)&twpl_example.static_voltage), NULL, 0.0, 0, "static component of voltage"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&twpl_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&twpl_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twpl_example.x_max), NULL, 0.0, 0, "x half-aperture"},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twpl_example.y_max), NULL, 0.0, 0, "y half-aperture"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&twpl_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&twpl_example.dy), NULL, 0.0, 0, "misalignment"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&twpl_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"N_STEPS", "", IS_LONG, 0, (long)((char *)&twpl_example.n_steps), NULL, 0.0, 100, "number of steps (for nonadaptive integration)"},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&twpl_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, non-adaptive runge-kutta, modified midpoint)"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&twpl_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0, "{t|p},{median|min|max|ave|first|light} (e.g., \"t,median\")"}
    } ;

MALIGN malign_example;
/* names for misaligment parameters */
PARAMETER malign_param[N_MALIGN_PARAMS] = {
    {"DXP", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&malign_example.dxp), NULL, 0.0, 0, "delta x'"},
    {"DYP", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&malign_example.dyp), NULL, 0.0, 0, "delta y'"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&malign_example.dx), NULL, 0.0, 0, "delta x"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&malign_example.dy), NULL, 0.0, 0, "delta y"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&malign_example.dz), NULL, 0.0, 0, "delta z"},
    {"DT", "S", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&malign_example.dt), NULL, 0.0, 0, "delta t"},
    {"DP", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&malign_example.dp), NULL, 0.0, 0, "delta p/pCentral"},
    {"DE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&malign_example.de), NULL, 0.0, 0, "delta gamma/gammaCentral"},
    {"ON_PASS", "", IS_LONG, 0, (long)((char *)&malign_example.on_pass), NULL, 0.0, -1, "pass on which to apply"},
    {"FORCE_MODIFY_MATRIX", "", IS_LONG, 0, (long)((char *)&malign_example.forceModifyMatrix), NULL, 0.0, 0, "modify the matrix even if on_pass>=0"},
    {"START_PID", "", IS_LONG, 0, (long)((char *)&malign_example.startPID), NULL, 0.0, -1, "starting particleID for particles to affect. By default, all particles are affected."},
    {"END_PID", "", IS_LONG, 0, (long)((char *)&malign_example.endPID), NULL, 0.0, -1, "ending particleID for particles to affect. By default, all particles are affected."},
    } ;

TW_LINAC twla_example;
/* names for traveling-wave linac parameters
 */
PARAMETER twla_param[N_TWLA_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.length), NULL, 0.0, 0, "length"},
    {"FREQUENCY", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.frequency), NULL, DEFAULT_FREQUENCY, 0, "frequency"},
    {"PHASE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.phase), NULL, 0.0, 0, "phase"},
    {"TIME_OFFSET", "S", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.time_offset), NULL, 0.0, 0, "time offset (adds to phase)"},
    {"EZ", "V/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.Ez), NULL, 0.0, 0, "electric field"},
    {"B_SOLENOID", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.B_solenoid), NULL, 0.0, 0, "solenoid field"},
    {"ACCURACY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twla_example.x_max), NULL, 0.0, 0, "x half-aperture"},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twla_example.y_max), NULL, 0.0, 0, "y half-aperture"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.dy), NULL, 0.0, 0, "misalignment"},
    {"BETA_WAVE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.beta_wave), NULL, DEFAULT_BETA_WAVE, 0, "(phase velocity)/c"},
    {"ALPHA", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.alpha), NULL, 0.0, 0, "field attenuation factor"},
    {"PHASE_REFERENCE", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"N_STEPS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.n_steps), NULL, 0.0, 100, "number of steps (for nonadaptive integration)"},
    {"FOCUSSING", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.focussing), NULL, 0.0, 1, "include focusing effects?"},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&twla_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, non-adaptive runge-kutta, modified midpoint)"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&twla_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0, "{t|p},{median|min|max|ave|first|light} (e.g., \"t,median\")"},
    {"CHANGE_P0", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.change_p0), NULL, 0.0, 0, "does element change central momentum?"},
    {"SUM_BN2", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twla_example.sum_bn2), NULL, 0.0, 0, "sum of squares of amplitudes of n!=0 space harmonics"},
    } ;

PEPPOT peppot_example;
/* names for pepper-pot plate */
PARAMETER peppot_param[N_PEPPOT_PARAMS] = {
    {"L"           , "M", IS_DOUBLE, 0, (long)((char *)&peppot_example.length), NULL, 0.0, 0, "length"},
    {"RADII"       , "M", IS_DOUBLE, 0, (long)((char *)&peppot_example.radii), NULL, 0.0, 0, "hole radius"},
    {"TRANSMISSION", "", IS_DOUBLE, 0, (long)((char *)&peppot_example.transmission), NULL, 0.0, 0, "transmission of material"},
    {"TILT"        , "RAD", IS_DOUBLE, 0, (long)((char *)&peppot_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"THETA_RMS"   , "RAD", IS_DOUBLE, 0, (long)((char *)&peppot_example.theta_rms), NULL, 0.0, 0, "rms scattering from material"},
    {"N_HOLES"     , "", IS_LONG  , 0, (long)((char *)&peppot_example.n_holes), NULL, 0.0, 0, "number of holes"},
    } ;

ENERGY energy_example;
/* names for energy */
PARAMETER energy_param[N_ENERGY_PARAMS] = {
    {"CENTRAL_ENERGY", "MC$a2$n", IS_DOUBLE, 0, (long)((char *)&energy_example.central_energy), NULL, 0.0, 0, "desired central gamma"},
    {"CENTRAL_MOMENTUM", "MC", IS_DOUBLE, 0, (long)((char *)&energy_example.central_momentum), NULL, 0.0, 0, "desired central beta*gamma"},
    {"MATCH_BEAMLINE", "", IS_LONG, 0, (long)((char *)&energy_example.match_beamline), NULL, 0.0, 0, "if nonzero, beamline reference momentum is set to beam average momentum"},
    {"MATCH_PARTICLES", "", IS_LONG, 0, (long)((char *)&energy_example.match_particles), NULL, 0.0, 0, "if nonzero, beam average momentum is set to beamline reference momentum"},
    } ;

MAXAMP maxamp_example;
/* names for max.amp. */
PARAMETER maxamp_param[N_MAXAMP_PARAMS] = {
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&maxamp_example.x_max), NULL, 0.0, 0, "x half-aperture"},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&maxamp_example.y_max), NULL, 0.0, 0, "y half-aperture"},
    {"ELLIPTICAL", "", IS_LONG, 0, (long)((char *)&maxamp_example.elliptical), NULL, 0.0, 0, "is aperture elliptical?"},
    {"EXPONENT", "", IS_LONG, 0, (long)((char *)&maxamp_example.exponent), NULL, 0.0, 2, "exponent for boundary equation in elliptical mode.  2 is a true ellipse."},
    {"YEXPONENT", "", IS_LONG, 0, (long)((char *)&maxamp_example.yExponent), NULL, 0.0, 0, "y exponent for boundary equation in elliptical mode.  If zero, defaults to EXPONENT."},
    {"OPEN_SIDE", "", IS_STRING, 0, (long)((char *)&maxamp_example.openSide), NULL, 0.0, 0, "which side, if any, is open (+x, -x, +y, -y)"},
    } ;

ROTATE rotate_example;
/* names for beam rotation */
PARAMETER rotate_param[N_ROTATE_PARAMS] = {
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rotate_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"EXCLUDE_FLOOR", "", IS_SHORT, 0, (long)((char *)&rotate_example.excludeFloor), NULL, 0.0, 0, "if non-zero, does not affect the floor coordinates"},
    {"EXCLUDE_OPTICS", "", IS_SHORT, 0, (long)((char *)&rotate_example.excludeOptics), NULL, 0.0, 0, "if non-zero, does not affect the optics (i.e., transfer matrix is unit matrix)"},
    } ;

/* names for transmission count */
PARAMETER trcount_param[N_TRCOUNT_PARAMS] = {
    {"DUMMY", "", IS_LONG, 0, 0, NULL, 0.0, 0, ""},
    } ;

/* names for reflection */
PARAMETER reflect_param[N_REFLECT_PARAMS] = {
    {"DUMMY", "", IS_LONG, 0, 0, NULL, 0.0, 0, ""},
    } ;

/* names for recirculation point */
PARAMETER recirc_param[N_RECIRC_PARAMS] = {
    {"I_RECIRC_ELEMENT", "", IS_LONG, 0, 0, NULL, 0.0, 0, ""},
    } ;

BRANCH branch_example;
/* names for branch instruction */
PARAMETER branch_param[N_BRANCH_PARAMS] = {
    {"COUNTER", "", IS_LONG, 0, (long)((char *)&branch_example.counter), NULL, 0.0, 0, "Counter, which is decremented by 1 for each pass. Set to negative value for unconditional branch."},
    {"INTERVAL", "", IS_LONG, 0, (long)((char *)&branch_example.interval), NULL, 0.0, 0, "Interval between branching. If non-positive, use COUNTER-based method instead."},
    {"OFFSET", "", IS_LONG, 0, (long)((char *)&branch_example.offset), NULL, 0.0, 0, "If INTERVAL method used, offset of branch passes."},
    {"VERBOSITY", "", IS_LONG, 0, (long)((char *)&branch_example.verbosity), NULL, 0.0, 0, "Larger values result in more output during running."},
    {"DEFAULT_TO_ELSE", "", IS_LONG, 0, (long)((char *)&branch_example.defaultToElse), NULL, 0.0, 0, "If non-zero, defaults to ELSE_TO when performing tracking for closed orbit, twiss_output, etc."},
    {"BRANCH_TO", "", IS_STRING, 0, (long)((char *)&branch_example.branchTo), NULL, 0.0, 0, "Optional name of element to which to jump when counter is non-positive."},
    {"ELSE_TO", "", IS_STRING, 0, (long)((char *)&branch_example.elseTo), NULL, 0.0, 0, "Optional name of element to which to jump when counter is positive."},
    } ;

QFRING qfring_example;
/* quadrupole fringe-field physical parameters */
PARAMETER qfring_param[N_QFRING_PARAMS]={
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&qfring_example.length), NULL, 0.0, 0, "length"},
    {"K1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&qfring_example.k1), NULL, 0.0, 0, "peak geometric strength"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&qfring_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&qfring_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&qfring_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&qfring_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&qfring_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"DIRECTION", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&qfring_example.direction), NULL, 0.0, 0, "1=entrance, -1=exit"},
    {"ORDER", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&qfring_example.order), NULL, 0.0, 0, "matrix order"}
    };

SCRAPER scraper_example;
/* scraper physical parameters */
PARAMETER scraper_param[N_SCRAPER_PARAMS]={
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&scraper_example.length), NULL, 0.0, 0, "length"},
    {"XO", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&scraper_example.Xo), NULL, 0.0, 0, "radiation length"},
    {"ENERGY_DECAY", "", IS_LONG, 0, (long)((char *)&scraper_example.energyDecay), NULL, 0.0, 0, "If nonzero, then particles will lose energy due to material using a simple exponential model."},
    {"ENERGY_STRAGGLE", "", IS_LONG, 0, (long)((char *)&scraper_example.energyStraggle), NULL, 0.0, 0, "Use simple-minded energy straggling model coupled with ENERGY_DECAY=1?"},
    {"NUCLEAR_BREMSSTRAHLUNG", "", IS_LONG, 0, (long)((char *)&scraper_example.nuclearBremsstrahlung), NULL, 0.0, 0, "Model energy loss to nuclear bremsstrahlung? If enabled, set ENERGY_DECAY=0 to disable simpler model."},
    {"ELECTRON_RECOIL", "", IS_LONG, 0, (long)((char *)&scraper_example.electronRecoil), NULL, 0.0, 0, "If non-zero, electron recoil during Coulomb scattering is included (results in energy change)."},
    {"Z", "", IS_LONG, 0, (long)((char *)&scraper_example.Z), NULL, 0.0, 0, "Atomic number"},
    {"A", "AMU", IS_DOUBLE, 0, (long)((char *)&scraper_example.A), NULL, 0.0, 0, "Atomic mass"},
    {"RHO", "KG/M^3", IS_DOUBLE, 0, (long)((char *)&scraper_example.rho), NULL, 0.0, 0, "Density"},       
    {"PLIMIT", "", IS_DOUBLE, 0, (long)((char *)&scraper_example.pLimit), NULL, 0.05, 0, "Probability cutoff for each slice"},
    {"POSITION", "M", IS_DOUBLE, 0, (long)((char *)&scraper_example.position), NULL, 0.0, 0, "position of edge"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&scraper_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&scraper_example.dy), NULL, 0.0, 0, "misalignment"},
    {"INSERT_FROM", "", IS_STRING, 0, (long)((char *)&scraper_example.insert_from), NULL, 0.0, 0, "direction from which inserted (+x, -x, x, +y, -y, y)"},
    {"DIRECTION", "", IS_LONG, 0, (long)((char *)&scraper_example.oldDirection), NULL, 0.0, -1, "Deprecated. use INSERT_FROM."},
    };

CENTER center_example;
/* beam centering physical parameters */
PARAMETER center_param[N_CENTER_PARAMS]={
    {"X" , "", IS_LONG, 0, (long)((char *)&center_example.doCoord[0]), NULL, 0.0, 1, "center x coordinates?"},
    {"XP", "", IS_LONG, 0, (long)((char *)&center_example.doCoord[1]), NULL, 0.0, 1, "center x' coordinates?"},
    {"Y" , "", IS_LONG, 0, (long)((char *)&center_example.doCoord[2]), NULL, 0.0, 1, "center y coordinates?" },
    {"YP", "", IS_LONG, 0, (long)((char *)&center_example.doCoord[3]), NULL, 0.0, 1, "center y' coordinates?"},
    {"S" , "", IS_LONG, 0, (long)((char *)&center_example.doCoord[4]), NULL, 0.0, 0, "center s coordinates?" },
    {"DELTA", "", IS_LONG, 0, (long)((char *)&center_example.doCoord[5]), NULL, 0.0, 0, "center delta coordinates?"},
    {"T", "", IS_LONG, 0, (long)((char *)&center_example.doCoord[6]), NULL, 0.0, 0, "center t coordinates?"},
    {"ONCE_ONLY", "", IS_LONG, 0, (long)((char *)&center_example.onceOnly), NULL, 0.0, 0, "compute centering offsets for first beam only, apply to all?"},
    {"ON_PASS", "", IS_LONG, 0, (long)((char *)&center_example.onPass), NULL, 0.0, -1, "If nonnegative, do centering on the nth pass only."},
    };

KICKER kicker_example;
/* kicker physical parameters */
PARAMETER kicker_param[N_KICKER_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&kicker_example.length), NULL, 0.0, 0, "length"},
    {"ANGLE", "RAD", IS_DOUBLE, 0, (long)((char *)&kicker_example.angle), NULL, 0.0, 0, "kick angle"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&kicker_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kicker_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kicker_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kicker_example.dz), NULL, 0.0, 0, "misalignment"},
    {"B2", "1/M^2", IS_DOUBLE, 0, (long)((char *)&kicker_example.b2), NULL, 0.0, 0, "Sextupole term: By=Bo*(1+b2*x^2)"},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&kicker_example.time_offset), NULL, 0.0, 0, "time offset of waveform"},
    {"PERIODIC", "", IS_LONG, 0, (long)((char *)&kicker_example.periodic), NULL, 0.0, 0, "is waveform periodic?"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&kicker_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"FIRE_ON_PASS", "", IS_LONG, 0, (long)((char *)&kicker_example.fire_on_pass), NULL, 0.0, 0, "pass number to fire on"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&kicker_example.n_kicks), NULL, 0.0, 0, "Number of kicks to use for simulation. 0 uses an exact result but ignores b2."},
    {"WAVEFORM", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&kicker_example.waveform), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving kick factor vs time"},
    {"DEFLECTION_MAP", "", IS_STRING, IS_STRING, (long)((char *)&kicker_example.deflectionMap), NULL, 0.0, 0, "optional filename giving the spatial variation of the deflection"},
    } ;

MKICKER mkicker_example;
/* multipole kicker physical parameters */
PARAMETER mkicker_param[N_MKICKER_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&mkicker_example.length), NULL, 0.0, 0, "length"},
    {"STRENGTH", "", IS_DOUBLE, 0, (long)((char *)&mkicker_example.strength), NULL, 0.0, 0, "geometric strength in 1/m^order"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&mkicker_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mkicker_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mkicker_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mkicker_example.dz), NULL, 0.0, 0, "misalignment"},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&mkicker_example.time_offset), NULL, 0.0, 0, "time offset of waveform"},
    {"ORDER", "", IS_LONG, 0, (long)((char *)&mkicker_example.order), NULL, 0.0, 0, "multipole order, where 1 is quadrupole, 2 is sextupole, etc."},
    {"PERIODIC", "", IS_LONG, 0, (long)((char *)&mkicker_example.periodic), NULL, 0.0, 0, "is waveform periodic?"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&mkicker_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"FIRE_ON_PASS", "", IS_LONG, 0, (long)((char *)&mkicker_example.fire_on_pass), NULL, 0.0, 0, "pass number to fire on"},
    {"N_KICKS", "", IS_LONG, 1, (long)((char *)&mkicker_example.n_kicks), NULL, 0.0, 0, "Number of kicks to use for simulation."},
    {"WAVEFORM", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&mkicker_example.waveform), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving kick factor vs time"},
    } ;

KSEXT ksext_example;
/* kick sextupole physical parameters */
PARAMETER ksext_param[N_KSEXT_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&ksext_example.length), NULL, 0.0, 0, "length"},
    {"K2", "1/M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.k2), NULL, 0.0, 0, "geometric strength"},
    {"K1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.k1), NULL, 0.0, 0, "geometric quadrupole strength error. See notes below!"},
    {"J1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.j1), NULL, 0.0, 0, "geometric skew quadrupole strength error. See notes below!"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"BORE", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.bore), NULL, 0.0, 0, "bore radius"},
    {"B", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.B), NULL, 0.0, 0, "field at pole tip (used if bore nonzero)"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&ksext_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks (rounded up to next multipole of 4 if INTEGRATION_ORDER=4)"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"HKICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&ksext_example.xkick), NULL, 0.0, 0, "horizontal correction kick"},
    {"VKICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&ksext_example.ykick), NULL, 0.0, 0, "vertical correction kick"},
    {"HCALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.xKickCalibration), NULL, 1.0, 0, "calibration factor for horizontal correction kick"},
    {"VCALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.yKickCalibration), NULL, 1.0, 0, "calibration factor for vertical correction kick"},
    {"HSTEERING", "", IS_SHORT, 0, (long)((char *)&ksext_example.xSteering), NULL, 0.0, 0, "use for horizontal correction?"},
    {"VSTEERING", "", IS_SHORT, 0, (long)((char *)&ksext_example.ySteering), NULL, 0.0, 0, "use for vertical correction?"},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&ksext_example.synch_rad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"SYSTEMATIC_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ksext_example.systematic_multipoles), NULL, 0.0, 0, "input file for systematic multipoles"},
    {"EDGE_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ksext_example.edge_multipoles), NULL, 0.0, 0, "input file for systematic edge multipoles"},
    {"RANDOM_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ksext_example.random_multipoles), NULL, 0.0, 0, "input file for random multipoles"},
    {"STEERING_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ksext_example.steering_multipoles), NULL, 0.0, 0, "input file for multipole content of steering kicks"},
    {"SYSTEMATIC_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.systematicMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply systematic and edge multipoles"},
    {"RANDOM_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.randomMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply random multipoles"},
    {"STEERING_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksext_example.steeringMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply steering multipoles"},
    {"INTEGRATION_ORDER", "", IS_SHORT, 0, (long)((char *)&ksext_example.integration_order), NULL, 0.0, 4, "integration order (2 or 4)"},
    {"SQRT_ORDER", "", IS_SHORT, 0, (long)((char *)&ksext_example.sqrtOrder), NULL, 0.0, 0, "Ignored, kept for backward compatibility only."},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&ksext_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"ISR1PART", "", IS_SHORT, 0, (long)((char *)&ksext_example.isr1Particle), NULL, 0.0, 1, "Include ISR for single-particle beam only if ISR=1 and ISR1PART=1"},
    {"EXPAND_HAMILTONIAN", "", IS_SHORT, 0, (long)((char *)&ksext_example.expandHamiltonian), NULL, 0.0, 0, "If 1, Hamiltonian is expanded to leading order."},
    };

KOCT koct_example;
/* kick octupole physical parameters */
PARAMETER koct_param[N_KOCT_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&koct_example.length), NULL, 0.0, 0, "length"},
    {"K3", "1/M$a4$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&koct_example.k3), NULL, 0.0, 0, "geometric strength"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&koct_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"BORE", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&koct_example.bore), NULL, 0.0, 0, "bore radius"},
    {"B", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&koct_example.B), NULL, 0.0, 0, "field at pole tip (used if bore nonzero)"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&koct_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&koct_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&koct_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&koct_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&koct_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks (rounded up to next multipole of 4 if INTEGRATION_ORDER=4)"},
    {"SYSTEMATIC_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&koct_example.systematic_multipoles), NULL, 0.0, 0, "input file for systematic multipoles"},
    {"RANDOM_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&koct_example.random_multipoles), NULL, 0.0, 0, "input file for random multipoles"},
    {"INTEGRATION_ORDER", "", IS_SHORT, 0, (long)((char *)&koct_example.integration_order), NULL, 0.0, 4, "integration order (2 or 4)"},
    {"SQRT_ORDER", "", IS_SHORT, 0, (long)((char *)&koct_example.sqrtOrder), NULL, 0.0, 0, "Ignored, kept for backward compatibility only."},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&koct_example.synch_rad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&koct_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"ISR1PART", "", IS_SHORT, 0, (long)((char *)&koct_example.isr1Particle), NULL, 0.0, 1, "Include ISR for single-particle beam only if ISR=1 and ISR1PART=1"},
    {"EXPAND_HAMILTONIAN", "", IS_SHORT, 0, (long)((char *)&koct_example.expandHamiltonian), NULL, 0.0, 0, "If 1, Hamiltonian is expanded to leading order."},
    };

KSBEND ksbend_example;
/* symplectic sector bending magnet physical parameters */
PARAMETER ksbend_param[N_KSBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.angle), NULL, 0.0, 0, "bend angle"},
    {"K1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.k1), NULL, 0.0, 0, "geometric quadrupole strength"},
    {"K2", "1/M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.k2), NULL, 0.0, 0, "geometric sextupole strength"},
    {"K3", "1/M$a4$n", IS_DOUBLE, 0, (long)((char *)&ksbend_example.k3), NULL, 0.0, 0, "geometric octupole strength"},
    {"K4", "1/M$a5$n", IS_DOUBLE, 0, (long)((char *)&ksbend_example.k4), NULL, 0.0, 0, "geometric decapole strength"},
    {"E1", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.e[0]), NULL, 0.0, 0, "entrance edge angle"},
    {"E2", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.e[1]), NULL, 0.0, 0, "exit edge angle"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"H1", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.h[0]), NULL, 0.0, 0, "entrance pole-face curvature"},
    {"H2", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.h[1]), NULL, 0.0, 0, "exit pole-face curvature"},
    {"HGAP", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.hgap), NULL, 0.0, 0, "half-gap between poles"},
    {"FINT", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.fint), NULL, DEFAULT_FINT, 0, "edge-field integral"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"ETILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.etilt), NULL, 0.0, 0, "error rotation about incoming longitudinal axis"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&ksbend_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks"},
    {"NONLINEAR", "", IS_LONG, 0, (long)((char *)&ksbend_example.nonlinear), NULL, 0.0, 1, "include nonlinear field components?"},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&ksbend_example.synch_rad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"EDGE1_EFFECTS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.edge_effects[0]), NULL, 0.0, 1, "include entrance edge effects?"},
    {"EDGE2_EFFECTS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.edge_effects[1]), NULL, 0.0, 1, "include exit edge effects?"},
    {"EDGE_ORDER", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.edge_order), NULL, 0.0, 1, "edge matrix order"},
    {"PARAXIAL", "", IS_LONG, 0, (long)((char *)&ksbend_example.paraxial), NULL, 0.0, 0, "use paraxial approximation?"},
    {"TRANSPORT", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&ksbend_example.TRANSPORT), NULL, 0.0, 0, "use (incorrect) TRANSPORT equations for T436 of edge?"},
    {"METHOD", "", IS_STRING, 0, (long)((char *)&ksbend_example.method), "modified-midpoint", 0.0, 0, "integration method (modified-midpoint, leap-frog"}
    };

KQUAD kquad_example;
/* kick quadrupole physical parameters */
PARAMETER kquad_param[N_KQUAD_PARAMS]={
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&kquad_example.length), NULL, 0.0, 0, "length"},
    {"K1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.k1), NULL, 0.0, 0, "geometric strength"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"BORE", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.bore), NULL, 0.0, 0, "bore radius"},
    {"B", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.B), NULL, 0.0, 0, "pole tip field (used if bore nonzero)"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"N_KICKS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks (rounded up to next multipole of 4 if INTEGRATION_ORDER=4)"},
    {"HKICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&kquad_example.xkick), NULL, 0.0, 0, "horizontal correction kick"},
    {"VKICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&kquad_example.ykick), NULL, 0.0, 0, "vertical correction kick"},
    {"HCALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.xKickCalibration), NULL, 1.0, 0, "calibration factor for horizontal correction kick"},
    {"VCALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.yKickCalibration), NULL, 1.0, 0, "calibration factor for vertical correction kick"},
    {"HSTEERING", "", IS_SHORT, 0, (long)((char *)&kquad_example.xSteering), NULL, 0.0, 0, "use for horizontal correction?"},
    {"VSTEERING", "", IS_SHORT, 0, (long)((char *)&kquad_example.ySteering), NULL, 0.0, 0, "use for vertical correction?"},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&kquad_example.synch_rad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"SYSTEMATIC_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&kquad_example.systematic_multipoles), NULL, 0.0, 0, "input file for systematic multipoles"},
    {"EDGE_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&kquad_example.edge_multipoles), NULL, 0.0, 0, "input file for systematic edge multipoles"},
    {"RANDOM_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&kquad_example.random_multipoles), NULL, 0.0, 0, "input file for random multipoles"},
    {"STEERING_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&kquad_example.steering_multipoles), NULL, 0.0, 0, "input file for multipole content of steering kicks"},
    {"SYSTEMATIC_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.systematicMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply systematic and edge multipoles"},
    {"RANDOM_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.randomMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply random multipoles"},
    {"STEERING_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.steeringMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply steering multipoles"},
    {"INTEGRATION_ORDER", "", IS_SHORT, 0, (long)((char *)&kquad_example.integration_order), NULL, 0.0, 4, "integration order (2 or 4)"},
    {"SQRT_ORDER", "", IS_SHORT, 0, (long)((char *)&kquad_example.sqrtOrder), NULL, 0.0, 0, "Ignored, kept for backward compatibility only."},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&kquad_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"ISR1PART", "", IS_SHORT, 0, (long)((char *)&kquad_example.isr1Particle), NULL, 0.0, 1, "Include ISR for single-particle beam only if ISR=1 and ISR1PART=1"},
    {"EDGE1_EFFECTS", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.edge1_effects), NULL, 0.0, 0, "include entrance edge effects?"},
    {"EDGE2_EFFECTS", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.edge2_effects), NULL, 0.0, 0, "include exit edge effects?"},
    {"LEFFECTIVE", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.lEffective), NULL, 0.0, 0, "Effective length. Ignored if non-positive."},
    {"I0P", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.fringeIntP[0]), NULL, 0.0, 0, "i0+ fringe integral"},
    {"I1P", "M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.fringeIntP[1]), NULL, 0.0, 0, "i1+ fringe integral"},
    {"I2P", "M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.fringeIntP[2]), NULL, 0.0, 0, "i2+ fringe integral"},
    {"I3P", "M$a4$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.fringeIntP[3]), NULL, 0.0, 0, "i3+ fringe integral"},
    {"LAMBDA2P", "M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.fringeIntP[4]), NULL, 0.0, 0, "lambda2+ fringe integral"},
    {"I0M", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.fringeIntM[0]), NULL, 0.0, 0, "i0- fringe integral"},
    {"I1M", "M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.fringeIntM[1]), NULL, 0.0, 0, "i1- fringe integral"},
    {"I2M", "M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.fringeIntM[2]), NULL, 0.0, 0, "i2- fringe integral"},
    {"I3M", "M$a4$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.fringeIntM[3]), NULL, 0.0, 0, "i3- fringe integral"},
    {"LAMBDA2M", "M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.fringeIntM[4]), NULL, 0.0, 0, "lambda2- fringe integral"},
    {"EDGE1_LINEAR", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.edge1Linear), NULL, 0.0, 1, "Use to selectively turn off linear part if EDGE1_EFFECTS nonzero."},
    {"EDGE2_LINEAR", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.edge2Linear), NULL, 0.0, 1, "Use to selectively turn off linear part if EDGE2_EFFECTS nonzero."},
    {"EDGE1_NONLINEAR_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.edge1NonlinearFactor), NULL, 1.0, 0, "Use to selectively scale nonlinear entrance edge effects if EDGE1_EFFECTS>1"},
    {"EDGE2_NONLINEAR_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.edge2NonlinearFactor), NULL, 1.0, 0, "Use to selectively scale nonlinear exit edge effects if EDGE2_EFFECTS>1"},
    {"RADIAL", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&kquad_example.radial), NULL, 0.0, 0, "If non-zero, converts the quadrupole into a radially-focusing lens"},
    {"EXPAND_HAMILTONIAN", "", IS_SHORT, 0, (long)((char *)&kquad_example.expandHamiltonian), NULL, 0.0, 0, "If 1, Hamiltonian is expanded to leading order."},
    {"TRACKING_MATRIX", "", IS_SHORT, 0, (long)((char *)&kquad_example.trackingBasedMatrix), NULL, 0.0, 0, "If nonzero, gives order of tracking-based matrix up to third order to be used for twiss parameters etc.  If zero, 2nd-order analytical matrix is used."},
    };

MAGNIFY magnify_example;
/* magnifier physical parameters */
PARAMETER magnify_param[N_MAGNIFY_PARAMS] = {
    {"MX", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&magnify_example.mx), NULL, 1.0, 0, "factor for x coordinates"},
    {"MXP", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&magnify_example.mxp), NULL, 1.0, 0, "factor for x' coordinates"},
    {"MY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&magnify_example.my), NULL, 1.0, 0, "factor for y coordinates"},
    {"MYP", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&magnify_example.myp), NULL, 1.0, 0, "factor for y' coordinates"},
    {"MS", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&magnify_example.ms), NULL, 1.0, 0, "factor for s coordinates"},
    {"MDP", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&magnify_example.mdp), NULL, 1.0, 0, "factor for (p-pCentral)/pCentral"},
    } ;
    
SAMPLE sample_example;
/* sample physical parameters */
PARAMETER sample_param[N_SAMPLE_PARAMS] = {
    {"FRACTION", "", IS_DOUBLE, 0, (long)((char *)&sample_example.fraction), NULL, 1.0, 0, "fraction to keep"},
    {"INTERVAL", "", IS_LONG, 0, (long)((char *)&sample_example.interval), NULL, 1.0, 1, "interval between sampled particles"},
    } ;
    
HVCOR hvcor_example;
/* horizontal/vertical corrector physical parameters */
PARAMETER hvcor_param[N_HVCOR_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hvcor_example.length), NULL, 0.0, 0, "length"},
    {"HKICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hvcor_example.xkick), NULL, 0.0, 0, "x kick angle"},
    {"VKICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hvcor_example.ykick), NULL, 0.0, 0, "y kick angle"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hvcor_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"B2", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hvcor_example.b2), NULL, 0.0, 0, "normalized sextupole strength (e.g., kick = KICK*(1+B2*x^2))"},
    {"HCALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hvcor_example.xcalibration), NULL, 1.0, 0, "factor applied to obtain x kick"},
    {"VCALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hvcor_example.ycalibration), NULL, 1.0, 0, "factor applied to obtain y kick"},
    {"EDGE_EFFECTS", "", IS_LONG, 0, (long)((char *)&hvcor_example.edge_effects), NULL, 0.0, 0, "include edge effects?"},
    {"ORDER", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&hvcor_example.order), NULL, 0.0, 0, "matrix order"},
    {"STEERING", "", IS_LONG, 0, (long)((char *)&hvcor_example.steering), NULL, 0.0, 1, "use for steering?"},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&hvcor_example.synchRad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"ISR", "", IS_LONG, 0, (long)((char *)&hvcor_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"LERAD", "", IS_DOUBLE, 0, (long)((char *)&hvcor_example.lEffRad), NULL, 0.0, 0, "if L=0, use this length for radiation computations"},
    };

SCATTER scatter_example;
/* scatter physical parameters */
PARAMETER scatter_param[N_SCATTER_PARAMS] = {
    {"X", "M", IS_DOUBLE, 0, (long)((char*)&scatter_example.x), NULL, 0.0, 0, "rms scattering level for x"},
    {"XP", "", IS_DOUBLE, 0, (long)((char*)&scatter_example.xp), NULL, 0.0, 0, "rms scattering level for x'"},
    {"Y", "M", IS_DOUBLE, 0, (long)((char*)&scatter_example.y), NULL, 0.0, 0, "rms scattering level for y"},
    {"YP", "", IS_DOUBLE, 0, (long)((char*)&scatter_example.yp), NULL, 0.0, 0, "rms scattering level for y'"},
    {"DP", "", IS_DOUBLE, 0, (long)((char*)&scatter_example.dp), NULL, 0.0, 0, "rms scattering level for (p-pCentral)/pCentral"},
    {"PROBABILITY", "", IS_DOUBLE, 0, (long)((char*)&scatter_example.probability), NULL, 1.0, 0, "Probability that any particle will be selected for scattering."},
    {"STARTONPASS", "", IS_LONG, 0, (long)((char*)&scatter_example.startOnPass), NULL, 0.0, 0, "Pass number to start on."},
    {"ENDONPASS", "", IS_LONG, 0, (long)((char*)&scatter_example.endOnPass), NULL, 0.0, -1, "Pass number to end on (inclusive).  Ignored if negative."},
    } ;

DSCATTER dscatter_example;
/* dscatter physical parameters */
PARAMETER dscatter_param[N_DSCATTER_PARAMS] = {
    {"PLANE", "", IS_STRING, 0, (long)((char*)&dscatter_example.plane), NULL, 0.0, 0, "Plane to scatter: xp, yp, dp (dp is deltaP/P)"},
    {"FILENAME", "", IS_STRING, 0, (long)((char*)&dscatter_example.fileName), NULL, 0.0, 0, "Name of SDDS file containing distribution function."},
    {"VALUENAME", "", IS_STRING, 0, (long)((char*)&dscatter_example.valueName), NULL, 0.0, 0, "Name of column containing the independent variable for the distribution function data."},
    {"CDFNAME", "", IS_STRING, 0, (long)((char*)&dscatter_example.cdfName), NULL, 0.0, 0, "Name of column containing the cumulative distribution function data."},
    {"PDFNAME", "", IS_STRING, 0, (long)((char*)&dscatter_example.pdfName), NULL, 0.0, 0, "Name of column containing the probability distribution function data."},
    {"ONCEPERPARTICLE", "", IS_LONG, 0, (long)((char*)&dscatter_example.oncePerParticle), NULL, 0.0, 0, "If nonzero, each particle can only get scattered once by this element."},
    {"FACTOR", "", IS_DOUBLE, 0, (long)((char*)&dscatter_example.factor), NULL, 1.0, 0, "Factor by which to multiply the independent variable values."},
    {"PROBABILITY", "", IS_DOUBLE, 0, (long)((char*)&dscatter_example.probability), NULL, 1.0, 0, "Probability that any particle will be selected for scattering."},
    {"GROUPID", "", IS_LONG, 0, (long)((char*)&dscatter_example.group), NULL, 0.0, -1, "Group ID number (nonnegative integer) for linking once-per-particle behavior of multiple elements."},
    {"RANDOMSIGN", "", IS_LONG, 0, (long)((char*)&dscatter_example.randomSign), NULL, 0.0, 0, "If non-zero, then the scatter is given a random sign.  Useful if distribution data is one-sided."},
    {"LIMITPERPASS", "", IS_LONG, 0, (long)((char*)&dscatter_example.limitPerPass), NULL, 0.0, -1, "Maximum number of particles that will be scattered on each pass."},
    {"LIMITTOTAL", "", IS_LONG, 0, (long)((char*)&dscatter_example.limitTotal), NULL, 0.0, -1, "Maximum number of particles that will be scatter for each step."},
    {"STARTONPASS", "", IS_LONG, 0, (long)((char*)&dscatter_example.startOnPass), NULL, 0.0, 0, "Pass number to start on."},
    {"ENDONPASS", "", IS_LONG, 0, (long)((char*)&dscatter_example.endOnPass), NULL, 0.0, -1, "Pass number to end on (inclusive).  Ignored if negative."},
    } ;
    
NIBEND nibend_example;
/* integrated bending magnet physical parameters */
PARAMETER nibend_param[N_NIBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.angle), NULL, 0.0, 0, "bending angle"},
    {"E1", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.e[0]), NULL, 0.0, 0, "entrance edge angle"},
    {"E2", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.e[1]), NULL, 0.0, 0, "exit edge angle"},
    {"TILT", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FINT", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.fint), NULL, DEFAULT_FINT, 0, "edge-field integral"},
    {"HGAP", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.hgap), NULL, 0.0, 0, "half-gap between poles"},
    {"FP1", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.fp1), NULL, 10.0, 0, "fringe parameter (tanh model)"},
    {"FP2", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.fp2), NULL, 0.0, 0, "not used"},
    {"FP3", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.fp3), NULL, 0.0, 0, "not used"},
    {"FP4", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.fp4), NULL, 0.0, 0, "not used"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"ETILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nibend_example.etilt), NULL, 0.0, 0, "error rotation about incoming longitudinal axis"},
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&nibend_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy (for nonadaptive integration, used as the step-size)"},
    {"MODEL", "", IS_STRING, 0, (long)((char *)&nibend_example.model), DEFAULT_NIBEND_TYPE, 0.0, 0, "fringe model (hard-edge, linear, cubic-spline, tanh, quintic, enge1, enge3, enge5)"},
    {"METHOD", "", IS_STRING, 0, (long)((char *)&nibend_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, modified-midpoint, two-pass modified-midpoint, leap-frog, non-adaptive runge-kutta)"},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&nibend_example.synch_rad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"ADJUST_BOUNDARY", "", IS_LONG, 0, (long)((char *)&nibend_example.adjustBoundary), NULL, 0.0, 1, "adjust fringe boundary position to make symmetric trajectory? (Not done if ADJUST_FIELD is nonzero.)"},
    {"ADJUST_FIELD", "", IS_LONG, 0, (long)((char *)&nibend_example.adjustField), NULL, 0.0, 0, "adjust central field strength to make symmetric trajectory?"},
    {"FUDGE_PATH_LENGTH", "", IS_LONG, 0, (long)((char *)&nibend_example.fudgePathLength), NULL, 0.0, 1, "fudge central path length to force it to equal the nominal length L?"},
    {"FRINGE_POSITION", "", IS_LONG, 0, (long)((char *)&nibend_example.fringePosition), NULL, 0.0, 0, "0=fringe centered on reference plane, -1=fringe inside, 1=fringe outside."},
    };

KPOLY kpoly_example;
/* kick-polynomial physical parameters */
PARAMETER kpoly_param[N_KPOLY_PARAMS] = {
    {"COEFFICIENT", "M$A-ORDER$N", IS_DOUBLE, 0, (long)((char *)&kpoly_example.coefficient), NULL, 0.0, 0, "coefficient of polynomial"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&kpoly_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&kpoly_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&kpoly_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 0, (long)((char *)&kpoly_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FACTOR", "", IS_DOUBLE, 0, (long)((char *)&kpoly_example.factor), NULL, 1.0, 0, "additional factor to apply"},
    {"ORDER", "", IS_LONG, 0, (long)((char *)&kpoly_example.order), NULL, 0.0, 0, "order of polynomial"},
    {"PLANE", "", IS_STRING, 0, (long)((char *)&kpoly_example.plane), "x", 0.0, 0, "plane to kick (x, y)"}
    };

HKPOLY hkpoly_example;
/* kick-polynomial physical parameters */
PARAMETER hkpoly_param[N_HKPOLY_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.length), NULL, 0.0, 0, "length for geometry only, ignored in tracking"},
{"K00", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[0][0]), NULL, 0.0, 0, "Coefficient of polynomial for kicks---ignored"},
{"K01", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[0][1]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K02", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[0][2]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K03", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[0][3]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K04", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[0][4]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K05", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[0][5]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K06", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[0][6]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K10", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[1][0]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K11", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[1][1]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K12", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[1][2]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K13", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[1][3]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K14", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[1][4]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K15", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[1][5]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K16", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[1][6]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K20", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[2][0]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K21", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[2][1]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K22", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[2][2]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K23", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[2][3]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K24", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[2][4]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K25", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[2][5]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K26", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[2][6]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K30", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[3][0]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K31", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[3][1]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K32", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[3][2]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K33", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[3][3]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K34", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[3][4]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K35", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[3][5]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K36", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[3][6]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K40", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[4][0]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K41", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[4][1]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K42", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[4][2]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K43", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[4][3]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K44", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.K[4][4]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K45", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[4][5]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K46", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[4][6]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K50", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[5][0]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K51", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[5][1]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K52", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[5][2]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K53", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[5][3]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K54", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[5][4]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K55", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[5][5]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K56", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[5][6]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K60", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[6][0]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K61", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[6][1]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K62", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[6][2]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K63", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[6][3]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K64", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[6][4]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K65", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[6][5]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"K66", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.K[6][6]), NULL, 0.0, 0, "Coefficient of polynomial for kicks"},
{"D00", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[0][0]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift---ignored"},
{"D01", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[0][1]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D02", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[0][2]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D03", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[0][3]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D04", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[0][4]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D05", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[0][5]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D06", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[0][6]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D10", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[1][0]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D11", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[1][1]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D12", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[1][2]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D13", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[1][3]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D14", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[1][4]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D15", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[1][5]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D16", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[1][6]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D20", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[2][0]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D21", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[2][1]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D22", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[2][2]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D23", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[2][3]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D24", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[2][4]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D25", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[2][5]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D26", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[2][6]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D30", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[3][0]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D31", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[3][1]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D32", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[3][2]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D33", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[3][3]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D34", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[3][4]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D35", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[3][5]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D36", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[3][6]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D40", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[4][0]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D41", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[4][1]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D42", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[4][2]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D43", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[4][3]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D44", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.D[4][4]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D45", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[4][5]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D46", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[4][6]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D50", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[5][0]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D51", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[5][1]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D52", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[5][2]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D53", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[5][3]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D54", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[5][4]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D55", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[5][5]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D56", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[5][6]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D60", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[6][0]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D61", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[6][1]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D62", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[6][2]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D63", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[6][3]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D64", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[6][4]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D65", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[6][5]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"D66", "", IS_DOUBLE, 0, (long)((char *)&hkpoly_example.D[6][6]), NULL, 0.0, 0, "Coefficient of polynomial for generalized drift"},
{"E000", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][0][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E001", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][0][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E002", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][0][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E003", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][0][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E004", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][0][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E005", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][0][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E006", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][0][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E010", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][1][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E011", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][1][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E012", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][1][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E013", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][1][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E014", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][1][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E015", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][1][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E016", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][1][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E020", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][2][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E021", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][2][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E022", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][2][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E023", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][2][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E024", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][2][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E025", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][2][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E026", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][2][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E030", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][3][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E031", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][3][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E032", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][3][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E033", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][3][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E034", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][3][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E035", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][3][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E036", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][3][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E040", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[0][4][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E041", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][4][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E042", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][4][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E043", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][4][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E044", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][4][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E045", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][4][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E046", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][4][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E050", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][5][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E051", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][5][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E052", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][5][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E053", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][5][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E054", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][5][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E055", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][5][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E056", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][5][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E060", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][6][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E061", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][6][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E062", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][6][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E063", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][6][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E064", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][6][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E065", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][6][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E066", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[0][6][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E100", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[1][0][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E101", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[1][0][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E102", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[1][0][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E103", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[1][0][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E104", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][0][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E105", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][0][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E106", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][0][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E110", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[1][1][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E111", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[1][1][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E112", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[1][1][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E113", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][1][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E114", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][1][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E115", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][1][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E116", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][1][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E120", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[1][2][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E121", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[1][2][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E122", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][2][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E123", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][2][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E124", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][2][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E125", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][2][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E126", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][2][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E130", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[1][3][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E131", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][3][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E132", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][3][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E133", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][3][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E134", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][3][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E135", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][3][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E136", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][3][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E140", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][4][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E141", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][4][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E142", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][4][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E143", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][4][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E144", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][4][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E145", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][4][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E146", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][4][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E150", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][5][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E151", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][5][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E152", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][5][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E153", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][5][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E154", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][5][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E155", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][5][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E156", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][5][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E160", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][6][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E161", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][6][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E162", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][6][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E163", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][6][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E164", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][6][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E165", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][6][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E166", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[1][6][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E200", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[2][0][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E201", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[2][0][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E202", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[2][0][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E203", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][0][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E204", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][0][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E205", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][0][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E206", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][0][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E210", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[2][1][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E211", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[2][1][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E212", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][1][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E213", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][1][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E214", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][1][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E215", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][1][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E216", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][1][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E220", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[2][2][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E221", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][2][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E222", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][2][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E223", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][2][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E224", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][2][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E225", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][2][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E226", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][2][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E230", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][3][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E231", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][3][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E232", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][3][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E233", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][3][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E234", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][3][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E235", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][3][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E236", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][3][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E240", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][4][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E241", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][4][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E242", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][4][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E243", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][4][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E244", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][4][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E245", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][4][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E246", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][4][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E250", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][5][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E251", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][5][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E252", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][5][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E253", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][5][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E254", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][5][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E255", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][5][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E256", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][5][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E260", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][6][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E261", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][6][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E262", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][6][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E263", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][6][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E264", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][6][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E265", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][6][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E266", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[2][6][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E300", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[3][0][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E301", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[3][0][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E302", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][0][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E303", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][0][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E304", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][0][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E305", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][0][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E306", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][0][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E310", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[3][1][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E311", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][1][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E312", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][1][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E313", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][1][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E314", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][1][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E315", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][1][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E316", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][1][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E320", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][2][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E321", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][2][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E322", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][2][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E323", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][2][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E324", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][2][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E325", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][2][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E326", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][2][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E330", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][3][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E331", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][3][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E332", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][3][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E333", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][3][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E334", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][3][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E335", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][3][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E336", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][3][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E340", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][4][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E341", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][4][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E342", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][4][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E343", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][4][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E344", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][4][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E345", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][4][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E346", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][4][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E350", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][5][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E351", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][5][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E352", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][5][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E353", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][5][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E354", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][5][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E355", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][5][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E356", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][5][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E360", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][6][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E361", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][6][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E362", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][6][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E363", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][6][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E364", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][6][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E365", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][6][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E366", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[3][6][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E400", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.E[4][0][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E401", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][0][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E402", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][0][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E403", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][0][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E404", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][0][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E405", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][0][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E406", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][0][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E410", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][1][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E411", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][1][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E412", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][1][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E413", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][1][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E414", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][1][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E415", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][1][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E416", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][1][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E420", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][2][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E421", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][2][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E422", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][2][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E423", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][2][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E424", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][2][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E425", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][2][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E426", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][2][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E430", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][3][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E431", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][3][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E432", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][3][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E433", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][3][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E434", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][3][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E435", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][3][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E436", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][3][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E440", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][4][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E441", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][4][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E442", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][4][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E443", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][4][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E444", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][4][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E445", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][4][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E446", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][4][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E450", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][5][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E451", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][5][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E452", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][5][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E453", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][5][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E454", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][5][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E455", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][5][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E456", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][5][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E460", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][6][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E461", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][6][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E462", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][6][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E463", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][6][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E464", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][6][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E465", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][6][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E466", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[4][6][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E500", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][0][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E501", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][0][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E502", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][0][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E503", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][0][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E504", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][0][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E505", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][0][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E506", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][0][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E510", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][1][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E511", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][1][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E512", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][1][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E513", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][1][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E514", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][1][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E515", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][1][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E516", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][1][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E520", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][2][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E521", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][2][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E522", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][2][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E523", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][2][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E524", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][2][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E525", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][2][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E526", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][2][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E530", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][3][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E531", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][3][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E532", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][3][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E533", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][3][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E534", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][3][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E535", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][3][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E536", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][3][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E540", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][4][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E541", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][4][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E542", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][4][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E543", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][4][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E544", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][4][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E545", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][4][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E546", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][4][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E550", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][5][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E551", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][5][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E552", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][5][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E553", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][5][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E554", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][5][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E555", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][5][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E556", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][5][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E560", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][6][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E561", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][6][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E562", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][6][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E563", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][6][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E564", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][6][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E565", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][6][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E566", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[5][6][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E600", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][0][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E601", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][0][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E602", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][0][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E603", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][0][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E604", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][0][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E605", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][0][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E606", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][0][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E610", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][1][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E611", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][1][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E612", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][1][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E613", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][1][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E614", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][1][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E615", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][1][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E616", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][1][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E620", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][2][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E621", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][2][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E622", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][2][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E623", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][2][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E624", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][2][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E625", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][2][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E626", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][2][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E630", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][3][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E631", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][3][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E632", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][3][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E633", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][3][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E634", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][3][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E635", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][3][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E636", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][3][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E640", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][4][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E641", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][4][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E642", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][4][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E643", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][4][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E644", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][4][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E645", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][4][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E646", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][4][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E650", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][5][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E651", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][5][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E652", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][5][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E653", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][5][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E654", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][5][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E655", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][5][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E656", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][5][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E660", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][6][0]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E661", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][6][1]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E662", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][6][2]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E663", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][6][3]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E664", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][6][4]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E665", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][6][5]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},
{"E666", "", IS_DOUBLE,                    0, (long)((char *)&hkpoly_example.E[6][6][6]), NULL, 0.0, 0, "Coefficient of polynomial for type 2 drifts"},

    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&hkpoly_example.factor), NULL, 1.0, 0, "additional factor to apply"},
    {"N_REPEATS", "", IS_LONG, 0, (long)((char *)&hkpoly_example.nRepeats), NULL, 0.0, 1, "Number of times to repeat the drift-kick-drift sequence. Strength of each application is reduced by this factor."},
    {"DRIFT_TYPE", "", IS_SHORT, 0, (long)((char *)&hkpoly_example.driftType), NULL, 0.0, 1, "If 1, then use D[i][j]. If 2, then use E[i][j][k]."},
    };

NISEPT nisept_example;
/* integrated septum magnet physical parameters */
PARAMETER nisept_param[N_NISEPT_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nisept_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nisept_example.angle), NULL, 0.0, 0, "bend angle"},
    {"E1", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nisept_example.e1), NULL, 0.0, 0, "entrance edge angle"},
    {"B1", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nisept_example.b1), NULL, 0.0, 0, "normalized gradient (K1=B1*L/ANGLE)"},
    {"Q1REF", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nisept_example.q1_ref), NULL, 0.0, 0, "distance from septum at which bending radius is L/ANGLE"},
    {"FLEN", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nisept_example.flen), NULL, 0.0, 0, "fringe field length"},
    {"ACCURACY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&nisept_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"METHOD", "", IS_STRING, 0, (long)((char *)&nisept_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, modified-midpoint, two-pass modified-midpoint, leap-frog, non-adaptive runge-kutta"},
    {"MODEL", "", IS_STRING, 0, (long)((char *)&nisept_example.model), DEFAULT_NIBEND_TYPE, 0.0, 0, "fringe model (hard-edge, linear, cubic-spline, tanh, quintic"},
    };

RAMPRF ramprf_example;
/* ramped rf cavity physical parameters */
PARAMETER ramprf_param[N_RAMPRF_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ramprf_example.length), NULL, 0.0, 0, "length"},
    {"VOLT", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ramprf_example.volt), NULL, 0.0, 0, "nominal voltage"},
    {"PHASE", "DEG", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ramprf_example.phase), NULL, 0.0, 0, "nominal phase"},
    {"FREQ", "Hz", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ramprf_example.freq), NULL, 500.0e6, 0, "nominal frequency"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&ramprf_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"VOLT_WAVEFORM", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&ramprf_example.vwaveform), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving voltage waveform factor vs time"},
    {"PHASE_WAVEFORM", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&ramprf_example.pwaveform), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving phase offset vs time (requires FREQ_WAVEFORM)"},
    {"FREQ_WAVEFORM", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&ramprf_example.fwaveform), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving frequencyfactor vs time (requires PHASE_WAVEFORM)"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&ramprf_example.fiducial), NULL, 0.0, 0, "mode for determining fiducial arrival time (light, tmean, first, pmaximum)"},
    };

/* momentum ramp physical parameters */
PARAMETER rampp_param[N_RAMPP_PARAMS] = {
    {"WAVEFORM", "", IS_STRING, PARAM_XY_WAVEFORM, 0, NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving momentum factor vs time"}
    };

STRAY stray_example;
/* stray field physical parameters */
PARAMETER stray_param[N_STRAY_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&stray_example.length), NULL, 0.0, 0, "length"},
    {"LBX", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&stray_example.lBx), NULL, 0.0, 0, "local Bx"},
    {"LBY", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&stray_example.lBy), NULL, 0.0, 0, "local By"},
    {"GBX", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&stray_example.gBx), NULL, 0.0, 0, "global Bx"},
    {"GBY", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&stray_example.gBy), NULL, 0.0, 0, "global By"},
    {"GBZ", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&stray_example.gBz), NULL, 0.0, 0, "global Bz"},
    {"ORDER", "", IS_LONG, 0, (long)((char *)&stray_example.order), NULL, 0.0, 0, "matrix order"},
    };

CSBEND csbend_example;
/* canonically-integrated sector bending magnet physical parameters */
PARAMETER csbend_param[N_CSBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&csbend_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&csbend_example.angle), NULL, 0.0, 0, "bend angle"},
    {"K1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.k1), NULL, 0.0, 0, "geometric quadrupole strength"},
    {"K2", "1/M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.k2), NULL, 0.0, 0, "geometric sextupole strength"},
    {"K3", "1/M$a4$n", IS_DOUBLE, 0, (long)((char *)&csbend_example.k3), NULL, 0.0, 0, "geometric octupole strength"},
    {"K4", "1/M$a5$n", IS_DOUBLE, 0, (long)((char *)&csbend_example.k4), NULL, 0.0, 0, "geometric decapole strength"},
    {"K5", "1/M$a6$n", IS_DOUBLE, 0, (long)((char *)&csbend_example.k5), NULL, 0.0, 0, "geometric 12-pole strength"},
    {"K6", "1/M$a7$n", IS_DOUBLE, 0, (long)((char *)&csbend_example.k6), NULL, 0.0, 0, "geometric 14-pole strength"},
    {"K7", "1/M$a8$n", IS_DOUBLE, 0, (long)((char *)&csbend_example.k7), NULL, 0.0, 0, "geometric 16-pole strength"},
    {"K8", "1/M$a9$n", IS_DOUBLE, 0, (long)((char *)&csbend_example.k8), NULL, 0.0, 0, "geometric 18-pole strength"},
    {"E1", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.e[0]), NULL, 0.0, 0, "entrance edge angle"},
    {"E2", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.e[1]), NULL, 0.0, 0, "exit edge angle"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"H1", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.h[0]), NULL, 0.0, 0, "entrance pole-face curvature"},
    {"H2", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.h[1]), NULL, 0.0, 0, "exit pole-face curvature"},
    {"HGAP", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.hgap), NULL, 0.0, 0, "half-gap between poles"},
    {"FINT", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.fintBoth), NULL, DEFAULT_FINT, 0, "edge-field integral"},
    {"FINT1", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.fint[0]), NULL, -1, 0, "edge-field integral. If negative, use FINT."},
    {"FINT2", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.fint[1]), NULL, -1, 0, "edge-field integral. If negative, use FINT."},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.fse), NULL, 0.0, 0, "fractional strength error of all components"},
    {"FSE_DIPOLE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.fseDipole), NULL, 0.0, 0, "fractional strength error of dipole component"},
    {"FSE_QUADRUPOLE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.fseQuadrupole), NULL, 0.0, 0, "fractional strength error of quadrupole component"},
    {"ETILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.etilt), NULL, 0.0, 0, "error rotation about incoming longitudinal axis"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&csbend_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks"},
    {"NONLINEAR", "", IS_SHORT, 0, (long)((char *)&csbend_example.nonlinear), NULL, 0.0, 1, "include nonlinear field components?"},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&csbend_example.synch_rad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"EDGE1_EFFECTS", "", IS_SHORT, 0, (long)((char *)&csbend_example.edge_effects[0]), NULL, 0.0, 1, "include entrance edge effects?"},
    {"EDGE2_EFFECTS", "", IS_SHORT, 0, (long)((char *)&csbend_example.edge_effects[1]), NULL, 0.0, 1, "include exit edge effects?"},
    {"EDGE_ORDER", "", IS_SHORT, 0, (long)((char *)&csbend_example.edge_order), NULL, 0.0, 1, "order to which to include edge effects"},
    {"INTEGRATION_ORDER", "", IS_SHORT, 0, (long)((char *)&csbend_example.integration_order), NULL, 0.0, 4, "integration order (2 or 4)"},
    {"EXPAND_HAMILTONIAN", "", IS_SHORT, 0, (long)((char *)&csbend_example.expandHamiltonian), NULL, 0.0, 0, "If 1, Hamiltonian is expanded to leading order."},
    {"EDGE1_KICK_LIMIT", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.edge_kick_limit[0]), NULL, -1., 0, "maximum kick entrance edge can deliver"},
    {"EDGE2_KICK_LIMIT", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.edge_kick_limit[1]), NULL, -1., 0, "maximum kick exit edge can deliver"},
    {"KICK_LIMIT_SCALING", "", IS_SHORT, 0, (long)((char *)&csbend_example.kick_limit_scaling), NULL, 0, 0, "scale maximum edge kick with FSE?"},
    {"USE_BN", "", IS_SHORT, 0, (long)((char *)&csbend_example.use_bn), NULL, 0.0, 0, "use b<n> instead of K<n>?"},
    {"EXPANSION_ORDER", "", IS_SHORT, 0, (long)((char *)&csbend_example.expansionOrder), NULL, 0.0, 0, "Order of field expansion. (0=auto)"},
    {"B1", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.b1), NULL, 0.0, 0, "K1 = b1/rho, where rho is bend radius"},
    {"B2", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.b2), NULL, 0.0, 0, "K2 = b2/rho"},
    {"B3", "1/M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.b3), NULL, 0.0, 0, "K3 = b3/rho"},
    {"B4", "1/M$a4$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.b4), NULL, 0.0, 0, "K4 = b4/rho"},
    {"B5", "1/M$a5$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.b5), NULL, 0.0, 0, "K5 = b5/rho"},
    {"B6", "1/M$a6$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.b6), NULL, 0.0, 0, "K6 = b6/rho"},
    {"B7", "1/M$a7$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.b7), NULL, 0.0, 0, "K7 = b7/rho"},
    {"B8", "1/M$a8$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csbend_example.b8), NULL, 0.0, 0, "K8 = b8/rho"},
    {"XREFERENCE", "M", IS_DOUBLE, 0, (long)((char *)&csbend_example.xReference), NULL, 0.0, 0, "reference x for interpretation of fn values"},
    {"F1", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.f1), NULL, 0.0, 0, "Fractional normal field error fn=bn*xr^n/n!, adds to K1 or b1."},
    {"F2", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.f2), NULL, 0.0, 0, "Fractional normal field error fn=bn*xr^n/n!, adds to K2 or b2."},
    {"F3", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.f3), NULL, 0.0, 0, "Fractional normal field error fn=bn*xr^n/n!, additive."},
    {"F4", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.f4), NULL, 0.0, 0, "Fractional normal field error fn=bn*xr^n/n!, additive."},
    {"F5", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.f5), NULL, 0.0, 0, "Fractional normal field error fn=bn*xr^n/n!, additive."},
    {"F6", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.f6), NULL, 0.0, 0, "Fractional normal field error fn=bn*xr^n/n!, additive."},
    {"F7", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.f7), NULL, 0.0, 0, "Fractional normal field error fn=bn*xr^n/n!, additive."},
    {"F8", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.f8), NULL, 0.0, 0, "Fractional normal field error fn=bn*xr^n/n!, additive."},
    {"G1", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.g1), NULL, 0.0, 0, "Fractional skew field error."},
    {"G2", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.g2), NULL, 0.0, 0, "Fractional skew field error."},
    {"G3", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.g3), NULL, 0.0, 0, "Fractional skew field error."},
    {"G4", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.g4), NULL, 0.0, 0, "Fractional skew field error."},
    {"G5", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.g5), NULL, 0.0, 0, "Fractional skew field error."},
    {"G6", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.g6), NULL, 0.0, 0, "Fractional skew field error."},
    {"G7", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.g7), NULL, 0.0, 0, "Fractional skew field error."},
    {"G8", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.g8), NULL, 0.0, 0, "Fractional skew field error."},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&csbend_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"ISR1PART", "", IS_SHORT, 0, (long)((char *)&csbend_example.isr1Particle), NULL, 0.0, 1, "Include ISR for single-particle beam only if ISR=1 and ISR1PART=1"},
    {"SQRT_ORDER", "", IS_SHORT, 0, (long)((char *)&csbend_example.sqrtOrder), NULL, 0.0, 0, "Ignored, kept for backward compatibility only."},
    {"USE_RAD_DIST", "", IS_SHORT, 0, (long)((char *)&csbend_example.distributionBasedRadiation), NULL, 0.0, 0, "If nonzero, overrides SYNCH_RAD and ISR, causing simulation of radiation from distributions, optionally including opening angle."},
    {"ADD_OPENING_ANGLE", "", IS_SHORT, 0, (long)((char *)&csbend_example.includeOpeningAngle), NULL, 0.0, 1, "If nonzero, radiation opening angle effects are added if USE_RAD_DIST is nonzero."},
    {"PHOTON_OUTPUT_FILE", "", IS_STRING, 0, (long)((char *)&csbend_example.photonOutputFile), NULL, 0.0, 0, "output file for photons, if USE_RAD_DIST=1"},
    {"PHOTON_LOW_ENERGY_CUTOFF", "eV", IS_DOUBLE, 0, (long)((char *)&csbend_example.photonLowEnergyCutoff), NULL, 0.0, 0, "Lower limit of photon energy to output."},
    {"REFERENCE_CORRECTION", "", IS_SHORT, 0, (long)((char *)&csbend_example.referenceCorrection), NULL, 0.0, 0, "If nonzero, reference trajectory is subtracted from particle trajectories to compensate for inaccuracy in integration."},
    {"TRACKING_MATRIX", "", IS_SHORT, 0, (long)((char *)&csbend_example.trackingMatrix), NULL, 0.0, 0, "If nonzero, gives order of tracking-based matrix up to third order to be used for twiss parameters etc.  If zero, 2nd-order analytical matrix is used."},
    {"FSE_CORRECTION", "", IS_SHORT, 0, (long)((char *)&csbend_example.fseCorrection), NULL, 0.0, 0, "If nonzero, FSE is adjusted to compensate for edge effects when EDGE1_EFFECTS or EDGE2_EFFECTS = 2"},
    };


CSRCSBEND csrcsbend_example;
/* canonically-integrated sector bending magnet with CSR physical parameters */
PARAMETER csrcsbend_param[N_CSRCSBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.angle), NULL, 0.0, 0, "bend angle"},
    {"K1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.k1), NULL, 0.0, 0, "geometric quadrupole strength"},
    {"K2", "1/M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.k2), NULL, 0.0, 0, "geometric sextupole strength"},
    {"K3", "1/M$a4$n", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.k3), NULL, 0.0, 0, "geometric octupole strength"},
    {"K4", "1/M$a5$n", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.k4), NULL, 0.0, 0, "geometric decapole strength"},
    {"K5", "1/M$a6$n", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.k5), NULL, 0.0, 0, "geometric 12-pole strength"},
    {"K6", "1/M$a7$n", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.k6), NULL, 0.0, 0, "geometric 14-pole strength"},
    {"K7", "1/M$a8$n", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.k7), NULL, 0.0, 0, "geometric 16-pole strength"},
    {"K8", "1/M$a9$n", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.k8), NULL, 0.0, 0, "geometric 18-pole strength"},
    {"E1", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.e[0]), NULL, 0.0, 0, "entrance edge angle"},
    {"E2", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.e[1]), NULL, 0.0, 0, "exit edge angle"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"H1", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.h[0]), NULL, 0.0, 0, "entrance pole-face curvature"},
    {"H2", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.h[1]), NULL, 0.0, 0, "exit pole-face curvature"},
    {"HGAP", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.hgap), NULL, 0.0, 0, "half-gap between poles"},
    {"FINT", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.fint), NULL, DEFAULT_FINT, 0, "edge-field integral"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"ETILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.etilt), NULL, 0.0, 0, "error rotation about incoming longitudinal axis"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks"},
    {"NONLINEAR", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.nonlinear), NULL, 0.0, 1, "include nonlinear field components?"},
    {"LINEARIZE", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.useMatrix), NULL, 0.0, 0, "use linear matrix instead of symplectic integrator?"},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.synch_rad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"EDGE1_EFFECTS", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.edge_effects[0]), NULL, 0.0, 1, "include entrance edge effects?"},
    {"EDGE2_EFFECTS", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.edge_effects[1]), NULL, 0.0, 1, "include exit edge effects?"},
    {"EDGE_ORDER", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.edge_order), NULL, 0.0, 1, "order to which to include edge effects"},
    {"INTEGRATION_ORDER", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.integration_order), NULL, 0.0, 4, "integration order (2 or 4)"},
    {"BINS", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.bins), NULL, 0.0, 0, "number of bins for CSR wake"},
    {"BIN_ONCE", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.binOnce), NULL, 0.0, 0, "bin only at the start of the dipole?"},
    {"BIN_RANGE_FACTOR", "", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.binRangeFactor), NULL, 1.2, 0, "Factor by which to increase the range of histogram compared to total bunch length.  Large value eliminates binning problems in CSRDRIFTs."},
    {"SG_HALFWIDTH", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.SGHalfWidth), NULL, 0.0, 0, "Savitzky-Golay filter half-width for smoothing current histogram.  If less than 1, no SG smoothing is performed."},
    {"SG_ORDER", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing current histogram"},
    {"SGDERIV_HALFWIDTH", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.SGDerivHalfWidth), NULL, 0.0, 0, "Savitzky-Golay filter half-width for taking derivative of current histogram.  Defaults to SG_HALFWIDTH (if positive) or else 1."},
    {"SGDERIV_ORDER", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.SGDerivOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for taking derivative of current histogram"},
    {"TRAPAZOID_INTEGRATION", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.trapazoidIntegration), NULL, 0.0, 1, "Select whether to use trapazoid-rule integration (default) or a simple sum."},
    {"OUTPUT_FILE", "", IS_STRING, 0, (long)((char *)&csrcsbend_example.histogramFile), NULL, 0.0, 0, "output file for CSR wakes"},
    {"OUTPUT_INTERVAL", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.outputInterval), NULL, 0.0, 1, "interval (in kicks) of output to OUTPUT_FILE"},
    {"OUTPUT_LAST_WAKE_ONLY", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.outputLastWakeOnly), NULL, 0.0, 0, "output final wake only?"},
    {"STEADY_STATE", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.steadyState), NULL, 0.0, 0, "use steady-state wake equations?"},
    {"IGF", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.integratedGreensFunction), NULL, 0.0, 0, "use integrated Greens function (requires STEADY_STATE=1)?"},
    {"USE_BN", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.use_bn), NULL, 0.0, 0, "use b<n> instead of K<n>?"},
    {"EXPANSION_ORDER", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.expansionOrder), NULL, 0.0, 0, "Order of field expansion. (0=auto)"},
    {"B1", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.b1), NULL, 0.0, 0, "K1 = b1/rho, where rho is bend radius"},
    {"B2", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.b2), NULL, 0.0, 0, "K2 = B2/rho"},
    {"B3", "1/M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.b3), NULL, 0.0, 0, "K3 = B3/rho"},
    {"B4", "1/M$a4$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.b4), NULL, 0.0, 0, "K4 = B4/rho"},
    {"B5", "1/M$a5$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.b5), NULL, 0.0, 0, "K5 = B5/rho"},
    {"B6", "1/M$a6$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.b6), NULL, 0.0, 0, "K6 = B6/rho"},
    {"B7", "1/M$a7$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.b7), NULL, 0.0, 0, "K7 = B7/rho"},
    {"B8", "1/M$a8$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&csrcsbend_example.b8), NULL, 0.0, 0, "K8 = B8/rho"},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"ISR1PART", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.isr1Particle), NULL, 0.0, 1, "Include ISR for single-particle beam only if ISR=1 and ISR1PART=1"},
    {"CSR", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.csr), NULL, 0.0, 1, "enable CSR computations?"},
    {"BLOCK_CSR", "", IS_SHORT, 0, (long)((char *)&csrcsbend_example.csrBlock), NULL, 0.0, 0, "block CSR from entering CSRDRIFT?"},
    {"DERBENEV_CRITERION_MODE", "", IS_STRING, 0, (long)((char *)&csrcsbend_example.derbenevCriterionMode), "disable", 0.0, 1, "disable, evaluate, or enforce"},
    {"PARTICLE_OUTPUT_FILE", "", IS_STRING, 0, (long)((char*)&csrcsbend_example.particleOutputFile), NULL, 0.0, 0, "name of file for phase-space output"},
    {"PARTICLE_OUTPUT_INTERVAL", "", IS_LONG, 0, (long)((char*)&csrcsbend_example.particleOutputInterval), NULL, 0.0, 0, "interval (in kicks) of output to PARTICLE_OUTPUT_FILE"},
    {"SLICE_ANALYSIS_INTERVAL", "", IS_LONG, 0, (long)((char*)&csrcsbend_example.sliceAnalysisInterval), NULL, 0.0, 0, "interval (in kicks) of output to slice analysis file (from slice_analysis command)"},    
    {"LOW_FREQUENCY_CUTOFF0", "", IS_DOUBLE, 0, (long)((char*)&csrcsbend_example.lowFrequencyCutoff0), NULL, -1.0, 0, "Highest spatial frequency at which low-frequency cutoff filter is zero.  If not positive, no low-frequency cutoff filter is applied. Frequency is in units of Nyquist (0.5/binsize)."},
    {"LOW_FREQUENCY_CUTOFF1", "", IS_DOUBLE, 0, (long)((char*)&csrcsbend_example.lowFrequencyCutoff1), NULL, -1.0, 0, "Lowest spatial frequency at which low-frequency cutoff filter is 1.  If not given, defaults to LOW_FREQUENCY_CUTOFF1."},
    {"HIGH_FREQUENCY_CUTOFF0", "", IS_DOUBLE, 0, (long)((char*)&csrcsbend_example.highFrequencyCutoff0), NULL, -1.0, 0, "Spatial frequency at which smoothing (high-frequency cutoff) filter begins.  If not positive, no frequency filter smoothing is done.  Frequency is in units of Nyquist (0.5/binsize)."},
    {"HIGH_FREQUENCY_CUTOFF1", "", IS_DOUBLE, 0, (long)((char*)&csrcsbend_example.highFrequencyCutoff1), NULL, -1.0, 0, "Spatial frequency at which smoothing (high-frequency cutoff) filter is 0.  If not given, defaults to HIGH_FREQUENCY_CUTOFF0."},
    {"CLIP_NEGATIVE_BINS", "", IS_SHORT, 0, (long)((char*)&csrcsbend_example.clipNegativeBins), NULL, 0.0, 1, "If non-zero, then any bins with negative counts after the filters are applied have the counts set to zero."},    
    {"WAKE_FILTER_FILE", "", IS_STRING, 0, (long)((char*)&csrcsbend_example.wakeFilterFile), NULL, 0.0, 0, "Name of file supplying wakefield filtering data."},
    {"WFF_FREQ_COLUMN", "", IS_STRING, 0, (long)((char*)&csrcsbend_example.wffFreqColumn), NULL, 0.0, 0, "Name of column supplying frequency values for wakefield filtering data."},
    {"WFF_REAL_COLUMN", "", IS_STRING, 0, (long)((char*)&csrcsbend_example.wffRealColumn), NULL, 0.0, 0, "Name of column supplying real values for wakefield filtering data."},
    {"WFF_IMAG_COLUMN", "", IS_STRING, 0, (long)((char*)&csrcsbend_example.wffImagColumn), NULL, 0.0, 0, "Name of column supplying imaginary values for wakefield filtering data."},
};

TUBEND tubend_example;
/* special bending magnet for top-up with entry through the side! */
PARAMETER tubend_param[N_TUBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&tubend_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&tubend_example.angle), NULL, 0.0, 0, "bend angle"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&tubend_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"OFFSET", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&tubend_example.offset), NULL, 0.0, 0, "horizontal offset of magnet center from arc center"},
    {"MAGNET_WIDTH", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&tubend_example.magnet_width), NULL, 0.0, 0, "horizontal width of the magnet pole"},
    {"MAGNET_ANGLE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&tubend_example.magnet_angle), NULL, 0.0, 0, "angle that the magnet was designed for"},
    };

TWMTA twmta_example;
/* names for traveling-wave muffin-tin linac parameters
 */
PARAMETER twmta_param[N_TWMTA_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.length), NULL, 0.0, 0, "length"},
    {"FREQUENCY", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.frequency), NULL, DEFAULT_FREQUENCY, 0, "frequency"},
    {"PHASE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.phase), NULL, 0.0, 0, "phase"},
    {"EZ", "V/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.Ez), NULL, 0.0, 0, "electric field"},
    {"ACCURACY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"X_MAX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.x_max), NULL, 0.0, 0, "x half-aperture"},
    {"Y_MAX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.y_max), NULL, 0.0, 0, "y half-aperture"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.dy), NULL, 0.0, 0, "misalignment"},
    {"KX", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.kx), NULL, 0.0, 0, "horizontal wave number"},
    {"BETA_WAVE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.beta_wave), NULL, DEFAULT_BETA_WAVE, 0, "(phase velocity)/c"},
    {"BSOL", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.Bsol), NULL, 0.0, 0, "solenoid field"},
    {"ALPHA", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.alpha), NULL, 0.0, 0, "field attenuation factor"},
    {"PHASE_REFERENCE", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"N_STEPS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&twmta_example.n_steps), NULL, 0.0, 100, "number of kicks"},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&twmta_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, non-adaptive runge-kutta, modified midpoint)"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&twmta_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0, "{t|p},{median|min|max|ave|first|light} (e.g., \"t,median\")"}
    } ;

MATTER matter_example;
/* matter physical parameters */
PARAMETER matter_param[N_MATTER_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&matter_example.length), NULL, 0.0, 0, "length"},
    {"LEFFECTIVE", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&matter_example.lEffective), NULL, 0.0, 0, "effective length (used if L=0)"},
    {"XO", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&matter_example.Xo), NULL, 0.0, 0, "radiation length"},
    {"ENERGY_DECAY", "", IS_LONG, 0, (long)((char *)&matter_example.energyDecay), NULL, 0.0, 0, "If nonzero, then particles will lose energy due to material using a simple exponential model."},
    {"ENERGY_STRAGGLE", "", IS_LONG, 0, (long)((char *)&matter_example.energyStraggle), NULL, 0.0, 0, "Use simple-minded energy straggling model coupled with ENERGY_DECAY=1?"},
    {"NUCLEAR_BREMSSTRAHLUNG", "", IS_LONG, 0, (long)((char *)&matter_example.nuclearBremsstrahlung), NULL, 0.0, 0, "Model energy loss to nuclear bremsstrahlung? If enabled, set ENERGY_DECAY=0 to disable simpler model."},
    {"ELECTRON_RECOIL", "", IS_LONG, 0, (long)((char *)&matter_example.electronRecoil), NULL, 0.0, 0, "If non-zero, electron recoil during Coulomb scattering is included (results in energy change)."},
    {"Z", "", IS_LONG, 0, (long)((char *)&matter_example.Z), NULL, 0.0, 0, "Atomic number"},
    {"A", "AMU", IS_DOUBLE, 0, (long)((char *)&matter_example.A), NULL, 0.0, 0, "Atomic mass"},
    {"RHO", "KG/M^3", IS_DOUBLE, 0, (long)((char *)&matter_example.rho), NULL, 0.0, 0, "Density"},       
    {"PLIMIT", "", IS_DOUBLE, 0, (long)((char *)&matter_example.pLimit), NULL, 0.05, 0, "Probability cutoff for each slice"},
    {"WIDTH", "M", IS_DOUBLE, 0, (long)((char *)&matter_example.width), NULL, 0.0, 0, "Full width of slots. If 0, no slots are present."},
    {"SPACING", "M", IS_DOUBLE, 0, (long)((char *)&matter_example.spacing), NULL, 0.0, 0, "Center-to-center spacing of slots. If 0, no slots are present."},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&matter_example.tilt), NULL, 0.0, 0, "Tilt of slot array about the longitudinal axis. "},
    {"CENTER", "M", IS_DOUBLE, 0, (long)((char *)&matter_example.center), NULL, 0.0, 0, "Position of center of slot array in rotated frame."},
    {"N_SLOTS", "", IS_LONG, 0, (long)((char *)&matter_example.nSlots), NULL, 0.0, 0, "Number of empty slots in material. If <=0, an infinite array is assumed."},
    {"START_PASS", "", IS_LONG, 0, (long)((char *)&matter_example.startPass), NULL, 0.0, -1, "If non-negative, pass on which to start interaction with beam."},    
    {"END_PASS", "", IS_LONG, 0, (long)((char *)&matter_example.endPass), NULL, 0.0, -1, "If non-negative, pass on which to end interaction with beam."},
    };

RFMODE rfmode_example;
/* RFMODE physical parameters */
PARAMETER rfmode_param[N_RFMODE_PARAMS] = {
    {"RA", "Ohm", IS_DOUBLE, 0, (long)((char *)&rfmode_example.Ra), NULL, 0.0, 0, "shunt impedance, Ra=V^2/P"},
    {"RS", "Ohm", IS_DOUBLE, 0, (long)((char *)&rfmode_example.Rs), NULL, 0.0, 0, "shunt impedance (Rs=Ra/2)"},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&rfmode_example.Q), NULL, 0.0, 1, "cavity Q"},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&rfmode_example.freq), NULL, 0.0, 0, "Resonant frequency of the cavity mode"},
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&rfmode_example.charge), NULL, 0.0, 0, "beam charge (or use CHARGE element)"},
    {"INITIAL_V", "V", IS_DOUBLE, 0, (long)((char *)&rfmode_example.initial_V), NULL, 0.0, 0, "initial beam-loading voltage"},
    {"INITIAL_PHASE", "RAD", IS_DOUBLE, 0, (long)((char *)&rfmode_example.initial_phase), NULL, 0.0, 0, "initial beam-loading phase"},
    {"INITIAL_T", "S", IS_DOUBLE, 0, (long)((char *)&rfmode_example.initial_t), NULL, 0.0, 0, "time at which INITIAL_V and INITIAL_PHASE held"},
    {"BETA", "", IS_DOUBLE, 0, (long)((char *)&rfmode_example.beta), NULL, 0.0, 0, "normalized load impedance"},
    {"BIN_SIZE", "S", IS_DOUBLE, 0, (long)((char *)&rfmode_example.bin_size), NULL, 0.0, 0, "bin size for current histogram (use 0 for autosize)"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&rfmode_example.n_bins), NULL, 0.0, 20, "number of bins for current histogram"},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&rfmode_example.interpolate), NULL, 0.0, 0, "if non-zero, interpolate voltage within bins"},
    {"PRELOAD", "", IS_LONG, 0, (long)((char *)&rfmode_example.preload), NULL, 0.0, 0, "preload cavity with steady-state field"},
    {"PRELOAD_CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&rfmode_example.preloadCharge), NULL, 0.0, 0, "beam charge used for preloading calculations"},
    {"PRELOAD_FACTOR", "", IS_DOUBLE, 0, (long)((char *)&rfmode_example.preload_factor), NULL, 1.0, 0, "multiply preloaded field by this value"},
    {"PRELOAD_HARMONIC", "", IS_LONG, 0, (long)((char *)&rfmode_example.preloadHarmonic), NULL, 0.0, 0, "If detuning from harmonic is greater than half the revolution frequency, automatic determination of the rf harmonic will fail. Give the harmonic explicitly with this parameter."},
    {"RIGID_UNTIL_PASS", "", IS_LONG, 0, (long)((char *)&rfmode_example.rigid_until_pass), NULL, 0.0, 0, "don't affect the beam until this pass"},
    {"DETUNED_UNTIL_PASS", "", IS_LONG, 0, (long)((char *)&rfmode_example.detuned_until_pass), NULL, 0.0, 0, "cavity is completely detuned until this pass"},
    {"SAMPLE_INTERVAL", "", IS_LONG, 0, (long)((char *)&rfmode_example.sample_interval), NULL, 0.0, 1, "passes between samples to RECORD file"},
    {"FLUSH_INTERVAL", "", IS_LONG, 0, (long)((char *)&rfmode_example.flush_interval), NULL, 0.0, 1000, "samples between flushing output to RECORD file"},
    {"RECORD", "", IS_STRING, 0, (long)((char *)&rfmode_example.record), NULL, 0.0, 0, "output file for cavity fields"},
    {"SINGLE_PASS", "", IS_LONG, 0, (long)((char *)&rfmode_example.single_pass), NULL, 0.0, 0, "if nonzero, don't accumulate field from pass to pass"},
    {"PASS_INTERVAL", "", IS_LONG, 0, (long)((char *)&rfmode_example.pass_interval), NULL, 0.0, 1, "interval in passes at which to apply PASS_INTERVAL times the field (may increase speed)"},
    {"FREQ_WAVEFORM", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&rfmode_example.fwaveform), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving frequency/f0 vs time, where f0 is the frequency given with the FREQ parameter"},
    {"Q_WAVEFORM", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&rfmode_example.Qwaveform), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving qualityFactor/Q0 vs time, where Q0 is the quality factor given the the Q parameter."},
    {"RAMP_PASSES", "", IS_LONG, 0, (long)((char *)&rfmode_example.rampPasses), NULL, 0.0, 0, "Number of passes over which to linearly ramp up the impedance to full strength."},
    {"BINLESS", "", IS_LONG, 0, (long)((char *)&rfmode_example.binless), NULL, 0.0, 0, "If nonzero, use algorithm that doesn't requiring binning.  Best for few particles, widely spaced."},
    {"RESET_FOR_EACH_STEP", "", IS_LONG, 0, (long)((char *)&rfmode_example.reset_for_each_step), NULL, 0.0, 1, "If nonzero, voltage and phase are reset for each simulation step."},
    {"LONG_RANGE_ONLY", "", IS_LONG, 0, (long)((char *)&rfmode_example.long_range_only), NULL, 0.0, 0, "If nonzero, induced voltage from present turn does not affect bunch. Results are not self-consistent!"},
    {"ALLOW_UNBINNED_PARTICLES", "", IS_LONG, 0, (long)((char *)&rfmode_example.allowUnbinnedParticles), NULL, 0.0, 0, "If nonzero, will keep running even if some particles fall outside the binning region. Use with caution!"},
    {"N_CAVITIES", "", IS_LONG, 0, (long)((char *)&rfmode_example.n_cavities), NULL, 0.0, 1, "effect is multiplied by this number, simulating N identical cavities"},
    {"BUNCHED_BEAM_MODE", "", IS_LONG, 0, (long)((char *)&rfmode_example.bunchedBeamMode), NULL, 0.0, 1, "If 1, then do calculations bunch-by-bunch. If >1, use pseudo bunches."},
    {"BUNCH_INTERVAL", "S", IS_DOUBLE, 0, (long)((char *)&rfmode_example.bunchInterval), NULL, 0.0, 0, "For pseudo-bunch mode, time between bunches."},
    {"DRIVE_FREQUENCY", "Hz", IS_DOUBLE, 0, (long)((char *)&rfmode_example.driveFrequency), NULL, 0.0, 0, "drive frequency from generator. If zero, no generator voltage is applied."},
    {"V_SETPOINT", "V", IS_DOUBLE, 0, (long)((char *)&rfmode_example.voltageSetpoint), NULL, 0.0, 0, "setpoint for total cavity voltage"},
    {"PHASE_SETPOINT", "DEG", IS_DOUBLE, 0, (long)((char *)&rfmode_example.phaseSetpoint), NULL, 0.0, 0, "setpoint for total cavity phase"},
    {"UPDATE_INTERVAL", "", IS_LONG, 0, (long)((char *)&rfmode_example.updateInterval), NULL, 0.0, 1, "update interval of feedback in units of rf period"},
    {"READ_OFFSET", "", IS_LONG, 0, (long)((char *)&rfmode_example.readOffset), NULL, 0.0, 0, "Offset in buckets of point at which voltage and phase are read for feedback relative to the first bunch passage. A positive value corresponds to reading before bunch passage."},
    {"ADJUSTMENT_START", "", IS_LONG, 0, (long)((char *)&rfmode_example.adjustmentStart), NULL, 0.0, 0, "Pass on which to begin adjustment of the effective voltage setpoint."},
    {"ADJUSTMENT_END", "", IS_LONG, 0, (long)((char *)&rfmode_example.adjustmentEnd), NULL, 0.0, 0, "Pass on which to stop adjustment of the effective voltage setpoint."},
    {"ADJUSTMENT_INTERVAL", "", IS_LONG, 0, (long)((char *)&rfmode_example.adjustmentInterval), NULL, 0.0, 100, "Interval in passes between adjustment of the effective voltage setpoint."},
    {"ADJUSTMENT_FRACTION", "", IS_DOUBLE, 0, (long)((char *)&rfmode_example.adjustmentFraction), NULL, 0.0, 0, "Fraction of voltage setpoint error taken out on each adjustment step"},
    {"AMPLITUDE_FILTER", "", IS_STRING, 0, (long)((char *)&rfmode_example.amplitudeFilterFile), NULL, 0.0, 0, "IIR filter specification for amplitude feedback"},
    {"PHASE_FILTER", "", IS_STRING, 0, (long)((char *)&rfmode_example.phaseFilterFile), NULL, 0.0, 0, "IIR filter specification for phase feedback"},
    {"IN_PHASE_FILTER", "", IS_STRING, 0, (long)((char *)&rfmode_example.IFilterFile), NULL, 0.0, 0, "IIR filter specification for in-phase component feedback"},
    {"QUADRATURE_FILTER", "", IS_STRING, 0, (long)((char *)&rfmode_example.QFilterFile), NULL, 0.0, 0, "IIR filter specification for quadrature component feedback"},
    {"FEEDBACK_RECORD", "", IS_STRING, 0, (long)((char *)&rfmode_example.feedbackRecordFile), NULL, 0.0, 0, "output file for feedback data"},
    {"MUTE_GENERATOR", "", IS_LONG, 0, (long)((char *)&rfmode_example.muteGenerator), NULL, 0.0, -1, "If nonnegative, gives the pass on which to mute the generator. This simulates an rf trip."},
    {"NOISE_ALPHA_GEN", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&rfmode_example.noiseData[I_NOISE_ALPHA_GEN]), NULL, 0.0, 0, "<filename>=<x>+<y> specifying alpha(t) for generator noise."},
    {"NOISE_PHI_GEN", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&rfmode_example.noiseData[I_NOISE_PHI_GEN]), NULL, 0.0, 0, "<filename>=<x>+<y> specifying dphi(t) for generator noise, in radians."},
    {"NOISE_ALPHA_V", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&rfmode_example.noiseData[I_NOISE_ALPHA_V]), NULL, 0.0, 0, "<filename>=<x>+<y> specifying alpha(t) for voltage noise."},
    {"NOISE_PHI_V", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&rfmode_example.noiseData[I_NOISE_PHI_V]), NULL, 0.0, 0, "<filename>=<x>+<y> specifying dphi(t) for voltage noise, in radians."},
    {"NOISE_I_GEN", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&rfmode_example.noiseData[I_NOISE_I_GEN]), NULL, 0.0, 0, "<filename>=<x>+<y> specifying ni(t) for in-phase generator noise."},
    {"NOISE_Q_GEN", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&rfmode_example.noiseData[I_NOISE_Q_GEN]), NULL, 0.0, 0, "<filename>=<x>+<y> specifying nq(t) for quadrature generator noise."},
    {"NOISE_I_V", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&rfmode_example.noiseData[I_NOISE_I_V]), NULL, 0.0, 0, "<filename>=<x>+<y> specifying ei(t) for in-phase voltage noise."},
    {"NOISE_Q_V", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&rfmode_example.noiseData[I_NOISE_Q_V]), NULL, 0.0, 0, "<filename>=<x>+<y> specifying eq(t) for quadrature voltage noise."},
    };

FRFMODE frfmode_example;
/* FRFMODE physical parameters */
PARAMETER frfmode_param[N_FRFMODE_PARAMS] = {
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&frfmode_example.filename), NULL, 0.0, 0, "input file"},
    {"BIN_SIZE", "S", IS_DOUBLE, 0, (long)((char *)&frfmode_example.bin_size), NULL, 0.0, 0, "bin size for current histogram (use 0 for autosize)"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&frfmode_example.n_bins), NULL, 0.0, 20, "number of bins for current histogram"},
    {"RIGID_UNTIL_PASS", "", IS_LONG, 0, (long)((char *)&frfmode_example.rigid_until_pass), NULL, 0.0, 0, "don't affect the beam until this pass"},
    {"USE_SYMM_DATA", "", IS_LONG, 0, (long)((char *)&frfmode_example.useSymmData), NULL, 0.0, 0, "use \"Symm\" columns from URMEL output file?"},
    {"FACTOR", "", IS_DOUBLE, 0, (long)((char *)&frfmode_example.factor), NULL, 1.0, 0, "factor by which to multiply shunt impedances"},
    {"CUTOFF", "HZ", IS_DOUBLE, 0, (long)((char *)&frfmode_example.cutoffFrequency), NULL, 0.0, 0, "If >0, cutoff frequency.  Modes above this frequency are ignored."},
    {"OUTPUT_FILE", "", IS_STRING, 0, (long)((char *)&frfmode_example.outputFile), NULL, 0.0, 0, "Output file for voltage in each mode."},
    {"FLUSH_INTERVAL", "", IS_LONG, 0, (long)((char *)&frfmode_example.flushInterval), NULL, 0.0, 1, "Interval in passes at which to flush output data."},
    {"RAMP_PASSES", "", IS_LONG, 0, (long)((char *)&frfmode_example.rampPasses), NULL, 0.0, 0, "Number of passes over which to linearly ramp up the impedance to full strength."},
    {"RESET_FOR_EACH_STEP", "", IS_LONG, 0, (long)((char *)&frfmode_example.reset_for_each_step), NULL, 0.0, 1, "If nonzero, voltage and phase are reset for each simulation step."},
    {"LONG_RANGE_ONLY", "", IS_LONG, 0, (long)((char *)&frfmode_example.long_range_only), NULL, 0.0, 0, "If nonzero, induced voltage from present turn does not affect bunch. Short range wake should be included via WAKE or ZLONGIT element."},
    {"N_CAVITIES", "", IS_LONG, 0, (long)((char *)&frfmode_example.n_cavities), NULL, 0.0, 1, "effect is multiplied by this number, simulating N identical cavities"},
    {"BUNCHED_BEAM_MODE", "", IS_LONG, 0, (long)((char *)&frfmode_example.bunchedBeamMode), NULL, 0.0, 1, "If non-zero, then do calculations bunch-by-bunch."},
    };

TRFMODE trfmode_example;
/* TRFMODE physical parameters */
PARAMETER trfmode_param[N_TRFMODE_PARAMS] = {
    {"RA", "Ohm/m", IS_DOUBLE, 0, (long)((char *)&trfmode_example.Ra), NULL, 0.0, 0, "shunt impedance, Ra=V^2/P"},
    {"RS", "Ohm/m", IS_DOUBLE, 0, (long)((char *)&trfmode_example.Rs), NULL, 0.0, 0, "shunt impedance (Rs=Ra/2)"},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&trfmode_example.Q), NULL, 0.0, 1, "cavity Q"},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&trfmode_example.freq), NULL, 0.0, 0, "frequency"},
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&trfmode_example.charge), NULL, 0.0, 0, "beam charge (or use CHARGE element)"},
    {"BETA", "", IS_DOUBLE, 0, (long)((char *)&trfmode_example.beta), NULL, 0.0, 0, "normalized load impedance"},
    {"BIN_SIZE", "S", IS_DOUBLE, 0, (long)((char *)&trfmode_example.bin_size), NULL, 0.0, 0, "bin size for current histogram (use 0 for autosize)"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&trfmode_example.n_bins), NULL, 0.0, 20, "number of bins for current histogram"},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&trfmode_example.interpolate), NULL, 0.0, 0, "if non-zero, interpolate voltage within bins"},
    {"PLANE", "", IS_STRING, 0, (long)((char *)&trfmode_example.plane), "both", 0.0, 0, "x, y, or both"},
    {"SAMPLE_INTERVAL", "", IS_LONG, 0, (long)((char *)&trfmode_example.sample_interval), NULL, 0.0, 1, "passes between output to RECORD file"},
    {"PER_PARTICLE_OUTPUT", "", IS_LONG, 0, (long)((char *)&trfmode_example.perParticleOutput), NULL, 0.0, 0, "If non-zero, then in BINLESS mode, provides per-particle output of RECORD data."},
    {"RECORD", "", IS_STRING, 0, (long)((char *)&trfmode_example.record), NULL, 0.0, 0, "output file for cavity data"},
    {"SINGLE_PASS", "", IS_LONG, 0, (long)((char *)&trfmode_example.single_pass), NULL, 0.0, 0, "if nonzero, don't accumulate field from pass to pass"},
    {"RIGID_UNTIL_PASS", "", IS_LONG, 0, (long)((char *)&trfmode_example.rigid_until_pass), NULL, 0.0, 0, "don't affect the beam until this pass"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&trfmode_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&trfmode_example.dy), NULL, 0.0, 0, "misalignment"},
    {"XFACTOR", "", IS_DOUBLE, 0, (long)((char *)&trfmode_example.xfactor), NULL, 1.0, 0, "factor by which to multiply shunt impedances"},
    {"YFACTOR", "", IS_DOUBLE, 0, (long)((char *)&trfmode_example.yfactor), NULL, 1.0, 0, "factor by which to multiply shunt impedances"},
    {"RAMP_PASSES", "", IS_LONG, 0, (long)((char *)&trfmode_example.rampPasses), NULL, 0.0, 0, "Number of passes over which to linearly ramp up the impedance to full strength."},
    {"BINLESS", "", IS_LONG, 0, (long)((char *)&trfmode_example.binless), NULL, 0.0, 0, "If nonzero, use algorithm that doesn't requiring binning.  Best for few particles, widely spaced."},
    {"RESET_FOR_EACH_STEP", "", IS_LONG, 0, (long)((char *)&trfmode_example.reset_for_each_step), NULL, 0.0, 1, "If nonzero, voltage and phase are reset for each simulation step."},
    {"LONG_RANGE_ONLY", "", IS_LONG, 0, (long)((char *)&trfmode_example.long_range_only), NULL, 0.0, 0, "If nonzero, induced voltage from present turn does not affect bunch. Short range wake should be included via TRWAKE or ZTRANSVERSE element."},
    {"N_CAVITIES", "", IS_LONG, 0, (long)((char *)&trfmode_example.n_cavities), NULL, 0.0, 1, "effect is multiplied by this number, simulating N identical cavities"},
    {"BUNCHED_BEAM_MODE", "", IS_LONG, 0, (long)((char *)&trfmode_example.bunchedBeamMode), NULL, 0.0, 1, "If non-zero, then do calculations bunch-by-bunch."},
    };

FTRFMODE ftrfmode_example;
/* FTRFMODE physical parameters */
PARAMETER ftrfmode_param[N_FTRFMODE_PARAMS] = {
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&ftrfmode_example.filename), NULL, 0.0, 0, "input file"},
    {"BIN_SIZE", "S", IS_DOUBLE, 0, (long)((char *)&ftrfmode_example.bin_size), NULL, 0.0, 0, "bin size for current histogram (use 0 for autosize)"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&ftrfmode_example.n_bins), NULL, 0.0, 20, "number of bins for current histogram"},
    {"RIGID_UNTIL_PASS", "", IS_LONG, 0, (long)((char *)&ftrfmode_example.rigid_until_pass), NULL, 0.0, 0, "don't affect the beam until this pass"},
    {"USE_SYMM_DATA", "", IS_LONG, 0, (long)((char *)&ftrfmode_example.useSymmData), NULL, 0.0, 0, "use \"Symm\" columns from URMEL output file?"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&ftrfmode_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&ftrfmode_example.dy), NULL, 0.0, 0, "misalignment"},
    {"XFACTOR", "", IS_DOUBLE, 0, (long)((char *)&ftrfmode_example.xfactor), NULL, 1.0, 0, "factor by which to multiply shunt impedances"},
    {"YFACTOR", "", IS_DOUBLE, 0, (long)((char *)&ftrfmode_example.yfactor), NULL, 1.0, 0, "factor by which to multiply shunt impedances"},
    {"CUTOFF", "HZ", IS_DOUBLE, 0, (long)((char *)&ftrfmode_example.cutoffFrequency), NULL, 0.0, 0, "If >0, cutoff frequency.  Modes above this frequency are ignored."},
    {"OUTPUT_FILE", "", IS_STRING, 0, (long)((char *)&ftrfmode_example.outputFile), NULL, 0.0, 0, "Output file for voltage in each mode."},
    {"FLUSH_INTERVAL", "", IS_LONG, 0, (long)((char *)&ftrfmode_example.flushInterval), NULL, 0.0, 1, "Interval in passes at which to flush output data."},
    {"RAMP_PASSES", "", IS_LONG, 0, (long)((char *)&ftrfmode_example.rampPasses), NULL, 0.0, 0, "Number of passes over which to linearly ramp up the impedance to full strength."},
    {"RESET_FOR_EACH_STEP", "", IS_LONG, 0, (long)((char *)&ftrfmode_example.reset_for_each_step), NULL, 0.0, 1, "If nonzero, voltage and phase are reset for each simulation step."},
    {"LONG_RANGE_ONLY", "", IS_LONG, 0, (long)((char *)&ftrfmode_example.long_range_only), NULL, 0.0, 0, "If nonzero, induced voltage from present turn does not affect bunch. Short range wake should be included via WAKE or ZLONGIT element."},
    {"N_CAVITIES", "", IS_LONG, 0, (long)((char *)&ftrfmode_example.n_cavities), NULL, 0.0, 1, "effect is multiplied by this number, simulating N identical cavities"},
    {"BUNCHED_BEAM_MODE", "", IS_LONG, 0, (long)((char *)&ftrfmode_example.bunchedBeamMode), NULL, 0.0, 1, "If non-zero, then do calculations bunch-by-bunch."},
    };

ZLONGIT zlongit_example;
/* ZLONGIT physical parameters */
PARAMETER zlongit_param[N_ZLONGIT_PARAMS] = {
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&zlongit_example.charge), NULL, 0.0, 0, "beam charge (or use CHARGE element)"},
    {"BROAD_BAND", "", IS_LONG, 0, (long)((char *)&zlongit_example.broad_band), NULL, 0.0, 0, "broad-band impedance?"},
    {"RA", "Ohm", IS_DOUBLE, 0, (long)((char *)&zlongit_example.Ra), NULL, 0.0, 0, "shunt impedance, Ra=V^2/P"},
    {"RS", "Ohm", IS_DOUBLE, 0, (long)((char *)&zlongit_example.Rs), NULL, 0.0, 0, "shunt impedance (Rs=Ra/2)"},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&zlongit_example.Q), NULL, 0.0, 1, "cavity Q"},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&zlongit_example.freq), NULL, 0.0, 0, "frequency (BROAD_BAND=1)"},
    {"ZREAL", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&zlongit_example.Zreal), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving real part of impedance vs f (BROAD_BAND=0)"},
    {"ZIMAG", "", IS_STRING, PARAM_XY_WAVEFORM, (long)((char *)&zlongit_example.Zimag), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving imaginary part of impedance vs f (BROAD_BAND=0)"},
    {"BIN_SIZE", "S", IS_DOUBLE, 0, (long)((char *)&zlongit_example.bin_size), NULL, 0.0, 0, "bin size for current histogram (use 0 for autosize)"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&zlongit_example.n_bins), NULL, 0.0, 128, "number of bins for current histogram"},
    {"MAX_N_BINS", "", IS_LONG, 0, (long)((char *)&zlongit_example.max_n_bins), NULL, 0.0, 0, "Maximum number of bins for current histogram"},
    {"WAKES", "", IS_STRING, 0, (long)((char *)&zlongit_example.wakes), NULL, 0.0, 0, "filename for output of wake"},
    {"WAKE_INTERVAL", "", IS_LONG, 0, (long)((char *)&zlongit_example.wake_interval), NULL, 0.0, 1, "interval in passes at which to output wake"},
    {"WAKE_START", "", IS_LONG, 0, (long)((char *)&zlongit_example.wake_start), NULL, 0.0, 0, "pass at which to start to output wake"},
    {"WAKE_END", "", IS_LONG, 0, (long)((char *)&zlongit_example.wake_end), NULL, 0.0, LONG_MAX, "pass at which to stop to output wake"},
    {"AREA_WEIGHT", "", IS_LONG, 0, (long)((char *)&zlongit_example.area_weight), NULL, 0.0, 0, "use area-weighting in assigning charge to histogram?"},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&zlongit_example.interpolate), NULL, 0.0, 0, "interpolate wake?"},
    {"SMOOTHING", "", IS_LONG, 0, (long)((char *)&zlongit_example.smoothing), NULL, 0.0, 0, "Use Savitzky-Golay filter to smooth current histogram?"},
    {"SG_ORDER", "", IS_LONG, 0, (long)((char *)&zlongit_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing"},
    {"SG_HALFWIDTH", "", IS_LONG, 0, (long)((char *)&zlongit_example.SGHalfWidth), NULL, 0.0, 4, "Savitzky-Golay filter halfwidth for smoothing"},
    {"REVERSE_TIME_ORDER", "", IS_LONG, 0, (long)((char *)&zlongit_example.reverseTimeOrder), NULL, 0.0, 0, "Reverse time-order of particles for wake computation?"},
    {"FACTOR", "", IS_DOUBLE, 0, (long)((char *)&zlongit_example.factor), NULL, 1.0, 0, "Factor by which to multiply impedance."},
    {"START_ON_PASS", "", IS_LONG, 0, (long)((char *)&zlongit_example.startOnPass), NULL, 1.0, 0, "The pass on which the impedance effects start."},
    {"RAMP_PASSES", "", IS_LONG, 0, (long)((char *)&zlongit_example.rampPasses), NULL, 0.0, 0, "Number of passes over which to linearly ramp up the impedance to full strength."},
    {"HIGH_FREQUENCY_CUTOFF0", "", IS_DOUBLE, 0, (long)((char*)&zlongit_example.highFrequencyCutoff0), NULL, -1.0, 0, "Frequency at which smoothing filter begins.  If not positive, no frequency filter smoothing is done.  Frequency is in units of Nyquist (0.5/binsize)."},
    {"HIGH_FREQUENCY_CUTOFF1", "", IS_DOUBLE, 0, (long)((char*)&zlongit_example.highFrequencyCutoff1), NULL, -1.0, 0, "Frequency at which smoothing filter is 0.  If not given, defaults to HIGH_FREQUENCY_CUTOFF0."},
    {"BUNCHED_BEAM_MODE", "", IS_LONG, 0, (long)((char *)&zlongit_example.bunchedBeamMode), NULL, 0.0, 1, "If non-zero, then do calculations bunch-by-bunch."},
    {"ALLOW_LONG_BEAM", "", IS_LONG, 0, (long)((char *)&zlongit_example.allowLongBeam), NULL, 0.0, 0, "Allow beam longer than covered by impedance data?"},
    };

SREFFECTS SReffects_example;
/* SREFFECTS physical parameters */
PARAMETER sreffects_param[N_SREFFECTS_PARAMS] = {
    {"JX", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.Jx), NULL, 1.0, 0, "x damping partition number"},
    {"JY", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.Jy), NULL, 1.0, 0, "y damping partition number"},
    {"JDELTA", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.Jdelta), NULL, 2.0, 0, "momentum damping partition number"},
    {"EXREF", "m", IS_DOUBLE, 0, (long)((char *)&SReffects_example.exRef), NULL, 0.0, 0, "reference equilibrium x emittance"},
    {"EYREF", "m", IS_DOUBLE, 0, (long)((char *)&SReffects_example.eyRef), NULL, 0.0, 0, "reference equilibrium y emittance"},
    {"SDELTAREF", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.SdeltaRef), NULL, 0.0, 0, "reference equilibrium fractional momentum spread"},
    {"DDELTAREF", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.DdeltaRef), NULL, 0.0, 0, "reference fractional momentum change per turn due to SR (negative value)"},
    {"PREF", "m$be$nc", IS_DOUBLE, 0, (long)((char *)&SReffects_example.pRef), NULL, 0.0, 0, "reference momentum (to which other reference values pertain)"},
    {"COUPLING", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.coupling), NULL, 0.0, 0, "x-y coupling"},
    {"FRACTION", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.fraction), NULL, 1.0, 0, "fraction of implied SR effect to simulate with each instance"},
    {"DAMPING", "", IS_LONG, 0, (long)((char *)&SReffects_example.damp), NULL, 0, 1, "include damping, less rf effects?"},
    {"QEXCITATION", "", IS_LONG, 0, (long)((char *)&SReffects_example.qExcite), NULL, 0, 1, "include quantum excitation?"},
    {"LOSSES", "", IS_LONG, 0, (long)((char *)&SReffects_example.loss), NULL, 0, 1, "include average losses?"},
    {"CUTOFF", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.cutoff), NULL, 100.0, 0, "cutoff (in sigmas) for gaussian random numbers"},
    {"INCLUDE_OFFSETS", "", IS_LONG, 0, (long)((char *)&SReffects_example.includeOffsets), NULL, 0, 1, "include orbit offsets in tracking (see below)?"},
    };

MODRF modrf_example;
PARAMETER modrf_param[N_MODRF_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&modrf_example.length), NULL, 0.0, 0, "length"},
    {"VOLT", "V", IS_DOUBLE, 0, (long)((char *)&modrf_example.volt), NULL, 0.0, 0, "nominal voltage"},
    {"PHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&modrf_example.phase), NULL, 0.0, 0, "nominal phase"},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&modrf_example.freq), NULL, 500.0e6, 0, "nominal frequency"},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&modrf_example.Q), NULL, 0.0, 0, "cavity Q"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&modrf_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"AMMAG", "", IS_DOUBLE, 0, (long)((char *)&modrf_example.amMag), NULL, 0.0, 0, "magnitude of amplitude modulation (fraction value)"},
    {"AMPHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&modrf_example.amPhase), NULL, 0.0, 0, "phase of amplitude modulation"},
    {"AMFREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&modrf_example.amFreq), NULL, 0.0, 0, "frequency of amplitude modulation"},
    {"AMDECAY", "1/s", IS_DOUBLE, 0, (long)((char *)&modrf_example.amDecay), NULL, 0.0, 0, "exponential decay rate of amplitude modulation"},
    {"PMMAG", "DEG", IS_DOUBLE, 0, (long)((char *)&modrf_example.pmMag), NULL, 0.0, 0, "magnitude of phase modulation"},
    {"PMPHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&modrf_example.pmPhase), NULL, 0.0, 0, "phase of phase modulation"},
    {"PMFREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&modrf_example.pmFreq), NULL, 0.0, 0, "frequency of phase modulation"},
    {"PMDECAY", "1/s", IS_DOUBLE, 0, (long)((char *)&modrf_example.pmDecay), NULL, 0.0, 0, "exponential decay rate of phase modulation"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&modrf_example.fiducial), NULL, 0.0, 0, "mode for determining fiducial arrival time (light, tmean, first, pmaximum)"},
    };    

BMAPXY bmapxy_example;
PARAMETER bmapxy_param[N_BMAPXY_PARAMS] = {
  {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bmapxy_example.length), NULL, 0.0, 0, "length"},
  {"STRENGTH", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bmapxy_example.strength), NULL, 0.0, 0, "factor by which to multiply field"},
  {"ACCURACY", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bmapxy_example.accuracy), NULL, 0.0, 0, "integration accuracy"},
  {"METHOD", NULL, IS_STRING, PARAM_CHANGES_MATRIX, (long)((char*)&bmapxy_example.method), NULL, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, modified-midpoint, two-pass modified-midpoint, leap-frog, non-adaptive runge-kutta"},
  {"FILENAME", NULL, IS_STRING, PARAM_CHANGES_MATRIX, (long)((char*)&bmapxy_example.filename), NULL, 0.0, 0, "name of file containing columns (x, y, Fx, Fy) giving normalized field (Fx, Fy) vs (x, y)"},
  {"FX", NULL, IS_STRING, PARAM_CHANGES_MATRIX, (long)((char*)&bmapxy_example.FxRpn), NULL, 0.0, 0, "rpn expression for Fx in terms of x and y"},
  {"FY", NULL, IS_STRING, PARAM_CHANGES_MATRIX, (long)((char*)&bmapxy_example.FyRpn), NULL, 0.0, 0, "rpn expression for Fy in terms of x and y"},
};  

ZTRANSVERSE ztransverse_example;
PARAMETER ztransverse_param[N_ZTRANSVERSE_PARAMS] = {
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.charge), NULL, 0.0, 0, "beam charge (or use CHARGE element)"},
    {"BROAD_BAND", "", IS_LONG, 0, (long)((char *)&ztransverse_example.broad_band), NULL, 0.0, 0, "broad-band impedance?"},
    {"RS", "Ohm/m", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.Rs), NULL, 0.0, 0, "shunt impedance (Rs=Ra/2=V^2/(2*P))"},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.Q), NULL, 0.0, 1, "cavity Q"},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.freq), NULL, 0.0, 0, "frequency (BROAD_BAND=1)"},
    {"INPUTFILE", "", IS_STRING, 0, (long)((char*)&ztransverse_example.inputFile), NULL, 0.0, 0, "name of file giving impedance (BROAD_BAND=0)"},
    {"FREQCOLUMN", "", IS_STRING, 0, (long)((char*)&ztransverse_example.freqColumn), NULL, 0.0, 0, "column in INPUTFILE containing frequency"},
    {"ZXREAL", "", IS_STRING, 0, (long)((char *)&ztransverse_example.ZxReal), NULL, 0.0, 0, "column in INPUTFILE containing real impedance for x plane"},
    {"ZXIMAG", "", IS_STRING, 0, (long)((char *)&ztransverse_example.ZxImag), NULL, 0.0, 0, "column in INPUTFILE containing imaginary impedance for x plane"},
    {"ZYREAL", "", IS_STRING, 0, (long)((char *)&ztransverse_example.ZyReal), NULL, 0.0, 0, "column in INPUTFILE containing real impedance for y plane"},
    {"ZYIMAG", "", IS_STRING, 0, (long)((char *)&ztransverse_example.ZyImag), NULL, 0.0, 0, "column in INPUTFILE containing imaginary impedance for y plane"},
    {"BIN_SIZE", "S", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.bin_size), NULL, 0.0, 0, "bin size for current histogram (use 0 for autosize)"},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&ztransverse_example.interpolate), NULL, 0.0, 0, "interpolate wake?"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&ztransverse_example.n_bins), NULL, 0.0, 128, "number of bins for current histogram"},
    {"MAX_N_BINS", "", IS_LONG, 0, (long)((char *)&ztransverse_example.max_n_bins), NULL, 0.0, 0, "Maximum number of bins for current histogram"},
    {"SMOOTHING", "", IS_LONG, 0, (long)((char *)&ztransverse_example.smoothing), NULL, 0.0, 0, "Use Savitzky-Golay filter to smooth current histogram?"},
    {"SG_ORDER", "", IS_LONG, 0, (long)((char *)&ztransverse_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing"},
    {"SG_HALFWIDTH", "", IS_LONG, 0, (long)((char *)&ztransverse_example.SGHalfWidth), NULL, 0.0, 4, "Savitzky-Golay filter halfwidth for smoothing"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.dy), NULL, 0.0, 0, "misalignment"},
    {"FACTOR", "", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.factor), NULL, 1.0, 0, "Factor by which to multiply x and y impedances."},
    {"XFACTOR", "", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.xfactor), NULL, 1.0, 0, "Factor by which to multiply x impedance."},
    {"YFACTOR", "", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.yfactor), NULL, 1.0, 0, "Factor by which to multiply y impedance."},
    {"WAKES", "", IS_STRING, 0, (long)((char *)&ztransverse_example.wakes), NULL, 0.0, 0, "filename for output of wake"},
    {"WAKE_INTERVAL", "", IS_LONG, 0, (long)((char *)&ztransverse_example.wake_interval), NULL, 0.0, 1, "interval in passes at which to output wake"},
    {"WAKE_START", "", IS_LONG, 0, (long)((char *)&ztransverse_example.wake_start), NULL, 0.0, 0, "pass at which to start to output wake"},
    {"WAKE_END", "", IS_LONG, 0, (long)((char *)&ztransverse_example.wake_end), NULL, 0.0, LONG_MAX, "pass at which to stop to output wake"},
    {"START_ON_PASS", "", IS_LONG, 0, (long)((char *)&ztransverse_example.startOnPass), NULL, 1.0, 0, "The pass on which the impedance effects start."},
    {"RAMP_PASSES", "", IS_LONG, 0, (long)((char *)&ztransverse_example.rampPasses), NULL, 1.0, 0, "Number of passes over which to linearly ramp up the impedance to full strength."},
    {"HIGH_FREQUENCY_CUTOFF0", "", IS_DOUBLE, 0, (long)((char*)&ztransverse_example.highFrequencyCutoff0), NULL, -1.0, 0, "Frequency at which smoothing filter begins.  If not positive, no frequency filter smoothing is done.  Frequency is in units of Nyquist (0.5/binsize)."},
    {"HIGH_FREQUENCY_CUTOFF1", "", IS_DOUBLE, 0, (long)((char*)&ztransverse_example.highFrequencyCutoff1), NULL, -1.0, 0, "Frequency at which smoothing filter is 0.  If not given, defaults to HIGH_FREQUENCY_CUTOFF0."},
    {"X_DRIVE_EXPONENT", "", IS_LONG, 0, (long)((char *)&ztransverse_example.xDriveExponent), NULL, 0.0, 1, "Exponent applied to x coordinates of drive particles"},
    {"Y_DRIVE_EXPONENT", "", IS_LONG, 0, (long)((char *)&ztransverse_example.yDriveExponent), NULL, 0.0, 1, "Exponent applied to y coordinates of drive particles"},
    {"X_PROBE_EXPONENT", "", IS_LONG, 0, (long)((char *)&ztransverse_example.xProbeExponent), NULL, 0.0, 0, "Exponent applied to x coordinates of probe particles"},
    {"Y_PROBE_EXPONENT", "", IS_LONG, 0, (long)((char *)&ztransverse_example.yProbeExponent), NULL, 0.0, 0, "Exponent applied to y coordinates of probe particles"},
    {"BUNCHED_BEAM_MODE", "", IS_LONG, 0, (long)((char *)&ztransverse_example.bunchedBeamMode), NULL, 0.0, 1, "If non-zero, then do calculations bunch-by-bunch."},
    {"ALLOW_LONG_BEAM", "", IS_LONG, 0, (long)((char *)&ztransverse_example.allowLongBeam), NULL, 0.0, 0, "Allow beam longer than covered by impedance data?"},
};

IBSCATTER ibs_example;
PARAMETER ibscatter_param[N_IBSCATTER_PARAMS] = {
  {"FACTOR", "", IS_DOUBLE, 0, (long)((char *)&ibs_example.factor), NULL, 1.0, 0, "factor by which to multiply growth rates before using"},
  {"DO_X", "", IS_LONG, 0, (long)((char *)&ibs_example.do_x), NULL, 0.0, 1, "do x-plane scattering?"},
  {"DO_Y", "", IS_LONG, 0, (long)((char *)&ibs_example.do_y), NULL, 0.0, 1, "do y-plane scattering?"},
  {"DO_Z", "", IS_LONG, 0, (long)((char *)&ibs_example.do_z), NULL, 0.0, 1, "do z-plane scattering?"},
  {"NSLICE", "", IS_LONG, 0, (long)((char *)&ibs_example.nslice), NULL, 0.0, 1, "The number of slices per bunch"},
  {"SMOOTH", "", IS_LONG, 0, (long)((char *)&ibs_example.smooth), NULL, 0.0, 1, "Use smooth method instead of random numbers?"},
  {"FORCE_MATCHED_TWISS", "", IS_LONG, 0, (long)((char *)&ibs_example.forceMatchedTwiss), NULL, 0.0, 0, "Force computations to be done with twiss parameters of the beamline, not the beam."},
  {"ISRING", "", IS_LONG, 0, (long)((char *)&ibs_example.isRing), NULL, 0.0, 1, "Is it storage ring?"},
  {"INTERVAL", "", IS_LONG, 0, (long)((char *)&ibs_example.interval), NULL, 0.0, 1, "Interval in passes at which to update output file."},
  {"FILENAME", "", IS_STRING, 0, (long)((char *)&ibs_example.filename), NULL, 0.0, 0, "Output filename."},
  {"BUNCHED_BEAM_MODE", "", IS_LONG, 0, (long)((char *)&ibs_example.bunchedBeamMode), NULL, 0.0, 1, "If non-zero, then do calculations bunch-by-bunch."},
  {"VERBOSE", "", IS_LONG, 0, (long)((char *)&ibs_example.verbose), NULL, 0.0, 0, "If non-zero, then print updates during calculations."},
};

WAKE wake_example;
/* WAKE physical parameters */
PARAMETER wake_param[N_WAKE_PARAMS] = {
    {"INPUTFILE", "", IS_STRING, 0, (long)((char *)&wake_example.inputFile), NULL, 0.0, 0, "name of file giving Green function"},
    {"TCOLUMN", "", IS_STRING, 0, (long)((char *)&wake_example.tColumn), NULL, 0.0, 0, "column in INPUTFILE containing time data"},
    {"WCOLUMN", "", IS_STRING, 0, (long)((char *)&wake_example.WColumn), NULL, 0.0, 0, "column in INPUTFILE containing Green function"},
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&wake_example.charge), NULL, 0.0, 0, "beam charge (or use CHARGE element)"},
    {"FACTOR", "C", IS_DOUBLE, 0, (long)((char *)&wake_example.factor), NULL, 1.0, 0, "factor by which to multiply wake"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&wake_example.n_bins), NULL, 0.0, 0, "number of bins for current histogram"},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&wake_example.interpolate), NULL, 0.0, 0, "interpolate wake?"},
    {"SMOOTHING", "", IS_LONG, 0, (long)((char *)&wake_example.smoothing), NULL, 0.0, 0, "Use Savitzky-Golay filter to smooth current histogram?"},
    {"SG_HALFWIDTH", "", IS_LONG, 0, (long)((char *)&wake_example.SGHalfWidth), NULL, 0.0, 4, "Savitzky-Golay filter half-width for smoothing"},
    {"SG_ORDER", "", IS_LONG, 0, (long)((char *)&wake_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing"},
    {"CHANGE_P0", "", IS_LONG, 0, (long)((char *)&wake_example.change_p0), NULL, 0.0, 0, "change central momentum?"},
    {"ALLOW_LONG_BEAM", "", IS_LONG, 0, (long)((char *)&wake_example.allowLongBeam), NULL, 0.0, 0, "allow beam longer than wake data?"},
    {"RAMP_PASSES", "", IS_LONG, 0, (long)((char *)&wake_example.rampPasses), NULL, 0.0, 0, "Number of passes over which to linearly ramp up the wake to full strength."},
    {"BUNCHED_BEAM_MODE", "", IS_LONG, 0, (long)((char *)&wake_example.bunchedBeamMode), NULL, 0.0, 1, "If non-zero, then do calculations bunch-by-bunch."},
    {"ACAUSAL_ALLOWED", "", IS_LONG, 0, (long)((char *)&wake_example.acausalAllowed), NULL, 0.0, 0, "If non-zero, then an acausal wake is allowed."},
    };

CORGPIPE corgpipe_example;
/* CORGPIPE physical parameters */
PARAMETER corgpipe_param[N_CORGPIPE_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&corgpipe_example.length), NULL, 0.0, 0, "length"},
    {"RADIUS", "M", IS_DOUBLE, 0, (long)((char *)&corgpipe_example.radius), NULL, 0.0, 0, "pipe radius"},
    {"PERIOD", "M", IS_DOUBLE, 0, (long)((char *)&corgpipe_example.period), NULL, 0.0, 0, "period of corrugations (<< radius recommended)"},
    {"GAP", "M", IS_DOUBLE, 0, (long)((char *)&corgpipe_example.gap), NULL, 0.0, 0, "gap in corrugations (< period required)"},
    {"DEPTH", "M", IS_DOUBLE, 0, (long)((char *)&corgpipe_example.depth), NULL, 0.0, 0, "depth of corrugations (<< radius, >~ period recommended)"},
    {"DT", "S", IS_DOUBLE, 0, (long)((char *)&corgpipe_example.dt), NULL, 0.0, 0, "maximum time duration of wake (0 for autoscale)"},
    {"TMAX", "S", IS_DOUBLE, 0, (long)((char *)&corgpipe_example.tmax), NULL, 0.0, 0, "maximum time duration of wake (0 for autoscale)"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&corgpipe_example.n_bins), NULL, 0.0, 0, "number of bins for charge histogram (0 for autoscale)"},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&corgpipe_example.interpolate), NULL, 0.0, 0, "interpolate wake?"},
    {"SMOOTHING", "", IS_LONG, 0, (long)((char *)&corgpipe_example.smoothing), NULL, 0.0, 0, "Use Savitzky-Golay filter to smooth current histogram?"},
    {"SG_HALFWIDTH", "", IS_LONG, 0, (long)((char *)&corgpipe_example.SGHalfWidth), NULL, 0.0, 4, "Savitzky-Golay filter half-width for smoothing"},
    {"SG_ORDER", "", IS_LONG, 0, (long)((char *)&corgpipe_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing"},
    {"CHANGE_P0", "", IS_LONG, 0, (long)((char *)&corgpipe_example.change_p0), NULL, 0.0, 0, "change central momentum?"},
    {"ALLOW_LONG_BEAM", "", IS_LONG, 0, (long)((char *)&corgpipe_example.allowLongBeam), NULL, 0.0, 0, "allow beam longer than wake data?"},
    {"RAMP_PASSES", "", IS_LONG, 0, (long)((char *)&corgpipe_example.rampPasses), NULL, 0.0, 0, "Number of passes over which to linearly ramp up the wake to full strength."},
    };

LRWAKE lrwake_example;
/* LRWAKE physical parameters */
PARAMETER lrwake_param[N_LRWAKE_PARAMS] = {
    {"INPUTFILE", "", IS_STRING, 0, (long)((char *)&lrwake_example.inputFile), NULL, 0.0, 0, "name of file giving Green function"},
    {"TCOLUMN", "", IS_STRING, 0, (long)((char *)&lrwake_example.WColumn[0]), NULL, 0.0, 0, "column in INPUTFILE containing time data"},
    {"WXCOLUMN", "", IS_STRING, 0, (long)((char *)&lrwake_example.WColumn[1]), NULL, 0.0, 0, "column in INPUTFILE containing horizontal dipole Green function"},
    {"WYCOLUMN", "", IS_STRING, 0, (long)((char *)&lrwake_example.WColumn[2]), NULL, 0.0, 0, "column in INPUTFILE containing vertical dipole Green function"},
    {"WZCOLUMN", "", IS_STRING, 0, (long)((char *)&lrwake_example.WColumn[3]), NULL, 0.0, 0, "column in INPUTFILE containing longitudinal Green function"},
    {"QXCOLUMN", "", IS_STRING, 0, (long)((char *)&lrwake_example.WColumn[4]), NULL, 0.0, 0, "column in INPUTFILE containing horizontal quadrupole Green function"},
    {"QYCOLUMN", "", IS_STRING, 0, (long)((char *)&lrwake_example.WColumn[5]), NULL, 0.0, 0, "column in INPUTFILE containing vertical quadrupole Green function"},
    {"FACTOR", "", IS_DOUBLE, 0, (long)((char *)&lrwake_example.factor), NULL, 1.0, 0, "factor by which to multiply wakes"},
    {"XFACTOR", "", IS_DOUBLE, 0, (long)((char *)&lrwake_example.xFactor), NULL, 1.0, 0, "factor by which to multiply longitudinal wake"},
    {"YFACTOR", "", IS_DOUBLE, 0, (long)((char *)&lrwake_example.yFactor), NULL, 1.0, 0, "factor by which to multiply horizontal dipole wake"},
    {"ZFACTOR", "", IS_DOUBLE, 0, (long)((char *)&lrwake_example.zFactor), NULL, 1.0, 0, "factor by which to multiply vertical dipole wake"},
    {"QXFACTOR", "", IS_DOUBLE, 0, (long)((char *)&lrwake_example.qxFactor), NULL, 1.0, 0, "factor by which to multiply horizontal quadrupole wake"},
    {"QYFACTOR", "", IS_DOUBLE, 0, (long)((char *)&lrwake_example.qyFactor), NULL, 1.0, 0, "factor by which to multiply vertical quadrupole wake"},
    {"TURNS_TO_KEEP", "", IS_LONG, 0, (long)((char *)&lrwake_example.turnsToKeep), NULL, 0.0, 128, "number of turns of data to retain"},
    {"RAMP_PASSES", "", IS_LONG, 0, (long)((char *)&lrwake_example.rampPasses), NULL, 0.0, 0, "Number of passes over which to linearly ramp up the wake to full strength."},
};

TRWAKE trwake_example;
/* TRWAKE physical parameters */
PARAMETER trwake_param[N_TRWAKE_PARAMS] = {
    {"INPUTFILE", "", IS_STRING, 0, (long)((char *)&trwake_example.inputFile), NULL, 0.0, 0, "name of file giving Green functions"},
    {"TCOLUMN", "", IS_STRING, 0, (long)((char *)&trwake_example.tColumn), NULL, 0.0, 0, "column in INPUTFILE containing time data"},
    {"WXCOLUMN", "", IS_STRING, 0, (long)((char *)&trwake_example.WxColumn), NULL, 0.0, 0, "column in INPUTFILE containing x Green function"},
    {"WYCOLUMN", "", IS_STRING, 0, (long)((char *)&trwake_example.WyColumn), NULL, 0.0, 0, "column in INPUTFILE containing y Green function"},
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&trwake_example.charge), NULL, 0.0, 0, "beam charge (or use CHARGE element)"},
    {"FACTOR", "", IS_DOUBLE, 0, (long)((char *)&trwake_example.factor), NULL, 1.0, 0, "factor by which to multiply both wakes"},
    {"XFACTOR", "", IS_DOUBLE, 0, (long)((char *)&trwake_example.xfactor), NULL, 1.0, 0, "factor by which to multiply x wake"},
    {"YFACTOR", "", IS_DOUBLE, 0, (long)((char *)&trwake_example.yfactor), NULL, 1.0, 0, "factor by which to multiply y wake"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&trwake_example.n_bins), NULL, 0.0, 0, "number of bins for current histogram"},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&trwake_example.interpolate), NULL, 0.0, 0, "interpolate wake?"},
    {"SMOOTHING", "", IS_LONG, 0, (long)((char *)&trwake_example.smoothing), NULL, 0.0, 0, "Use Savitzky-Golay filter to smooth current histogram?"},
    {"SG_HALFWIDTH", "", IS_LONG, 0, (long)((char *)&trwake_example.SGHalfWidth), NULL, 0.0, 4, "Savitzky-Golay filter half-width for smoothing"},
    {"SG_ORDER", "", IS_LONG, 0, (long)((char *)&trwake_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&trwake_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&trwake_example.dy), NULL, 0.0, 0, "misalignment"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&trwake_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"X_DRIVE_EXPONENT", "", IS_LONG, 0, (long)((char *)&trwake_example.xDriveExponent), NULL, 0.0, 1, "Exponent applied to x coordinates of drive particles"},
    {"Y_DRIVE_EXPONENT", "", IS_LONG, 0, (long)((char *)&trwake_example.yDriveExponent), NULL, 0.0, 1, "Exponent applied to y coordinates of drive particles"},
    {"X_PROBE_EXPONENT", "", IS_LONG, 0, (long)((char *)&trwake_example.xProbeExponent), NULL, 0.0, 0, "Exponent applied to x coordinates of probe particles"},
    {"Y_PROBE_EXPONENT", "", IS_LONG, 0, (long)((char *)&trwake_example.yProbeExponent), NULL, 0.0, 0, "Exponent applied to y coordinates of probe particles"},
    {"RAMP_PASSES", "", IS_LONG, 0, (long)((char *)&trwake_example.rampPasses), NULL, 0.0, 0, "Number of passes over which to linearly ramp up the wake to full strength."},
    {"BUNCHED_BEAM_MODE", "", IS_LONG, 0, (long)((char *)&trwake_example.bunchedBeamMode), NULL, 0.0, 1, "If non-zero, then do calculations bunch-by-bunch."},
    {"ACAUSAL_ALLOWED", "", IS_LONG, 0, (long)((char *)&trwake_example.acausalAllowed), NULL, 0.0, 0, "If non-zero, then an acausal wake is allowed."},
    };

CHARGE charge_example;
/* CHARGE physical parameters */
PARAMETER charge_param[N_CHARGE_PARAMS] = {
    {"TOTAL", "C", IS_DOUBLE, 0, (long)((char *)&charge_example.charge), NULL, 0.0, 0, "total charge in beam"},
    {"PER_PARTICLE", "C", IS_DOUBLE, 0, (long)((char *)&charge_example.chargePerParticle), NULL, 0.0, 0, "charge per macroparticle"},
    {"ALLOW_TOTAL_CHANGE", NULL, IS_LONG, 0, (long)((char *)&charge_example.allowChangeWhileRunning), NULL, 0.0, 0, "If nonzero, allow total charge to change while tracking even if number of particles does not change.  Useful for ramping of charge."},
};

PFILTER pfilter_example;
/* PFILTER physical parameters */
PARAMETER pfilter_param[N_PFILTER_PARAMS] = {
    {"DELTALIMIT", "", IS_DOUBLE, 0, (long)((char *)&pfilter_example.deltaLimit), NULL, -1.0, 0, "maximum fractional momentum deviation"},
    {"LOWERFRACTION", "", IS_DOUBLE, 0, (long)((char *)&pfilter_example.lowerFraction), NULL, 0.0, 0, "fraction of lowest-momentum particles to remove"},
    {"UPPERFRACTION", "", IS_DOUBLE, 0, (long)((char *)&pfilter_example.upperFraction), NULL, 0.0, 0, "fraction of highest-momentum particles to remove"},
    {"FIXPLIMITS", "", IS_LONG, 0, (long)((char *)&pfilter_example.fixPLimits), NULL, 0.0, 0, "fix the limits in p from LOWERFRACTION and UPPERFRACTION applied to first beam"},
    {"BEAMCENTERED", "", IS_LONG, 0, (long)((char *)&pfilter_example.beamCentered), NULL, 0.0, 0, "if nonzero, center for DELTALIMIT is average beam momentum"},
    {"BINS", "", IS_LONG, 0, (long)((char*)&pfilter_example.bins), NULL, 0.0, 1024, "number of bins"},
};

HISTOGRAM histogram_example;
/* HISTOGRAM physical parameters */
PARAMETER histogram_param[N_HISTOGRAM_PARAMS] = {
  {"FILENAME", "", IS_STRING, 0, (long)((char *)&histogram_example.filename), "", 0.0, 0, "filename for histogram output, possibly incomplete (see below)"},
  {"INTERVAL", "", IS_LONG, 0, (long)((char *)&histogram_example.interval), NULL, 0.0, 1, "interval in passes between output"},
  {"START_PASS", "", IS_LONG, 0, (long)((char*)&histogram_example.startPass), NULL, 0.0, 0, "starting pass for output"},
  {"BINS", "", IS_LONG, 0, (long)((char*)&histogram_example.bins), NULL, 0.0, 50, "number of bins"},
  {"FIXED_BIN_SIZE", "", IS_SHORT, 0, (long)((char*)&histogram_example.fixedBinSize), NULL, 0.0, 0, "if nonzero, bin size is fixed after the first histogram is made"},
  {"X_DATA", "", IS_SHORT, 0, (long)((char*)&histogram_example.xData), NULL, 0.0, 1, "histogram x and x'?"},
  {"Y_DATA", "", IS_SHORT, 0, (long)((char*)&histogram_example.yData), NULL, 0.0, 1, "histogram y and y'?"},
  {"LONGIT_DATA", "", IS_SHORT, 0, (long)((char*)&histogram_example.longitData), NULL, 0.0, 1, "histogram t and p?"},
  {"BIN_SIZE_FACTOR", "", IS_DOUBLE, 0, (long)((char*)&histogram_example.binSizeFactor), NULL, 1.0, 0, "multiply computed bin size by this factor before histogramming"},
  {"NORMALIZE", "", IS_SHORT, 0, (long)((char*)&histogram_example.normalize), NULL, 0.0, 1, "normalize histogram with bin size and number of particles?"},
  {"DISABLE", "", IS_SHORT, 0, (long)((char *)&histogram_example.disable), NULL, 0.0, 0, "If nonzero, no output will be generated."},    
  {"SPARSE", "", IS_SHORT, 0, (long)((char *)&histogram_example.sparse), NULL, 0.0, 0, "If nonzero, only bins with non-zero counts will be output."},    
  {"START_PID", "", IS_LONG, 0, (long)((char *)&histogram_example.startPID), NULL, 0.0, -1, "starting particleID for particles to include"},
  {"END_PID", "", IS_LONG, 0, (long)((char *)&histogram_example.endPID), NULL, 0.0, -1, "ending particleID for particles to include"},
};

MHISTOGRAM mhistogram_example;
/* MHISTOGRAM physical parameters */
PARAMETER mhistogram_param[N_MHISTOGRAM_PARAMS] = {
  {"FILE1D", "", IS_STRING, 0, (long)((char *)&mhistogram_example.file1d), NULL, 0.0, 0, "filename for 1d histogram output, possibly incomplete (see below)"},
  {"FILE2DH", "", IS_STRING, 0, (long)((char *)&mhistogram_example.file2dH), NULL, 0.0, 0, "filename for 2d x-x' histogram output, possibly incomplete (see below)"},
  {"FILE2DV", "", IS_STRING, 0, (long)((char *)&mhistogram_example.file2dV), NULL, 0.0, 0, "filename for 2d y-y' histogram output, possibly incomplete (see below)"},
  {"FILE2DL", "", IS_STRING, 0, (long)((char *)&mhistogram_example.file2dL), NULL, 0.0, 0, "filename for 2d dt-deltaP histogram output, possibly incomplete (see below)"},
  {"FILE4D", "", IS_STRING, 0, (long)((char *)&mhistogram_example.file4d), NULL, 0.0, 0, "filename for 4d x-x'-y-y' histogram output, possibly incomplete (see below)"},
  {"FILE6D", "", IS_STRING, 0, (long)((char *)&mhistogram_example.file6d), NULL, 0.0, 0, "filename for 6d x-x'-y-y'-dt-deltaP histogram output, possibly incomplete (see below)"},
  {"INPUT_BINS", " ", IS_STRING, 0, (long)((char *)&mhistogram_example.inputBinFile), NULL, 0.0, 0, "Name of SDDS file contains input bin number."},
  {"INTERVAL", "", IS_LONG, 0, (long)((char *)&mhistogram_example.interval), NULL, 0.0, 1, "interval in passes between output."},
  {"START_PASS", "", IS_LONG, 0, (long)((char*)&mhistogram_example.startPass), NULL, 0.0, 0, "starting pass for output"},
  {"NORMALIZE", "", IS_SHORT, 0, (long)((char*)&mhistogram_example.normalize), NULL, 0.0, 1, "normalize histogram with number of particles?"},
  {"DISABLE", "", IS_SHORT, 0, (long)((char *)&mhistogram_example.disable), NULL, 0.0, 0, "If nonzero, no output will be generated."},    
  {"LUMPED", "", IS_SHORT, 0, (long)((char *)&mhistogram_example.lumped), NULL, 0.0, 0, "If nonzero, then results at elements with same name will be output to a single multipage SDDS file."},    
};
  
CSRDRIFT csrdrift_example;
/* CSR drift length physical parameters */
PARAMETER csrdrift_param[N_CSRDRIFT_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&csrdrift_example.length), NULL, 0.0, 0, "length"},
    {"ATTENUATION_LENGTH", "M", IS_DOUBLE, 0, (long)((char *)&csrdrift_example.attenuationLength), NULL, 0.0, 0, "exponential attenuation length for wake"},
    {"DZ", "", IS_DOUBLE, 0, (long)((char *)&csrdrift_example.dz), NULL, 0.0, 0, "interval between kicks"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&csrdrift_example.nKicks), NULL, 0.0, 1, "number of kicks (if DZ is zero)"},
    {"SPREAD", "", IS_SHORT, 0, (long)((char *)&csrdrift_example.spread), NULL, 0.0, 0, "use spreading function?"},
    {"USE_OVERTAKING_LENGTH", "", IS_SHORT, 0, (long)((char *)&csrdrift_example.useOvertakingLength), NULL, 0.0, 0, "use overtaking length for ATTENUATION_LENGTH?"},
    {"OL_MULTIPLIER", "", IS_DOUBLE, 0, (long)((char *)&csrdrift_example.overtakingLengthMultiplier), NULL, 1.0, 0, "factor by which to multiply the overtaking length to get the attenuation length"},
    {"CSR", "", IS_SHORT, 0, (long)((char *)&csrdrift_example.csr), NULL, 0.0, 1, "do CSR calcuations"},
    {"USE_SALDIN54", "", IS_SHORT, 0, (long)((char*)&csrdrift_example.useSaldin54), NULL, 0, 0, "Use Saldin et al eq. 54 (NIM A 398 (1997) 373-394 for decay vs z?"},
    {"SALDIN54POINTS", "", IS_LONG, 0, (long)((char*)&csrdrift_example.nSaldin54Points), NULL, 0.0, 1000, "Number of values of position inside bunch to average for Saldin eq 54."},
    {"SALDIN54NORM_MODE", "", IS_STRING, 0, (long)((char *)&csrdrift_example.normMode), "peak", 0.0, 0, "peak or first"},
    {"SPREAD_MODE", "", IS_STRING, 0, (long)((char *)&csrdrift_example.spreadMode), "full", 0.0, 0, "full, simple, or radiation-only"},
    {"WAVELENGTH_MODE", "", IS_STRING, 0, (long)((char *)&csrdrift_example.wavelengthMode), "sigmaz", 0.0, 0, "sigmaz or peak-to-peak"},
    {"BUNCHLENGTH_MODE", "", IS_STRING, 0, (long)((char *)&csrdrift_example.bunchlengthMode), "68-percentile", 0.0, 0, "rms, 68-percentile, or 90-percentile"},
    {"SALDIN54_OUTPUT", "", IS_STRING, 0, (long)((char*)&csrdrift_example.Saldin54Output), NULL, 0.0, 0, "Filename for output of CSR intensity vs. z as computed using Saldin eq 54."},
    {"USE_STUPAKOV", "", IS_SHORT, 0, (long)((char *)&csrdrift_example.useStupakov), NULL, 0.0, 0, "Use treatment from G. Stupakov's note of 9/12/2001?"},
    {"STUPAKOV_OUTPUT", "", IS_STRING, 0, (long)((char*)&csrdrift_example.StupakovOutput), NULL, 0.0, 0, "Filename for output of CSR wake vs. s as computed using Stupakov's equations."},
    {"STUPAKOV_OUTPUT_INTERVAL", "", IS_LONG, 0, (long)((char*)&csrdrift_example.StupakovOutputInterval), NULL, 0.0, 1, "Interval (in kicks) between output of Stupakov wakes."},
    {"SLICE_ANALYSIS_INTERVAL", "", IS_LONG, 0, (long)((char*)&csrdrift_example.sliceAnalysisInterval), NULL, 0.0, 0, "interval (in kicks) of output to slice analysis file (from slice_analysis command)"},
    {"LINEARIZE", "", IS_SHORT, 0, (long)((char*)&csrdrift_example.linearOptics), NULL, 0.0, 0, "use linear optics for drift pieces?"},
    {"LSC_INTERPOLATE", "", IS_SHORT, 0, (long)((char *)&csrdrift_example.LSCInterpolate), NULL, 0.0, 1, "Interpolate computed LSC wake?"},
    {"LSC_BINS", "", IS_LONG, 0, (long)((char*)&csrdrift_example.LSCBins), NULL, 0.0, 0, "If non-zero, include LSC with given number of bins."},
    {"LSC_LOW_FREQUENCY_CUTOFF0", "", IS_DOUBLE, 0, (long)((char*)&csrdrift_example.LSCLowFrequencyCutoff0), NULL, -1.0, 0, "Highest spatial frequency at which low-frequency cutoff filter is zero.  If not positive, no low-frequency cutoff filter is applied. Frequency is in units of Nyquist (0.5/binsize)."},
    {"LSC_LOW_FREQUENCY_CUTOFF1", "", IS_DOUBLE, 0, (long)((char*)&csrdrift_example.LSCLowFrequencyCutoff1), NULL, -1.0, 0, "Lowest spatial frequency at which low-frequency cutoff filter is 1.  If not given, defaults to LOW_FREQUENCY_CUTOFF1."},
    {"LSC_HIGH_FREQUENCY_CUTOFF0", "", IS_DOUBLE, 0, (long)((char*)&csrdrift_example.LSCHighFrequencyCutoff0), NULL, -1.0, 0, "Spatial frequency at which smoothing filter begins for LSC.  If not positive, no frequency filter smoothing is done.  Frequency is in units of Nyquist (0.5/binsize)."},
    {"LSC_HIGH_FREQUENCY_CUTOFF1", "", IS_DOUBLE, 0, (long)((char*)&csrdrift_example.LSCHighFrequencyCutoff1), NULL, -1.0, 0, "Spatial frequency at which smoothing filter is 0 for LSC.  If not given, defaults to HIGH_FREQUENCY_CUTOFF0."},
    {"LSC_RADIUS_FACTOR", "", IS_DOUBLE, 0, (long)((char*)&csrdrift_example.LSCRadiusFactor), NULL, 1.7, 0, "Radius factor for LSC computation."},
    };

RFCW rfcw_example;
/* rf cavity with wakes physical parameters */
PARAMETER rfcw_param[N_RFCW_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfcw_example.length), NULL, 0.0, 0, "length"},
    {"CELL_LENGTH", "M", IS_DOUBLE, 0, (long)((char *)&rfcw_example.cellLength), NULL, 0.0, 0, "cell length (used to scale wakes, which are assumed to be given for a cell, according to L/CELL_LENGTH)"},
    {"VOLT", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfcw_example.volt), NULL, 0.0, 0, "voltage"},
    {"PHASE", "DEG", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfcw_example.phase), NULL, 0.0, 0, "phase"},
    {"FREQ", "Hz", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfcw_example.freq), NULL, 500.0e6, 0, "frequency"},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&rfcw_example.Q), NULL, 0.0, 0, "cavity Q (for cavity that charges up to voltage from 0)"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&rfcw_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"CHANGE_P0", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&rfcw_example.change_p0), NULL, 0.0, 0, "does element change central momentum?"}, 
    {"CHANGE_T", "", IS_LONG, 0, (long)((char *)&rfcw_example.change_t), NULL, 0.0, 0, "see RFCA documentation"}, 
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&rfcw_example.fiducial), NULL, 0.0, 0, "mode for determining fiducial arrival time (light, tmean, first, pmaximum)"},
    {"END1_FOCUS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&rfcw_example.end1Focus), NULL, 0.0, 0, "include focusing at entrance?"},
    {"END2_FOCUS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&rfcw_example.end2Focus), NULL, 0.0, 0, "include focusing at exit?"},
    {"BODY_FOCUS_MODEL", "", IS_STRING, PARAM_CHANGES_MATRIX, (long)((char *)&rfcw_example.bodyFocusModel), NULL, 0.0, 0, "None (default) or SRS (simplified Rosenzweig/Serafini for standing wave)"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&rfcw_example.nKicks), NULL, 0.0, 0, "Number of kicks to use for kick method.  Set to zero for matrix method."},
    {"ZWAKE", "", IS_LONG, 0, (long)((char *)&rfcw_example.includeZWake), NULL, 0.0, 1, "If zero, longitudinal wake is turned off."},
    {"TRWAKE", "", IS_LONG, 0, (long)((char *)&rfcw_example.includeTrWake), NULL, 0.0, 1, "If zero, transverse wakes are turned off."},
    {"WAKEFILE", "", IS_STRING, 0, (long)((char *)&rfcw_example.wakeFile), NULL, 0.0, 0, "name of file containing Green functions"},
    {"ZWAKEFILE", "", IS_STRING, 0, (long)((char *)&rfcw_example.zWakeFile), NULL, 0.0, 0, "if WAKEFILE=NULL, optional name of file containing longitudinal Green function"},
    {"TRWAKEFILE", "", IS_STRING, 0, (long)((char *)&rfcw_example.trWakeFile), NULL, 0.0, 0, "if WAKEFILE=NULL, optional name of file containing transverse Green functions"},
    {"TCOLUMN", "", IS_STRING, 0, (long)((char *)&rfcw_example.tColumn), NULL, 0.0, 0, "column containing time data"},
    {"WXCOLUMN", "", IS_STRING, 0, (long)((char *)&rfcw_example.WxColumn), NULL, 0.0, 0, "column containing x Green function"},
    {"WYCOLUMN", "", IS_STRING, 0, (long)((char *)&rfcw_example.WyColumn), NULL, 0.0, 0, "column containing y Green function"},
    {"WZCOLUMN", "", IS_STRING, 0, (long)((char *)&rfcw_example.WzColumn), NULL, 0.0, 0, "column containing longitudinal Green function"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&rfcw_example.n_bins), NULL, 0.0, 0, "number of bins for current histogram"},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&rfcw_example.interpolate), NULL, 0.0, 0, "interpolate wake?"},
    {"SMOOTHING", "", IS_LONG, 0, (long)((char *)&rfcw_example.smoothing), NULL, 0.0, 0, "Use Savitzky-Golay filter to smooth current histogram?"},
    {"SG_HALFWIDTH", "", IS_LONG, 0, (long)((char *)&rfcw_example.SGHalfWidth), NULL, 0.0, 4, "Savitzky-Golay filter half-width for smoothing"},
    {"SG_ORDER", "", IS_LONG, 0, (long)((char *)&rfcw_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfcw_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&rfcw_example.dy), NULL, 0.0, 0, "misalignment"},
    {"LINEARIZE", "", IS_LONG, 0, (long)((char *)&rfcw_example.linearize), NULL, 0.0, 0, "Linearize phase dependence?"},
    {"LSC", "", IS_LONG, 0, (long)((char *)&rfcw_example.doLSC), NULL, 0.0, 0, "Include longitudinal space-charge impedance?"},
    {"LSC_BINS", "", IS_LONG, 0, (long)((char *)&rfcw_example.LSCBins), NULL, 0.0, 1024, "Number of bins for LSC calculations"},
    {"LSC_INTERPOLATE", "", IS_LONG, 0, (long)((char *)&rfcw_example.LSCInterpolate), NULL, 0.0, 1, "Interpolate computed LSC wake?"},
    {"LSC_LOW_FREQUENCY_CUTOFF0", "", IS_DOUBLE, 0, (long)((char*)&rfcw_example.LSCLowFrequencyCutoff0), NULL, -1.0, 0, "Highest spatial frequency at which low-frequency cutoff filter is zero.  If not positive, no low-frequency cutoff filter is applied. Frequency is in units of Nyquist (0.5/binsize)."},
    {"LSC_LOW_FREQUENCY_CUTOFF1", "", IS_DOUBLE, 0, (long)((char*)&rfcw_example.LSCLowFrequencyCutoff1), NULL, -1.0, 0, "Lowest spatial frequency at which low-frequency cutoff filter is 1.  If not given, defaults to LOW_FREQUENCY_CUTOFF1."},
    {"LSC_HIGH_FREQUENCY_CUTOFF0", "", IS_DOUBLE, 0, (long)((char*)&rfcw_example.LSCHighFrequencyCutoff0), NULL, -1.0, 0, "Spatial frequency at which smoothing filter begins for LSC.  If not positive, no frequency filter smoothing is done.  Frequency is in units of Nyquist (0.5/binsize)."},
    {"LSC_HIGH_FREQUENCY_CUTOFF1", "", IS_DOUBLE, 0, (long)((char*)&rfcw_example.LSCHighFrequencyCutoff1), NULL, -1.0, 0, "Spatial frequency at which smoothing filter is 0 for LSC.  If not given, defaults to HIGH_FREQUENCY_CUTOFF0."},
    {"LSC_RADIUS_FACTOR", "", IS_DOUBLE, 0, (long)((char*)&rfcw_example.LSCRadiusFactor), NULL, 1.7, 0, "LSC radius is (Sx+Sy)/2*RADIUS_FACTOR"},
    {"WAKES_AT_END", "", IS_LONG, 0, (long)((char *)&rfcw_example.wakesAtEnd), NULL, 0.0, 0, "Do wake kicks at end of segment (for backward compatibility)?"},
    };
   
REMCOR remcor_example;
/* beam centering physical parameters */
PARAMETER remcor_param[N_REMCOR_PARAMS]={
    {"X" , "", IS_SHORT, 0, (long)((char *)&remcor_example.x), NULL, 0.0, 1, "remove correlations in x?"},
    {"XP", "", IS_SHORT, 0, (long)((char *)&remcor_example.xp), NULL, 0.0, 1, "remove correlations in x'?"},
    {"Y" , "", IS_SHORT, 0, (long)((char *)&remcor_example.y), NULL, 0.0, 1, "remove correlations in y?"},
    {"YP", "", IS_SHORT, 0, (long)((char *)&remcor_example.yp), NULL, 0.0, 1, "remove correlations in y'?"},
    {"WITH", "", IS_SHORT, 0, (long)((char *)&remcor_example.with), NULL, 0.0, 6, "coordinate to remove correlations with (1,2,3,4,5,6)=(x,x',y,y',s,dP/Po)"},
    {"ONCE_ONLY", "", IS_SHORT, 0, (long)((char *)&remcor_example.onceOnly), NULL, 0.0, 0, "compute correction only for first beam, apply to all?"},
    };

MAP_SOLENOID mapSol_example;

PARAMETER mapSolenoid_param[N_MAPSOLENOID_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mapSol_example.length), NULL, 0.0, 0, "length"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mapSol_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mapSol_example.dy), NULL, 0.0, 0, "misalignment"},
    {"ETILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mapSol_example.eTilt), NULL, 0.0, 0, "misalignment"},
    {"EYAW", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mapSol_example.eYaw), NULL, 0.0, 0, "misalignment"},
    {"EPITCH", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mapSol_example.ePitch), NULL, 0.0, 0, "misalignment"},
    {"N_STEPS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&mapSol_example.n_steps), NULL, 0.0, 100, "number of steps (for nonadaptive integration)"},
    {"INPUTFILE", "", IS_STRING, 0, (long)((char *)&mapSol_example.inputFile), NULL, 0.0, 0, "SDDS file containing (Br, Bz) vs (r, z).  Each page should have values for a fixed r."},
    {"RCOLUMN", "", IS_STRING, 0, (long)((char *)&mapSol_example.rColumn), NULL, 0.0, 0, "column containing r values"},
    {"ZCOLUMN", "", IS_STRING, 0, (long)((char *)&mapSol_example.zColumn), NULL, 0.0, 0, "column containing z values"},
    {"BRCOLUMN", "", IS_STRING, 0, (long)((char *)&mapSol_example.BrColumn), NULL, 0.0, 0, "column containing Br values"},
    {"BZCOLUMN", "", IS_STRING, 0, (long)((char *)&mapSol_example.BzColumn), NULL, 0.0, 0, "column containing Bz values"},
    {"FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mapSol_example.factor), NULL, DEFAULT_ACCURACY, 0, "factor by which to multiply fields in file"},
    {"BXUNIFORM", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mapSol_example.BxUniform), NULL, 0.0, 0, "uniform horizontal field to superimpose on solenoid field"},
    {"BYUNIFORM", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mapSol_example.ByUniform), NULL, 0.0, 0, "uniform vertical field to superimpose on solenoid field"},
    {"LUNIFORM", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mapSol_example.lUniform), NULL, 0.0, 0, "length of uniform field superimposed on solenoid field"},
    {"ACCURACY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mapSol_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&mapSol_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, non-adaptive runge-kutta, modified midpoint)"},
    } ;

TWISSELEMENT twissElem_example;

PARAMETER twissElement_param[N_TWISSELEMENT_PARAMS] = {
  {"BETAX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss.betax), NULL, 1.0, 0, "horizontal beta function"},
  {"ALPHAX", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss.alphax), NULL, 0.0, 0, "horizontal alpha function"},
  {"ETAX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss.etax), NULL, 0.0, 0, "horizontal eta function"},
  {"ETAXP", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss.etapx), NULL, 0.0, 0, "slope of horizontal eta function"},
  {"BETAY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss.betay), NULL, 1.0, 0, "vertical beta function"},
  {"ALPHAY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss.alphay), NULL, 0.0, 0, "vertical alpha function"},
  {"ETAY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss.etay), NULL, 0.0, 0, "vertical eta function"},
  {"ETAYP", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss.etapy), NULL, 0.0, 0, "slope of vertical eta function"},
  {"FROM_BEAM", "", IS_SHORT, 0, (long)((char *)&twissElem_example.fromBeam), NULL, 0.0, 0, "compute transformation from tracked beam properties instead of Twiss parameters?"},
  {"FROM_0VALUES", "", IS_SHORT, 0, (long)((char *)&twissElem_example.from0Values), NULL, 0.0, 0, "if non-zero, transformation is from the \"0\" values provided in the element definition"},
  {"COMPUTE_ONCE", "", IS_SHORT, 0, (long)((char *)&twissElem_example.computeOnce), NULL, 0.0, 0, "compute transformation only for first beam or lattice functions?"},
  {"APPLY_ONCE", "", IS_SHORT, 0, (long)((char *)&twissElem_example.applyOnce), NULL, 0.0, 1, "apply correction only on first pass through for each beam?"},
  {"VERBOSE", "", IS_SHORT, 0, (long)((char *)&twissElem_example.verbose), NULL, 0.0, 0, "if non-zero, print extra information about transformations"},
  {"DISABLE", "", IS_SHORT, 0, (long)((char *)&twissElem_example.disable), NULL, 0.0, 0, "if non-zero, element is ignored"},
  {"BETAX0", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss0.betax), NULL, 1.0, 0, "initial horizontal beta function (if FROM_0VALUES nonzero)"},
  {"ALPHAX0", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss0.alphax), NULL, 0.0, 0, "initial horizontal alpha function (if FROM_0VALUES nonzero)"},
  {"ETAX0", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss0.etax), NULL, 0.0, 0, "initial horizontal eta function (if FROM_0VALUES nonzero)"},
  {"ETAXP0", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss0.etapx), NULL, 0.0, 0, "initial slope of horizontal eta function (if FROM_0VALUES nonzero)"},
  {"BETAY0", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss0.betay), NULL, 1.0, 0, "initial vertical beta function (if FROM_0VALUES nonzero)"},
  {"ALPHAY0", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss0.alphay), NULL, 0.0, 0, "initial vertical alpha function (if FROM_0VALUES nonzero)"},
  {"ETAY0", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss0.etay), NULL, 0.0, 0, "initial vertical eta function (if FROM_0VALUES nonzero)"},
  {"ETAYP0", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&twissElem_example.twiss0.etapy), NULL, 0.0, 0, "initial slope of vertical eta function (if FROM_0VALUES nonzero)"},
};

WIGGLER wiggler_example;

PARAMETER wiggler_param[N_WIGGLER_PARAMS] = {
  {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&wiggler_example.length), NULL, 0.0, 0, "length"},
  {"RADIUS", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&wiggler_example.radius), NULL, 0.0, 0, "Peak bending radius.  Ignored if K or B is non-negative."},
  {"K", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&wiggler_example.K), NULL, 0.0, 0, "Dimensionless strength parameter."},
  {"B", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&wiggler_example.B), NULL, 0.0, 0, "Peak vertical magnetic field. Ignored if K is non-negative"},
  {"DX", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&wiggler_example.dx), NULL, 0.0, 0, "Misaligment."},
  {"DY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&wiggler_example.dy), NULL, 0.0, 0, "Misaligment."},
  {"DZ", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&wiggler_example.dz), NULL, 0.0, 0, "Misaligment."},
  {"TILT", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&wiggler_example.tilt), NULL, 0.0, 0, "Rotation about beam axis."},
  {"POLES", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&wiggler_example.poles), NULL, 0.0, 0, "Number of wiggler poles"},
  {"FOCUSING", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&wiggler_example.focusing), NULL, 0.0, 1, "If 0, turn off vertical focusing (this is unphysical!)"},
} ;

CWIGGLER cwiggler_example;

PARAMETER cwiggler_param[N_CWIGGLER_PARAMS] = {
  {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.length), NULL, 0.0, 0, "Total length"},
  {"B_MAX", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.BMax), NULL, 0.0, 0, "Maximum on-axis magnetic field."},
  {"BX_MAX", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.BxMax), NULL, 0.0, 0, "Maximum on-axis magnetic field.  Ignored if B_MAX is nonzero."},
  {"BY_MAX", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.ByMax), NULL, 0.0, 0, "Maximum on-axis magnetic field.  Ignored if B_MAX is nonzero."},
  {"TGU_GRADIENT", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.tguGradient), NULL, 0.0, 0, "Transverse gradient divided by maximum on-axis field, used if TGU=1."},
  {"TGU_COMP_FACTOR", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.tguCompFactor), NULL, 1.0, 0, "Use to adjust constant field component to reduce trajectory error."},
  {"POLE1_FACTOR", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.poleFactor[0]), NULL, 1.0, 0, "Use to adjust first and last pole strength, e.g., to reduce trajectory error."},
  {"POLE2_FACTOR", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.poleFactor[1]), NULL, 1.0, 0, "Use to adjust second and penultimate pole strength, e.g., to reduce trajectory error."},
  {"POLE3_FACTOR", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.poleFactor[2]), NULL, 1.0, 0, "Use to adjust third and third-from=last pole strength, e.g., to reduce trajectory error."},
  {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.dx), NULL, 0.0, 0, "Misaligment."},
  {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.dy), NULL, 0.0, 0, "Misaligment."},
  {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.dz), NULL, 0.0, 0, "Misaligment."},
  {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.tilt), NULL, 0.0, 0, "Rotation about beam axis."},
  {"PERIODS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.periods), NULL, 0.0, 0, "Number of wiggler periods."},
  {"STEPS_PER_PERIOD", "", IS_LONG, 0, (long)((char *)&cwiggler_example.stepsPerPeriod), NULL, 0.0, 10, "Integration steps per period."},
  {"INTEGRATION_ORDER", "", IS_SHORT, 0, (long)((char *)&cwiggler_example.integrationOrder), NULL, 0.0, 4, "Integration order (2 or 4)."},
  {"BY_FILE", " ", IS_STRING, 0, (long)((char *)&cwiggler_example.ByFile), NULL, 0.0, 0, "Name of SDDS file with By harmonic data."},
  {"BX_FILE", " ", IS_STRING, 0, (long)((char *)&cwiggler_example.BxFile), NULL, 0.0, 0, "Name of SDDS file with Bx harmonic data."},
  {"BY_SPLIT_POLE", "", IS_SHORT, 0, (long)((char *)&cwiggler_example.BySplitPole), NULL, 0.0, 0, "Use \"split-pole\" expansion for By?"},
  {"BX_SPLIT_POLE", "", IS_SHORT, 0, (long)((char *)&cwiggler_example.BxSplitPole), NULL, 0.0, 0, "Use \"split-pole\" expansion for Bx?"},
  {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&cwiggler_example.sr), NULL, 0.0, 0, "Include classical, single-particle synchrotron radiation?"},
  {"ISR", "", IS_SHORT, 0, (long)((char *)&cwiggler_example.isr), NULL, 0.0, 0, "Include incoherent synchrotron radiation (quantum excitation)?"},
  {"ISR1PART", "", IS_SHORT, 0, (long)((char *)&cwiggler_example.isr1Particle), NULL, 0.0, 1, "Include ISR for single-particle beam only if ISR=1 and ISR1PART=1"},
  {"SINUSOIDAL", "", IS_SHORT, 0, (long)((char *)&cwiggler_example.sinusoidal), NULL, 0.0, 0, "Ideal sinusoidal wiggler?  If non-zero, BX_FILE and BY_FILE are not used."},
  {"VERTICAL", "", IS_SHORT,  PARAM_CHANGES_MATRIX, (long)((char *)&cwiggler_example.vertical), NULL, 0.0, 0, "If SINUSOIDAL is non-zero, then setting this to non-zero gives a vertical wiggler.  Default is horizontal."},
  {"HELICAL", "", IS_SHORT, 0, (long)((char *)&cwiggler_example.helical), NULL, 0.0, 0, "Ideal helical wiggler?  If non-zero and SINUSOIDAL is also non-zero, BX_FILE and BY_FILE are not used."},
  {"TGU", "", IS_SHORT, 0, (long)((char *)&cwiggler_example.tgu), NULL, 0.0, 0, "Ideal transverse gradient undulator? If non-zero and SINUSOIDAL is also non-zero, BX_FILE and BY_FILE are not used. Give gradient in TGU_GRADIENT."},
  {"FORCE_MATCHED", "", IS_SHORT, 0, (long)((char *)&cwiggler_example.forceMatched), NULL, 0.0, 1, "Force matched dispersion for first harmonics?  If non-zero, start and end of magnetic field will be inset from the ends of the device if phase is not 0 or $\\pi$."},
  {"FIELD_OUTPUT", "", IS_STRING, 0, (long)((char *)&cwiggler_example.fieldOutput), NULL, 0.0, 0, "Name of file to which field samples will be written.  Slow, so use only for debugging."},
  {"VERBOSITY", "", IS_SHORT, 0, (long)((char *)&cwiggler_example.verbosity), NULL, 0.0, 0, "A higher value requires more detailed printouts related to computations."},
} ;

APPLE apple_example;

PARAMETER apple_param[N_APPLE_PARAMS] = {
  {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.length), NULL, 0.0, 0, "Total length"},
  {"B_MAX", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.BMax), NULL, 0.0, 0, "Maximum on-axis magnetic field at gap=GAP0 and equal longitudinal phases of PHASE_1,2,3,4"},
  {"SHIM_SCALE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.shimScale), NULL, 1.0, 0, "Scaling factor of shim correction field."},
  {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.dx), NULL, 0.0, 0, "Misaligment."},
  {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.dy), NULL, 0.0, 0, "Misaligment."},
  {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.dz), NULL, 0.0, 0, "Misaligment."},
  {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.tilt), NULL, 0.0, 0, "Rotation about beam axis."},
  {"PERIODS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.periods), NULL, 0.0, 0, "Total number of wiggler periods. Include end poles"},
  {"STEP", "", IS_LONG, 0, (long)((char *)&apple_example.step), NULL, 0.0, 1, "Number of normal periods to track for each step"},
  {"ORDER", "", IS_SHORT, 0, (long)((char *)&apple_example.order), NULL, 0.0, 0, "Order=3 including the 3rd order terms. Otherwise using 2nd order formula."},
  {"END_POLE", "", IS_SHORT, 0, (long)((char *)&apple_example.End_Pole), NULL, 0.0, 1, "The ending poles are treated as 2 half periods at each sides of the wiggler with reducing field strength, such as 0.25, -0.75, ..., 0.75, -0.25. Periods has to > 2"},
  {"SHIM_ON", "", IS_SHORT, 0, (long)((char *)&apple_example.shimOn), NULL, 0.0, 0, "Include shim correction"},
  {"INPUT_FILE", " ", IS_STRING, 0, (long)((char *)&apple_example.Input), NULL, 0.0, 0, "Name of SDDS file with By harmonic data given at GAP0 and equal longitudinal phases."},
  {"SHIM_INPUT", " ", IS_STRING, 0, (long)((char *)&apple_example.shimInput), NULL, 0.0, 0, "Name of SDDS file with shim field integral harmonic data given at GAP0."},
  {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&apple_example.sr), NULL, 0.0, 0, "Include classical, single-particle synchrotron radiation?"},
  {"ISR", "", IS_SHORT, 0, (long)((char *)&apple_example.isr), NULL, 0.0, 0, "Include incoherent synchrotron radiation (quantum excitation)?"},
  {"ISR1PART", "", IS_SHORT, 0, (long)((char *)&apple_example.isr1Particle), NULL, 0.0, 1, "Include ISR for single-particle beam only if ISR=1 and ISR1PART=1"},
  {"X0", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.x0), NULL, 0.0, 0, "Offset of magnet row center in meter."},
  {"GAP0", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.gap0), NULL, 0.0, 0, "Nominal magnetic gap."},
  {"D_GAP", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.dgap), NULL, 0.0, 0, "Delta gap: actual gap - nominal gap"},
  {"PHASE_1", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.phi1), NULL, 0.0, 0, "Longitudinal phase of the first row (top right)"},
  {"PHASE_2", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.phi2), NULL, 0.0, 0, "Longitudinal phase of the second row (top left)"},
  {"PHASE_3", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.phi3), NULL, 0.0, 0, "Longitudinal phase of the third row (bottom left)"},
  {"PHASE_4", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&apple_example.phi4), NULL, 0.0, 0, "Longitudinal phase of the fourth row (bottom right) "},
  {"VERBOSITY", "", IS_SHORT, 0, (long)((char *)&apple_example.verbosity), NULL, 0.0, 0, "A higher value requires more detailed printouts related to computations."},
} ;

SCRIPT script_example;

PARAMETER script_param[N_SCRIPT_PARAMS] = {
  {"L", "M", IS_DOUBLE, 0, (long)((char *)&script_example.length), NULL, 0.0, 0, "Length to be used for matrix-based operations such as twiss parameter computation."},
  {"COMMAND", "", IS_STRING, 0, (long)((char *)&script_example.command), NULL, 0.0, 0, "SDDS-compliant command to apply to the beam.  Use the sequence %i to represent the input filename and %o to represent the output filename."},
  {"USE_CSH", "", IS_SHORT, 0, (long)((char *)&script_example.useCsh), NULL, 0.0, 1, "Use C-shell for execution (may be slower)?"},
  {"VERBOSITY", "", IS_SHORT, 0, (long)((char *)&script_example.verbosity), NULL, 0.0, 0, "Set the verbosity level."},
  {"RPN_PARAMETERS", "", IS_SHORT, 0, (long)((char *)&script_example.rpnParameters), NULL, 0.0, 0, "If nonzero, then parameters from the script output file are loaded into RPN variables."},
  {"START_PASS", "", IS_LONG, 0, (long)((char *)&script_example.startPass), NULL, 0.0, -1, "Start script action on this pass.  Before that, behaves like a drift space."},
  {"END_PASS", "", IS_LONG, 0, (long)((char *)&script_example.endPass), NULL, 0.0, -1, "End script action after this pass.  Before that, behaves like a drift space."},
  {"PASS_INTERVAL", "", IS_LONG, 0, (long)((char *)&script_example.passInterval), NULL, 0.0, -1, "Execute script only every Nth pass following START_PASS, including START_PASS. Otherwise, behaves like a drift space."},
  {"ON_PASS", "", IS_LONG, 0, (long)((char *)&script_example.onPass), NULL, 0.0, -1, "Perform script action only on this pass, overriding other pass controls. Other than that, behaves like a drift space."},
  {"DIRECTORY", "", IS_STRING, 0, (long)((char *)&script_example.directory), NULL, 0.0, 0, "Directory in which to place input and output files.  If blank, the present working directory is used."},
  {"ROOTNAME", "", IS_STRING, 0, (long)((char *)&script_example.rootname), NULL, 0.0, 0, "Rootname for use in naming input and output files.  %s may be used to represent the run rootname."},
  {"INPUT_EXTENSION", "", IS_STRING, 0, (long)((char *)&script_example.inputExtension), "in", 0.0, 0, "Extension for the script input file."},
  {"OUTPUT_EXTENSION", "", IS_STRING, 0, (long)((char *)&script_example.outputExtension), "out", 0.0, 0, "Extension for the script output file."},
  {"KEEP_FILES", "", IS_SHORT, 0, (long)((char *)&script_example.keepFiles), NULL, 0.0, 0, "If nonzero, then script input and output files are not deleted after use.  By default, they are deleted."},
  {"DRIFT_MATRIX", "", IS_SHORT, 0, (long)((char *)&script_example.driftMatrix), NULL, 0.0, 0, "If nonzero, then for non-tracking calculations the element is treated as a drift space."},
  {"USE_PARTICLE_ID", "", IS_SHORT, 0, (long)((char *)&script_example.useParticleID), NULL, 0.0, 1, "If nonzero, then the output file will supply particle IDs. Otherwise, particles are renumbered."},
  {"NO_NEW_PARTICLES", "", IS_SHORT, 0, (long)((char *)&script_example.noNewParticles), NULL, 0.0, 1, "If nonzero, then no new particles will be added in the script output file."},
  {"DETERMINE_LOSSES_FROM_PID", "", IS_SHORT, 0, (long)((char *)&script_example.determineLossesFromParticleID), NULL, 0.0, 1, "If nonzero and if USE_PARTICLE_ID is nonzero, then particleID data from script output is used to determine which particles were lost."},
  {"SOFT_FAILURE", "", IS_SHORT, 0, (long)((char *)&script_example.softFailure), NULL, 0.0, 1, "If output file does not exist or can't be read, consider all particles lost."},
  {"NP0", "", IS_DOUBLE, 0, (long)((char *)&script_example.NP[0]), NULL, 0.0, 0, "User-defined numerical parameter for command substitution for sequence %np0"},
  {"NP1", "", IS_DOUBLE, 0, (long)((char *)&script_example.NP[1]), NULL, 0.0, 0, "User-defined numerical parameter for command substitution for sequence %np1"},
  {"NP2", "", IS_DOUBLE, 0, (long)((char *)&script_example.NP[2]), NULL, 0.0, 0, "User-defined numerical parameter for command substitution for sequence %np2"},
  {"NP3", "", IS_DOUBLE, 0, (long)((char *)&script_example.NP[3]), NULL, 0.0, 0, "User-defined numerical parameter for command substitution for sequence %np3"},
  {"NP4", "", IS_DOUBLE, 0, (long)((char *)&script_example.NP[4]), NULL, 0.0, 0, "User-defined numerical parameter for command substitution for sequence %np4"},
  {"NP5", "", IS_DOUBLE, 0, (long)((char *)&script_example.NP[5]), NULL, 0.0, 0, "User-defined numerical parameter for command substitution for sequence %np5"},
  {"NP6", "", IS_DOUBLE, 0, (long)((char *)&script_example.NP[6]), NULL, 0.0, 0, "User-defined numerical parameter for command substitution for sequence %np6"},
  {"NP7", "", IS_DOUBLE, 0, (long)((char *)&script_example.NP[7]), NULL, 0.0, 0, "User-defined numerical parameter for command substitution for sequence %np7"},
  {"NP8", "", IS_DOUBLE, 0, (long)((char *)&script_example.NP[8]), NULL, 0.0, 0, "User-defined numerical parameter for command substitution for sequence %np8"},
  {"NP9", "", IS_DOUBLE, 0, (long)((char *)&script_example.NP[9]), NULL, 0.0, 0, "User-defined numerical parameter for command substitution for sequence %np9"},
  {"SP0", "", IS_STRING, 0, (long)((char *)&script_example.SP[0]), NULL, 0.0, 0, "User-defined string parameter for command substitution for sequence %sp0"},
  {"SP1", "", IS_STRING, 0, (long)((char *)&script_example.SP[1]), NULL, 0.0, 0, "User-defined string parameter for command substitution for sequence %sp1"},
  {"SP2", "", IS_STRING, 0, (long)((char *)&script_example.SP[2]), NULL, 0.0, 0, "User-defined string parameter for command substitution for sequence %sp2"},
  {"SP3", "", IS_STRING, 0, (long)((char *)&script_example.SP[3]), NULL, 0.0, 0, "User-defined string parameter for command substitution for sequence %sp3"},
  {"SP4", "", IS_STRING, 0, (long)((char *)&script_example.SP[4]), NULL, 0.0, 0, "User-defined string parameter for command substitution for sequence %sp4"},
  {"SP5", "", IS_STRING, 0, (long)((char *)&script_example.SP[5]), NULL, 0.0, 0, "User-defined string parameter for command substitution for sequence %sp5"},
  {"SP6", "", IS_STRING, 0, (long)((char *)&script_example.SP[6]), NULL, 0.0, 0, "User-defined string parameter for command substitution for sequence %sp6"},
  {"SP7", "", IS_STRING, 0, (long)((char *)&script_example.SP[7]), NULL, 0.0, 0, "User-defined string parameter for command substitution for sequence %sp7"},
  {"SP8", "", IS_STRING, 0, (long)((char *)&script_example.SP[8]), NULL, 0.0, 0, "User-defined string parameter for command substitution for sequence %sp8"},
  {"SP9", "", IS_STRING, 0, (long)((char *)&script_example.SP[9]), NULL, 0.0, 0, "User-defined string parameter for command substitution for sequence %sp9"},
};
  
FLOORELEMENT floorElem_example;

PARAMETER floor_param[N_FLOORELEMENT_PARAMS] = {
  {"X", "", IS_DOUBLE, 0, (long)((char *)&floorElem_example.position[0]), NULL, 0.0, 0, "X coordinate"},
  {"Y", "", IS_DOUBLE, 0, (long)((char *)&floorElem_example.position[1]), NULL, 0.0, 0, "Y coordinate"},
  {"Z", "", IS_DOUBLE, 0, (long)((char *)&floorElem_example.position[2]), NULL, 0.0, 0, "Z coordinate"},
  {"THETA", "", IS_DOUBLE, 0, (long)((char *)&floorElem_example.angle[0]), NULL, 0.0, 0, "theta value"},
  {"PHI", "", IS_DOUBLE, 0, (long)((char *)&floorElem_example.angle[1]), NULL, 0.0, 0, "phi value"},
  {"PSI", "", IS_DOUBLE, 0, (long)((char *)&floorElem_example.angle[2]), NULL, 0.0, 0, "psi value"},
};

LTHINLENS lthinlens_example;
PARAMETER lthinlens_param[N_LTHINLENS_PARAMS]={
    {"FX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lthinlens_example.fx), NULL, 0.0, 0, "focal length in horizontal plane"},
    {"FY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lthinlens_example.fy), NULL, 0.0, 0, "focal length in vertical plane"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lthinlens_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lthinlens_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lthinlens_example.dz), NULL, 0.0, 0, "misalignment"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lthinlens_example.tilt), NULL, 0.0, 0, "misalignment rotation about longitudinal axis"},
    {"YAW", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lthinlens_example.yaw), NULL, 0.0, 0, "misalignment rotation about vertical axis"},
    {"PITCH", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lthinlens_example.pitch), NULL, 0.0, 0, "misalignment rotation about transverse horizontal axis"},
    };

LMIRROR lmirror_example;
PARAMETER lmirror_param[N_LMIRROR_PARAMS]={
    {"RX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lmirror_example.Rx), NULL, 0.0, 0, "radius in horizontal plane"},
    {"RY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lmirror_example.Ry), NULL, 0.0, 0, "radius in vertical plane"},
    {"THETA", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lmirror_example.theta), NULL, 0.0, 0, "angle of incidence (in horizontal plane)"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lmirror_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lmirror_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lmirror_example.dz), NULL, 0.0, 0, "misalignment"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lmirror_example.tilt), NULL, 0.0, 0, "misalignment rotation about longitudinal axis"},
    {"YAW", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lmirror_example.yaw), NULL, 0.0, 0, "misalignment rotation about vertical axis"},
    {"PITCH", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lmirror_example.pitch), NULL, 0.0, 0, "misalignment rotation about transverse horizontal axis"},
    };

EMATRIX ematrix_example;
PARAMETER ematrix_param[N_EMATRIX_PARAMS]={
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&ematrix_example.length), NULL, 0.0, 0, "Length (used only for position computation)"},
    {"ANGLE", "RAD", IS_DOUBLE, 0, (long)((char *)&ematrix_example.angle), NULL, 0.0, 0, "Angle (used only for position computation)"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.dz), NULL, 0.0, 0, "misalignment"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.tilt), NULL, 0.0, 0, "Tilt angle"},
    {"YAW", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.yaw), NULL, 0.0, 0, "Yaw angle"},
    {"PITCH", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.pitch), NULL, 0.0, 0, "Pitch angle"},
    {"ORDER", "", IS_SHORT, 0, (long)((char *)&ematrix_example.order), NULL, 0.0, 0, ""},
    {"C1", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.C[0]), NULL, 0.0, 0, ""},
    {"C2", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.C[1]), NULL, 0.0, 0, ""},
    {"C3", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.C[2]), NULL, 0.0, 0, ""},
    {"C4", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.C[3]), NULL, 0.0, 0, ""},
    {"C5", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.C[4]), NULL, 0.0, 0, ""},
    {"C6", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.C[5]), NULL, 0.0, 0, "Change in momentum offset"},
    {"DELTAP", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.deltaP), NULL, 0.0, 0, "Change in central momentum (beta*gamma)"},
    {"R11", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[0][0]), NULL, 0.0, 0, ""},
    {"R12", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[0][1]), NULL, 0.0, 0, ""},
    {"R13", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[0][2]), NULL, 0.0, 0, ""},
    {"R14", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[0][3]), NULL, 0.0, 0, ""},
    {"R15", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[0][4]), NULL, 0.0, 0, ""},
    {"R16", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[0][5]), NULL, 0.0, 0, ""},
    {"R21", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[1][0]), NULL, 0.0, 0, ""},
    {"R22", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[1][1]), NULL, 0.0, 0, ""},
    {"R23", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[1][2]), NULL, 0.0, 0, ""},
    {"R24", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[1][3]), NULL, 0.0, 0, ""},
    {"R25", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[1][4]), NULL, 0.0, 0, ""},
    {"R26", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[1][5]), NULL, 0.0, 0, ""},
    {"R31", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[2][0]), NULL, 0.0, 0, ""},
    {"R32", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[2][1]), NULL, 0.0, 0, ""},
    {"R33", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[2][2]), NULL, 0.0, 0, ""},
    {"R34", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[2][3]), NULL, 0.0, 0, ""},
    {"R35", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[2][4]), NULL, 0.0, 0, ""},
    {"R36", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[2][5]), NULL, 0.0, 0, ""},
    {"R41", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[3][0]), NULL, 0.0, 0, ""},
    {"R42", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[3][1]), NULL, 0.0, 0, ""},
    {"R43", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[3][2]), NULL, 0.0, 0, ""},
    {"R44", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[3][3]), NULL, 0.0, 0, ""},
    {"R45", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[3][4]), NULL, 0.0, 0, ""},
    {"R46", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[3][5]), NULL, 0.0, 0, ""},
    {"R51", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[4][0]), NULL, 0.0, 0, ""},
    {"R52", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[4][1]), NULL, 0.0, 0, ""},
    {"R53", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[4][2]), NULL, 0.0, 0, ""},
    {"R54", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[4][3]), NULL, 0.0, 0, ""},
    {"R55", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[4][4]), NULL, 0.0, 0, ""},
    {"R56", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[4][5]), NULL, 0.0, 0, ""},
    {"R61", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[5][0]), NULL, 0.0, 0, ""},
    {"R62", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[5][1]), NULL, 0.0, 0, ""},
    {"R63", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[5][2]), NULL, 0.0, 0, ""},
    {"R64", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[5][3]), NULL, 0.0, 0, ""},
    {"R65", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[5][4]), NULL, 0.0, 0, ""},
    {"R66", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.R[5][5]), NULL, 0.0, 0, ""},
    {"T111", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][0][0]), NULL, 0.0, 0, ""},
    {"T121", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][1][0]), NULL, 0.0, 0, ""},
    {"T122", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][1][1]), NULL, 0.0, 0, ""},
    {"T131", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][2][0]), NULL, 0.0, 0, ""},
    {"T132", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][2][1]), NULL, 0.0, 0, ""},
    {"T133", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][2][2]), NULL, 0.0, 0, ""},
    {"T141", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][3][0]), NULL, 0.0, 0, ""},
    {"T142", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][3][1]), NULL, 0.0, 0, ""},
    {"T143", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][3][2]), NULL, 0.0, 0, ""},
    {"T144", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][3][3]), NULL, 0.0, 0, ""},
    {"T151", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][4][0]), NULL, 0.0, 0, ""},
    {"T152", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][4][1]), NULL, 0.0, 0, ""},
    {"T153", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][4][2]), NULL, 0.0, 0, ""},
    {"T154", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][4][3]), NULL, 0.0, 0, ""},
    {"T155", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][4][4]), NULL, 0.0, 0, ""},
    {"T161", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][5][0]), NULL, 0.0, 0, ""},
    {"T162", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][5][1]), NULL, 0.0, 0, ""},
    {"T163", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][5][2]), NULL, 0.0, 0, ""},
    {"T164", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][5][3]), NULL, 0.0, 0, ""},
    {"T165", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][5][4]), NULL, 0.0, 0, ""},
    {"T166", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[0][5][5]), NULL, 0.0, 0, ""},
    {"T211", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][0][0]), NULL, 0.0, 0, ""},
    {"T221", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][1][0]), NULL, 0.0, 0, ""},
    {"T222", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][1][1]), NULL, 0.0, 0, ""},
    {"T231", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][2][0]), NULL, 0.0, 0, ""},
    {"T232", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][2][1]), NULL, 0.0, 0, ""},
    {"T233", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][2][2]), NULL, 0.0, 0, ""},
    {"T241", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][3][0]), NULL, 0.0, 0, ""},
    {"T242", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][3][1]), NULL, 0.0, 0, ""},
    {"T243", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][3][2]), NULL, 0.0, 0, ""},
    {"T244", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][3][3]), NULL, 0.0, 0, ""},
    {"T251", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][4][0]), NULL, 0.0, 0, ""},
    {"T252", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][4][1]), NULL, 0.0, 0, ""},
    {"T253", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][4][2]), NULL, 0.0, 0, ""},
    {"T254", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][4][3]), NULL, 0.0, 0, ""},
    {"T255", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][4][4]), NULL, 0.0, 0, ""},
    {"T261", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][5][0]), NULL, 0.0, 0, ""},
    {"T262", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][5][1]), NULL, 0.0, 0, ""},
    {"T263", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][5][2]), NULL, 0.0, 0, ""},
    {"T264", "1", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][5][3]), NULL, 0.0, 0, ""},
    {"T265", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][5][4]), NULL, 0.0, 0, ""},
    {"T266", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[1][5][5]), NULL, 0.0, 0, ""},
    {"T311", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][0][0]), NULL, 0.0, 0, ""},
    {"T321", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][1][0]), NULL, 0.0, 0, ""},
    {"T322", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][1][1]), NULL, 0.0, 0, ""},
    {"T331", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][2][0]), NULL, 0.0, 0, ""},
    {"T332", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][2][1]), NULL, 0.0, 0, ""},
    {"T333", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][2][2]), NULL, 0.0, 0, ""},
    {"T341", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][3][0]), NULL, 0.0, 0, ""},
    {"T342", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][3][1]), NULL, 0.0, 0, ""},
    {"T343", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][3][2]), NULL, 0.0, 0, ""},
    {"T344", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][3][3]), NULL, 0.0, 0, ""},
    {"T351", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][4][0]), NULL, 0.0, 0, ""},
    {"T352", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][4][1]), NULL, 0.0, 0, ""},
    {"T353", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][4][2]), NULL, 0.0, 0, ""},
    {"T354", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][4][3]), NULL, 0.0, 0, ""},
    {"T355", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][4][4]), NULL, 0.0, 0, ""},
    {"T361", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][5][0]), NULL, 0.0, 0, ""},
    {"T362", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][5][1]), NULL, 0.0, 0, ""},
    {"T363", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][5][2]), NULL, 0.0, 0, ""},
    {"T364", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][5][3]), NULL, 0.0, 0, ""},
    {"T365", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][5][4]), NULL, 0.0, 0, ""},
    {"T366", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[2][5][5]), NULL, 0.0, 0, ""},
    {"T411", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][0][0]), NULL, 0.0, 0, ""},
    {"T421", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][1][0]), NULL, 0.0, 0, ""},
    {"T422", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][1][1]), NULL, 0.0, 0, ""},
    {"T431", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][2][0]), NULL, 0.0, 0, ""},
    {"T432", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][2][1]), NULL, 0.0, 0, ""},
    {"T433", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][2][2]), NULL, 0.0, 0, ""},
    {"T441", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][3][0]), NULL, 0.0, 0, ""},
    {"T442", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][3][1]), NULL, 0.0, 0, ""},
    {"T443", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][3][2]), NULL, 0.0, 0, ""},
    {"T444", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][3][3]), NULL, 0.0, 0, ""},
    {"T451", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][4][0]), NULL, 0.0, 0, ""},
    {"T452", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][4][1]), NULL, 0.0, 0, ""},
    {"T453", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][4][2]), NULL, 0.0, 0, ""},
    {"T454", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][4][3]), NULL, 0.0, 0, ""},
    {"T455", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][4][4]), NULL, 0.0, 0, ""},
    {"T461", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][5][0]), NULL, 0.0, 0, ""},
    {"T462", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][5][1]), NULL, 0.0, 0, ""},
    {"T463", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][5][2]), NULL, 0.0, 0, ""},
    {"T464", "1", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][5][3]), NULL, 0.0, 0, ""},
    {"T465", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][5][4]), NULL, 0.0, 0, ""},
    {"T466", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[3][5][5]), NULL, 0.0, 0, ""},
    {"T511", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][0][0]), NULL, 0.0, 0, ""},
    {"T521", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][1][0]), NULL, 0.0, 0, ""},
    {"T522", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][1][1]), NULL, 0.0, 0, ""},
    {"T531", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][2][0]), NULL, 0.0, 0, ""},
    {"T532", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][2][1]), NULL, 0.0, 0, ""},
    {"T533", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][2][2]), NULL, 0.0, 0, ""},
    {"T541", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][3][0]), NULL, 0.0, 0, ""},
    {"T542", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][3][1]), NULL, 0.0, 0, ""},
    {"T543", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][3][2]), NULL, 0.0, 0, ""},
    {"T544", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][3][3]), NULL, 0.0, 0, ""},
    {"T551", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][4][0]), NULL, 0.0, 0, ""},
    {"T552", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][4][1]), NULL, 0.0, 0, ""},
    {"T553", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][4][2]), NULL, 0.0, 0, ""},
    {"T554", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][4][3]), NULL, 0.0, 0, ""},
    {"T555", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][4][4]), NULL, 0.0, 0, ""},
    {"T561", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][5][0]), NULL, 0.0, 0, ""},
    {"T562", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][5][1]), NULL, 0.0, 0, ""},
    {"T563", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][5][2]), NULL, 0.0, 0, ""},
    {"T564", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][5][3]), NULL, 0.0, 0, ""},
    {"T565", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][5][4]), NULL, 0.0, 0, ""},
    {"T566", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[4][5][5]), NULL, 0.0, 0, ""},
    {"T611", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][0][0]), NULL, 0.0, 0, ""},
    {"T621", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][1][0]), NULL, 0.0, 0, ""},
    {"T622", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][1][1]), NULL, 0.0, 0, ""},
    {"T631", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][2][0]), NULL, 0.0, 0, ""},
    {"T632", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][2][1]), NULL, 0.0, 0, ""},
    {"T633", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][2][2]), NULL, 0.0, 0, ""},
    {"T641", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][3][0]), NULL, 0.0, 0, ""},
    {"T642", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][3][1]), NULL, 0.0, 0, ""},
    {"T643", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][3][2]), NULL, 0.0, 0, ""},
    {"T644", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][3][3]), NULL, 0.0, 0, ""},
    {"T651", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][4][0]), NULL, 0.0, 0, ""},
    {"T652", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][4][1]), NULL, 0.0, 0, ""},
    {"T653", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][4][2]), NULL, 0.0, 0, ""},
    {"T654", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][4][3]), NULL, 0.0, 0, ""},
    {"T655", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][4][4]), NULL, 0.0, 0, ""},
    {"T661", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][5][0]), NULL, 0.0, 0, ""},
    {"T662", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][5][1]), NULL, 0.0, 0, ""},
    {"T663", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][5][2]), NULL, 0.0, 0, ""},
    {"T664", "1", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][5][3]), NULL, 0.0, 0, ""},
    {"T665", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][5][4]), NULL, 0.0, 0, ""},
    {"T666", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ematrix_example.T[5][5][5]), NULL, 0.0, 0, ""},
    };

TFBPICKUP tfbPickup_example;

PARAMETER tfbpickup_param[N_TFBPICKUP_PARAMS] = {
   {"ID", "", IS_STRING, 0, (long)((char*)&tfbPickup_example.ID), NULL, 0.0, 0, "System identifier"},
   {"PLANE", "", IS_STRING, 0, (long)((char*)&tfbPickup_example.plane), "x", 0.0, 0, "\"x\", \"y\", \"delta\", or \"phase\""},
   {"RMS_NOISE", "M", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.rmsNoise), NULL, 0.0, 0, "RMS noise to add to position readings."},
   {"A0", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[0]), NULL, 0.0, 0, "Filter coefficient"},
   {"A1", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[1]), NULL, 0.0, 0, "Filter coefficient"},
   {"A2", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[2]), NULL, 0.0, 0, "Filter coefficient"},
   {"A3", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[3]), NULL, 0.0, 0, "Filter coefficient"},
   {"A4", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[4]), NULL, 0.0, 0, "Filter coefficient"},
   {"A5", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[5]), NULL, 0.0, 0, "Filter coefficient"},
   {"A6", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[6]), NULL, 0.0, 0, "Filter coefficient"},
   {"A7", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[7]), NULL, 0.0, 0, "Filter coefficient"},
   {"A8", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[8]), NULL, 0.0, 0, "Filter coefficient"},
   {"A9", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[9]), NULL, 0.0, 0, "Filter coefficient"},
   {"A10", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[10]), NULL, 0.0, 0, "Filter coefficient"},
   {"A11", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[11]), NULL, 0.0, 0, "Filter coefficient"},
   {"A12", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[12]), NULL, 0.0, 0, "Filter coefficient"},
   {"A13", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[13]), NULL, 0.0, 0, "Filter coefficient"},
   {"A14", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[14]), NULL, 0.0, 0, "Filter coefficient"},
   {"A15", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[15]), NULL, 0.0, 0, "Filter coefficient"},
   {"A16", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[16]), NULL, 0.0, 0, "Filter coefficient"},
   {"A17", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[17]), NULL, 0.0, 0, "Filter coefficient"},
   {"A18", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[18]), NULL, 0.0, 0, "Filter coefficient"},
   {"A19", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[19]), NULL, 0.0, 0, "Filter coefficient"},
   {"A20", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[20]), NULL, 0.0, 0, "Filter coefficient"},
   {"A21", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[21]), NULL, 0.0, 0, "Filter coefficient"},
   {"A22", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[22]), NULL, 0.0, 0, "Filter coefficient"},
   {"A23", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[23]), NULL, 0.0, 0, "Filter coefficient"},
   {"A24", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[24]), NULL, 0.0, 0, "Filter coefficient"},
   {"A25", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[25]), NULL, 0.0, 0, "Filter coefficient"},
   {"A26", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[26]), NULL, 0.0, 0, "Filter coefficient"},
   {"A27", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[27]), NULL, 0.0, 0, "Filter coefficient"},
   {"A28", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[28]), NULL, 0.0, 0, "Filter coefficient"},
   {"A29", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.a[29]), NULL, 0.0, 0, "Filter coefficient"},
   {"UPDATE_INTERVAL", "", IS_LONG, 0, (long)((char*)&tfbPickup_example.updateInterval), NULL, 0.0, 0, "Interval in turns for sampling data and updating filter output."},
   {"START_PASS", "", IS_LONG, 0, (long)((char*)&tfbPickup_example.startPass), NULL, 0.0, -1, "If positive, first pass on which to perform computations."},
   {"END_PASS", "", IS_LONG, 0, (long)((char*)&tfbPickup_example.endPass), NULL, 0.0, -1, "If positive, last pass on which to perform computations."},
   {"REFERENCE_FREQUENCY", "", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.referenceFrequency), NULL, 0.0, 0, "Reference frequency for computing phase offsets."},
   {"DX", "M", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.dx), NULL, 0.0, 0, "Horizontal offset (subtracted from pickup signal)."},
   {"DY", "M", IS_DOUBLE, 0, (long)((char*)&tfbPickup_example.dy), NULL, 0.0, 0, "Vertical offset (subtracted from pickup signal)"},
   {"BUNCHED_BEAM_MODE", "", IS_SHORT, 0, (long)((char*)&tfbPickup_example.bunchedBeamMode), NULL, 0.0, 1, "If non-zero, run in bunched beam mode."},
} ;

TFBDRIVER tfbDriver_example;

PARAMETER tfbdriver_param[N_TFBDRIVER_PARAMS] = {
   {"ID", "", IS_STRING, 0, (long)((char*)&tfbDriver_example.ID), NULL, 0.0, 0, "System identifier"},
   {"STRENGTH", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.strength), NULL, 0.0, 0, "Strength factor"},
   {"KICK_LIMIT", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.kickLimit), NULL, 0.0, 0, "Limit on applied kick, nominally in radians."},
   {"FREQUENCY", "Hz", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.frequency), NULL, 0.0, 0, "Resonant frequency of the unloaded kicker cavity."},
   {"DRIVE_FREQUENCY", "Hz", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.driveFrequency), NULL, 0.0, 0, "Drive frequency. If zero, defaults to resonant frequency of the loaded cavity."},
   {"CLOCK_FREQUENCY", "Hz", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.clockFrequency), NULL, 0.0, 0, "Clock frequency used for timing of the changes to generator current. Typically the rf or bunch frequency is used."},
   {"CLOCK_OFFSET", "s", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.clockOffset), NULL, 0.0, 0, "Offset of the generator current change relative to clock tick. Clock tick is nominally aligned to the bunch center."},
   {"PHASE", "Deg", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.phase), NULL, 0.0, 0, "Phase of the applied voltage relative to the bunch center, with 0 being on-crest.x2"},
   {"RAOVERQ", "Ohm", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.RaOverQ), NULL, 0.0, 0, "Shunt impedance, Ra/Q=V^2/(P*Q)."},
   {"QLOADED", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.QLoaded), NULL, 0.0, 0, "Loaded Q of the cavity."},
   {"OUTPUT_FILE", "", IS_STRING, 0, (long)((char*)&tfbDriver_example.outputFile), NULL, 0.0, 0, "File for logging filter output and driver output"},
   {"DELAY", "", IS_LONG, 0, (long)((char*)&tfbDriver_example.delay), NULL, 0.0, 0, "Delay (in turns)"},
   {"A0", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[0]), NULL, 1.0, 0, "Filter coefficient"},
   {"A1", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[1]), NULL, 0.0, 0, "Filter coefficient"},
   {"A2", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[2]), NULL, 0.0, 0, "Filter coefficient"},
   {"A3", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[3]), NULL, 0.0, 0, "Filter coefficient"},
   {"A4", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[4]), NULL, 0.0, 0, "Filter coefficient"},
   {"A5", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[5]), NULL, 0.0, 0, "Filter coefficient"},
   {"A6", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[6]), NULL, 0.0, 0, "Filter coefficient"},
   {"A7", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[7]), NULL, 0.0, 0, "Filter coefficient"},
   {"A8", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[8]), NULL, 0.0, 0, "Filter coefficient"},
   {"A9", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[9]), NULL, 0.0, 0, "Filter coefficient"},
   {"A10", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[10]), NULL, 0.0, 0, "Filter coefficient"},
   {"A11", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[11]), NULL, 0.0, 0, "Filter coefficient"},
   {"A12", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[12]), NULL, 0.0, 0, "Filter coefficient"},
   {"A13", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[13]), NULL, 0.0, 0, "Filter coefficient"},
   {"A14", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[14]), NULL, 0.0, 0, "Filter coefficient"},
   {"A15", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[15]), NULL, 0.0, 0, "Filter coefficient"},
   {"A16", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[16]), NULL, 0.0, 0, "Filter coefficient"},
   {"A17", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[17]), NULL, 0.0, 0, "Filter coefficient"},
   {"A18", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[18]), NULL, 0.0, 0, "Filter coefficient"},
   {"A19", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[19]), NULL, 0.0, 0, "Filter coefficient"},
   {"A20", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[20]), NULL, 0.0, 0, "Filter coefficient"},
   {"A21", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[21]), NULL, 0.0, 0, "Filter coefficient"},
   {"A22", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[22]), NULL, 0.0, 0, "Filter coefficient"},
   {"A23", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[23]), NULL, 0.0, 0, "Filter coefficient"},
   {"A24", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[24]), NULL, 0.0, 0, "Filter coefficient"},
   {"A25", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[25]), NULL, 0.0, 0, "Filter coefficient"},
   {"A26", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[26]), NULL, 0.0, 0, "Filter coefficient"},
   {"A27", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[27]), NULL, 0.0, 0, "Filter coefficient"},
   {"A28", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[28]), NULL, 0.0, 0, "Filter coefficient"},
   {"A29", "", IS_DOUBLE, 0, (long)((char*)&tfbDriver_example.a[29]), NULL, 0.0, 0, "Filter coefficient"},
   {"UPDATE_INTERVAL", "", IS_LONG, 0, (long)((char*)&tfbDriver_example.updateInterval), NULL, 0.0, 0, "Interval in units of pickup update interval for sampling pickup data and updating filter output."},
   {"OUTPUT_INTERVAL", "", IS_LONG, 0, (long)((char*)&tfbDriver_example.outputInterval), NULL, 0.0, 1024, "Number of samples to buffer between writing output file updates."},
   {"START_PASS", "", IS_LONG, 0, (long)((char*)&tfbDriver_example.startPass), NULL, 0.0, -1, "If positive, first pass on which to drive beam."},
   {"END_PASS", "", IS_LONG, 0, (long)((char*)&tfbDriver_example.endPass), NULL, 0.0, -1, "If positive, last pass on which to drive beam."},
   {"LONGITUDINAL", "", IS_SHORT, 0, (long)((char*)&tfbDriver_example.longitudinal), NULL, 0.0, 0, "If non-zero, kick is in the longituidinal plane. KICK_LIMIT is in fractional momentum deviation."},
   {"BUNCHED_BEAM_MODE", "", IS_SHORT, 0, (long)((char*)&tfbDriver_example.bunchedBeamMode), NULL, 0.0, 1, "If non-zero, run in bunched beam mode."},
} ;

LSCDRIFT lscdrift_example;

PARAMETER lscdrift_param[N_LSCDRIFT_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lscdrift_example.length), NULL, 0.0, 0, "length"},
    {"LEFFECTIVE", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lscdrift_example.lEffective), NULL, 0.0, 0, "effective length (used if L=0)"},
    {"BINS", "", IS_LONG, 0, (long)((char *)&lscdrift_example.bins), NULL, 0.0, 0, "number of bins for current histogram"},
    {"SMOOTHING", "", IS_SHORT, 0, (long)((char *)&lscdrift_example.smoothing), NULL, 0.0, 0, "Use Savitzky-Golay filter to smooth current histogram?"},
    {"SG_HALFWIDTH", "", IS_SHORT, 0, (long)((char *)&lscdrift_example.SGHalfWidth), NULL, 0.0, 1, "Savitzky-Golay filter half-width for smoothing current histogram"},
    {"SG_ORDER", "", IS_SHORT, 0, (long)((char *)&lscdrift_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing current histogram"},
    {"INTERPOLATE", "", IS_SHORT, 0, (long)((char *)&lscdrift_example.interpolate), NULL, 0.0, 1, "Interpolate wake?"},
    {"LSC", "", IS_SHORT, 0, (long)((char *)&lscdrift_example.lsc), NULL, 0.0, 1, "Include longitudinal space-charge impedance?  If zero, acts like ordinary drift."},
    {"AUTO_LEFFECTIVE", "", IS_SHORT, 0, (long)((char *)&lscdrift_example.autoLEffective), NULL, 0.0, 0, "In nonzero and if L=0, the LEFFECTIVE parameter is set to the length of the previous element."},
    {"LOW_FREQUENCY_CUTOFF0", "", IS_DOUBLE, 0, (long)((char*)&lscdrift_example.lowFrequencyCutoff0), NULL, -1.0, 0, "Highest spatial frequency at which low-frequency cutoff filter is zero.  If not positive, no low-frequency cutoff filter is applied. Frequency is in units of Nyquist (0.5/binsize)."},
    {"LOW_FREQUENCY_CUTOFF1", "", IS_DOUBLE, 0, (long)((char*)&lscdrift_example.lowFrequencyCutoff1), NULL, -1.0, 0, "Lowest spatial frequency at which low-frequency cutoff filter is 1.  If not given, defaults to LOW_FREQUENCY_CUTOFF1."},
    {"HIGH_FREQUENCY_CUTOFF0", "", IS_DOUBLE, 0, (long)((char*)&lscdrift_example.highFrequencyCutoff0), NULL, -1.0, 0, "Spatial frequency at which smoothing filter begins.  If not positive, no frequency filter smoothing is done.  Frequency is in units of Nyquist (0.5/binsize)."},
    {"HIGH_FREQUENCY_CUTOFF1", "", IS_DOUBLE, 0, (long)((char*)&lscdrift_example.highFrequencyCutoff1), NULL, -1.0, 0, "Spatial frequency at which smoothing filter is 0.  If not given, defaults to HIGH_FREQUENCY_CUTOFF0."},
    {"RADIUS_FACTOR", "", IS_DOUBLE, 0, (long)((char*)&lscdrift_example.radiusFactor), NULL, 1.7, 0, "LSC radius is (Sx+Sy)/2*RADIUS_FACTOR"},
};

LSRMDLTR lsrMdltr_example;
PARAMETER lsrMdltr_param[N_LSRMDLTR_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lsrMdltr_example.length), NULL, 0.0, 0, "length"},
    {"BU", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lsrMdltr_example.Bu), NULL, 0.0, 0, "Undulator peak field"},
    {"TGU_GRADIENT", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lsrMdltr_example.tguGradient), NULL, 0.0, 0, "Transverse gradient divided by maximum on-axis field."},
    {"TGU_COMP_FACTOR", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lsrMdltr_example.tguCompFactor), NULL, 1.0, 0, "Use to adjust constant field component to reduce trajectory error."},
    {"PERIODS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&lsrMdltr_example.periods), NULL, 0.0, 0, "Number of undulator periods."},
    {"METHOD", NULL, IS_STRING, 0, (long)((char*)&lsrMdltr_example.method), "non-adaptive runge-kutta", 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, modified-midpoint, two-pass modified-midpoint, leap-frog, non-adaptive runge-kutta)"},
    {"FIELD_EXPANSION", NULL, IS_STRING, 0, (long)((char*)&lsrMdltr_example.fieldExpansion), "leading terms", 0.0, 0, "ideal, exact, or \"leading terms\""},
    {"ACCURACY", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lsrMdltr_example.accuracy), NULL, 0.0, 0, "Integration accuracy for adaptive integration. (Not recommended)"},
    {"N_STEPS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&lsrMdltr_example.nSteps), NULL, 0.0, 0, "Number of integration steps for non-adaptive integration."},
    {"POLE_FACTOR1", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lsrMdltr_example.poleFactor1), NULL, 1.557150345503998e-01, 0, "Strength factor for the first and last pole."},
    {"POLE_FACTOR2", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lsrMdltr_example.poleFactor2), NULL, 3.806870122885810e-01, 0, "Strength factor for the second and second-to-last pole."},
    {"POLE_FACTOR3", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&lsrMdltr_example.poleFactor3), NULL, 8.028293373481792e-01, 0, "Strength factor for the third and third-to-last pole."},
    {"LASER_WAVELENGTH", "M", IS_DOUBLE, 0, (long)((char *)&lsrMdltr_example.usersLaserWavelength), NULL, 0.0, 0, "Laser wavelength. If zero, the wavelength is calculated from the resonance condition."},
    {"LASER_PEAK_POWER", "W", IS_DOUBLE, 0, (long)((char *)&lsrMdltr_example.laserPeakPower), NULL, 0.0, 0, "laser peak power"},
    {"LASER_W0", "M", IS_DOUBLE, 0, (long)((char *)&lsrMdltr_example.laserW0), NULL, 1.0, 0, "laser spot size at waist, $w_0 = \\sqrt{2}\\sigma_x = \\sqrt{2}\\sigma_y$"},
    {"LASER_PHASE", "RAD", IS_DOUBLE, 0, (long)((char *)&lsrMdltr_example.laserPhase), NULL, 0.0, 0, "laser phase"},
    {"LASER_X0", "M", IS_DOUBLE, 0, (long)((char *)&lsrMdltr_example.laserX0), NULL, 0.0, 0, "laser horizontal offset at center of wiggler"},
    {"LASER_Y0", "M", IS_DOUBLE, 0, (long)((char *)&lsrMdltr_example.laserY0), NULL, 0.0, 0, "laser vertical offset at center of wiggler"},
    {"LASER_Z0", "M", IS_DOUBLE, 0, (long)((char *)&lsrMdltr_example.laserZ0), NULL, 0.0, 0, "offset of waist position from center of wiggler"},
    {"LASER_TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&lsrMdltr_example.laserTilt), NULL, 0.0, 0, "laser tilt"},
    {"LASER_M", "", IS_SHORT, 0, (long)((char *)&lsrMdltr_example.laserM), NULL, 0.0, 0, "laser horizontal mode number (<5)"},
    {"LASER_N", "", IS_SHORT, 0, (long)((char *)&lsrMdltr_example.laserN), NULL, 0.0, 0, "laser vertical mode number (<5)"},
    {"SYNCH_RAD", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&lsrMdltr_example.synchRad), NULL, 0.0, 0, "Include classical, single-particle synchrotron radiation?"},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&lsrMdltr_example.isr), NULL, 0.0, 0, "Include quantum excitation?"},
    {"HELICAL", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&lsrMdltr_example.helical), NULL, 0.0, 0, "If non-zero, simulate helical undulator."},
    {"TIME_PROFILE", NULL, IS_STRING, PARAM_XY_WAVEFORM, (long)((char*)&lsrMdltr_example.timeProfileFile), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving time-dependent modulation of the laser electric and magnetic fields."},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&lsrMdltr_example.timeProfileOffset), NULL, 0.0, 0, "Time offset of the laser profile."},
};  

EDRIFT edrift_example;
PARAMETER edrift_param[N_EDRIFT_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&edrift_example.length), NULL, 0.0, 0, "length"},
};

SCMULT scmult_example;   
/*
PARAMETER scmult_param[N_SCMULT_PARAMS] = {
};
*/
/*#ifdef WIN32*/
PARAMETER scmult_param[1] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&edrift_example.length), NULL, 0.0, 0, "length"},
};
/*#else
PARAMETER scmult_param[] = {
};
#endif*/

ILMATRIX ilmatrix_example;
PARAMETER ilmatrix_param[N_ILMATRIX_PARAMS]={
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&ilmatrix_example.length), NULL, 0.0, 0, "Length (used for position and time-of-flight computation)"},
    {"NUX", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.tune[0]), NULL, 0.0, 0, "Horizontal tune"},
    {"NUY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.tune[1]), NULL, 0.0, 0, "Vertical tune"},
    {"NUX1M", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.chrom[0]), NULL, 0.0, 0, "First chromatic derivative of the horizontal tune"},
    {"NUY1M", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.chrom[1]), NULL, 0.0, 0, "First chromatic derivative of the vertical tune"},
    {"NUX2M", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.chrom2[0]), NULL, 0.0, 0, "Second chromatic derivative of the horizontal tune"},
    {"NUY2M", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.chrom2[1]), NULL, 0.0, 0, "Second chromatic derivative of the vertical tune"},
    {"NUX3M", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.chrom3[0]), NULL, 0.0, 0, "Third chromatic derivative of the horizontal tune"},
    {"NUY3M", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.chrom3[1]), NULL, 0.0, 0, "Third chromatic derivative of the vertical tune"},
    {"NUX1AX", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.tswax[0]), NULL, 0.0, 0, "First amplitude derivative of the horizontal tune wrt Ax"},
    {"NUY1AX", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.tswax[1]), NULL, 0.0, 0, "First amplitude derivative of the vertical tune wrt Ax"},
    {"NUX1AY", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.tsway[0]), NULL, 0.0, 0, "First amplitude derivative of the horizontal tune wrt Ay"},
    {"NUY1AY", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.tsway[1]), NULL, 0.0, 0, "First amplitude derivative of the vertical tune wrt Ay"},
    {"NUX2AX", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.tswax2[0]), NULL, 0.0, 0, "Second amplitude derivative of the horizontal tune wrt Ax"},
    {"NUY2AX", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.tswax2[1]), NULL, 0.0, 0, "Second amplitude derivative of the vertical tune wrt Ax"},
    {"NUX2AY", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.tsway2[0]), NULL, 0.0, 0, "Second amplitude derivative of the horizontal tune wrt Ay"},
    {"NUY2AY", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.tsway2[1]), NULL, 0.0, 0, "Second amplitude derivative of the vertical tune wrt Ay"},
    {"NUX1AX1AY", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.tswaxay[0]), NULL, 0.0, 0, "Amplitude derivative of the horizontal tune wrt Ax and Ay"},
    {"NUY1AX1AY", "1/M^2", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.tswaxay[1]), NULL, 0.0, 0, "Amplitude derivative of the vertical tune wrt Ax and Ay"},
    {"BETAX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.beta[0]), NULL, 0.0, 0, "On-momentum horizontal beta function"},
    {"BETAY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.beta[1]), NULL, 0.0, 0, "On-momentum vertical beta function"},
    {"BETAX1M", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.beta1[0]), NULL, 0.0, 0, "First chromatic derivative of horizontal beta function"},
    {"BETAY1M", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.beta1[1]), NULL, 0.0, 0, "First chromatic derivative of vertical beta function"},
    {"ALPHAX", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.alpha[0]), NULL, 0.0, 0, "On-momentum horizontal alpha function"},
    {"ALPHAY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.alpha[1]), NULL, 0.0, 0, "On-momentum vertical alpha function"},
    {"ALPHAX1M", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.alpha1[0]), NULL, 0.0, 0, "First chromatic derivative of horizontal alpha function"},
    {"ALPHAY1M", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.alpha1[1]), NULL, 0.0, 0, "First chromatic derivative of vertical alpha function"},
    {"ETAX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.eta[0]), NULL, 0.0, 0, "On-momentum horizontal eta function"},
    {"ETAPX", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.eta[1]), NULL, 0.0, 0, "On-momentum horizontal eta' function"},
    {"ETAY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.eta[2]), NULL, 0.0, 0, "On-momentum vertical eta function"},
    {"ETAPY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.eta[3]), NULL, 0.0, 0, "On-momentum vertical eta' function"},
    {"ETAX1", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.eta1[0]), NULL, 0.0, 0, "First chromatic derivative of horizontal eta function"},
    {"ETAPX1", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.eta1[1]), NULL, 0.0, 0, "First chromatic derivative of horizontal eta' function"},
    {"ETAY1", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.eta1[2]), NULL, 0.0, 0, "First chromatic derivative of vertical eta function"},
    {"ETAPY1", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.eta1[3]), NULL, 0.0, 0, "First chromatic derivative of vertical eta' function"},
    {"ALPHAC", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.alphac[0]), NULL, 0.0, 0, "First-order momentum compaction factor"},
    {"ALPHAC2", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.alphac[1]), NULL, 0.0, 0, "Second-order momentum compaction factor"},
    {"ALPHAC3", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.alphac[2]), NULL, 0.0, 0, "Third-order momentum compaction factor"},
    {"DS1AX", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.dsdA[0]), NULL, 0.0, 0, "First amplitude derivative of the path length wrt Ax"},
    {"DS1AY", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.dsdA[1]), NULL, 0.0, 0, "First amplitude derivative of the path length wrt Ay"},
    {"DS2AX", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.dsdA2[0]), NULL, 0.0, 0, "Second amplitude derivative of the path length wrt Ax"},
    {"DS2AY", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.dsdA2[1]), NULL, 0.0, 0, "Second amplitude derivative of the path length wrt Ay"},
    {"DS1AX1AY", "1/M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.dsdAxAy), NULL, 0.0, 0, "Amplitude derivative of the path length wrt Ax and Ay"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ilmatrix_example.tilt), NULL, 0.0, 0, "Rotation angle about the longitudinal axis."},
    {"CROSS_RESONANCE", "", IS_SHORT, 0, (long)((char *)&ilmatrix_example.allowResonanceCrossing), NULL, 0.0, 0, "If zero, then particles that cross an integer or half-integer resonance are considered lost."},
    {"VERBOSITY", "", IS_SHORT, 0, (long)((char *)&ilmatrix_example.verbosity), NULL, 0.0, 0, "If nonzero, then information about particle losses is printed out."}
};

TSCATTER tscatter_example;   
PARAMETER tscatter_param[N_TSCATTER_PARAMS] = {
    {"DUMMY", "", IS_LONG, 0, 0, NULL, 0.0, 0, ""},
};

KQUSE kquse_example;

/* kick quad+sext physical parameters */
PARAMETER kquse_param[N_KQUSE_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&kquse_example.length), NULL, 0.0, 0, "length"},
    {"K1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquse_example.k1), NULL, 0.0, 0, "geometric quadrupole strength"},
    {"K2", "1/M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquse_example.k2), NULL, 0.0, 0, "geometric sextupole strength"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquse_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquse_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquse_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquse_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE1", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquse_example.fse1), NULL, 0.0, 0, "fractional strength error for K1"},
    {"FSE2", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&kquse_example.fse2), NULL, 0.0, 0, "fractional strength error for K2"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&kquse_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks"},
    {"INTEGRATION_ORDER", "", IS_SHORT, 0, (long)((char *)&kquse_example.integration_order), NULL, 0.0, 4, "integration order (2 or 4)"},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&kquse_example.synch_rad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&kquse_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"ISR1PART", "", IS_SHORT, 0, (long)((char *)&kquse_example.isr1Particle), NULL, 0.0, 1, "Include ISR for single-particle beam only if ISR=1 and ISR1PART=1"},
    {"MATRIX_TRACKING", "", IS_SHORT, 0, (long)((char *)&kquse_example.matrixTracking), NULL, 0.0, 0, "For testing only."},
    {"EXPAND_HAMILTONIAN", "", IS_SHORT, 0, (long)((char *)&kquse_example.expandHamiltonian), NULL, 0.0, 0, "If 1, Hamiltonian is expanded to leading order."},
    };


UKICKMAP ukickmap_example;

/* undulator kick map physical parameters */
PARAMETER ukickmap_param[N_UKICKMAP_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ukickmap_example.length), NULL, 0.0, 0, "length"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ukickmap_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ukickmap_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ukickmap_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ukickmap_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FIELD_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ukickmap_example.fieldFactor), NULL, 1.0, 0, "Factor by which to multiply the magnetic fields."},
    {"XY_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ukickmap_example.xyFactor), NULL, 1.0, 0, "Factor by which to multiply the x and y values in the input file."},
    {"YAW", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ukickmap_example.yaw), NULL, 0.0, 0, "Yaw angle of the device. Meaningful only if N_KICKS is not 1."},
    {"INPUT_FILE", " ", IS_STRING, 0, (long)((char *)&ukickmap_example.inputFile), NULL, 0.0, 0, "Name of SDDS file with undulator kickmap data."},
    {"N_KICKS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&ukickmap_example.nKicks), NULL, 0.0, 1, "Number of kicks into which to split the element."},
    {"PERIODS", "", IS_LONG, 0, (long)((char *)&ukickmap_example.periods), NULL, 0.0, 0, "Number of periods (for radiation integral computations only)."},
    {"KREF", "", IS_DOUBLE, 0, (long)((char *)&ukickmap_example.Kreference), NULL, 0.0, 0, "Reference value of undulator parameter. K=KREF*FIELD_FACTOR is used for radiation integral calculations only assuming period=L/PERIODS."},
    {"KACTUAL", "", IS_DOUBLE, 0, (long)((char *)&ukickmap_example.Kactual), NULL, 0.0, 0, "Value of undulator parameter, used for radiation integral calculations only assuming period=L/PERIODS."},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&ukickmap_example.synchRad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&ukickmap_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"YAW_END", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&ukickmap_example.yawEnd), NULL, 0.0, 0, "-1=Entrance, 0=Center, 1=Exit"},
    {"SINGLE_PERIOD_MAP", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&ukickmap_example.singlePeriodMap), NULL, 0.0, 0, "if non-zero, the map file is for a single period. L still pertains to the full device. Set N_KICKS to the number of periods."},
    };

FTABLE ftable_example;

/* field table physical parameters */
PARAMETER ftable_param[N_FTABLE_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.l0), NULL, 0.0, 0, "The effective field length measured along a straight line."},
    {"ANGLE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.angle), NULL, 0.0, 0, "The designed bending angle"},
    {"L1", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.l1), NULL, 0.0, 0, "The left fringe field length."},
    {"L2", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.l2), NULL, 0.0, 0, "The right fringe field length. L1+L+L2=Total z span in the input field table."},
    {"E1", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.e1), NULL, 0.0, 0, "The designed entrance edge angle"},
    {"E2", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.e2), NULL, 0.0, 0, "The designed exit edge angle"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.factor), NULL, 1.0, 0, "Factor by which to multiply field data."},
    {"THRESHOLD", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.threshold), NULL, 1e-8, 0, "Fields smaller than this are considered 0."},
    {"INPUT_FILE", "", IS_STRING, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.inputFile), NULL, 0.0, 0, "Name of SDDS file which contains field data."},
    {"N_KICKS", "", IS_LONG, PARAM_CHANGES_MATRIX, (long)((char *)&ftable_example.nKicks), NULL, 0.0, 1, "Number of kicks into which to split the element."},
    {"VERBOSE", "", IS_SHORT, 0, (long)((char *)&ftable_example.verbose), NULL, 0.0, 0, "Used for debugging code. Not applicable to Pelegant"},
    {"SIMPLE_INPUT", "", IS_SHORT, 0, (long)((char *)&ftable_example.simpleInput), NULL, 0.0, 0, "If non-zero, use simple input format."},
    };

/* emittance scaling element physical parameters */

EMITTANCEELEMENT emittanceElem_example;

PARAMETER emittanceElement_param[N_EMITTANCEELEMENT_PARAMS] = {
  {"EMITX", "M", IS_DOUBLE, 0, (long)((char *)&emittanceElem_example.emit[0]), NULL, -1.0, 0, "horizontal emittance"},
  {"EMITY", "M", IS_DOUBLE, 0, (long)((char *)&emittanceElem_example.emit[1]), NULL, -1.0, 0, "vertical emittance"},
  {"EMITNX", "M", IS_DOUBLE, 0, (long)((char *)&emittanceElem_example.emitn[0]), NULL, -1.0, 0, "horizontal normalized emittance"},
  {"EMITNY", "M", IS_DOUBLE, 0, (long)((char *)&emittanceElem_example.emitn[1]), NULL, -1.0, 0, "vertical normalized emittance"},
};

/* radiation integral scaling element  */
MRADINTEGRALS mRadIntegrals_example;
PARAMETER mRadIntegrals_param[N_MRADITEGRALS_PARAMS] = {
  {"FACTOR", "", IS_DOUBLE, 0, (long)((char *)&mRadIntegrals_example.factor), NULL, 1.0, 0, "factor"},
} ;  

MRFDF mrfdf_example;
/* names for multipole rf deflector parameters */
PARAMETER mrfdf_param[N_MRFDF_PARAMS] = {
    {"FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.factor), NULL, 1.0, 0, "A factor by which to multiply all components."},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"A1", "V/m", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.a[0]), NULL, 0.0, 0, "Vertically-deflecting dipole"},
    {"A2", "V/m$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.a[1]), NULL, 0.0, 0, "Skew quadrupole"},
    {"A3", "V/m$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.a[2]), NULL, 0.0, 0, "Skew sextupole"},
    {"A4", "V/m$a4$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.a[3]), NULL, 0.0, 0, "Skew octupole"},
    {"A5", "V/m$a5$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.a[4]), NULL, 0.0, 0, "Skew decapole"},
    {"B1", "V/m", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.b[0]), NULL, 0.0, 0, "Horizontally-deflecting dipole"},
    {"B2", "V/m$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.b[1]), NULL, 0.0, 0, "Normal quadrupole"},
    {"B3", "V/m$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.b[2]), NULL, 0.0, 0, "Normal sextupole"},
    {"B4", "V/m$a4$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.b[3]), NULL, 0.0, 0, "Normal octupole"},
    {"B5", "V/m$a5$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.b[4]), NULL, 0.0, 0, "Normal decapole"},
    {"FREQUENCY1", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.frequency[0]), NULL, DEFAULT_FREQUENCY, 0, "Dipole frequency"},
    {"FREQUENCY2", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.frequency[1]), NULL, DEFAULT_FREQUENCY, 0, "Quadrupole frequency"},
    {"FREQUENCY3", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.frequency[2]), NULL, DEFAULT_FREQUENCY, 0, "Sextupole frequency"},
    {"FREQUENCY4", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.frequency[3]), NULL, DEFAULT_FREQUENCY, 0, "Octupole frequency"},
    {"FREQUENCY5", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.frequency[4]), NULL, DEFAULT_FREQUENCY, 0, "Decapole frequency"},
    {"PHASE1", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.phase[0]), NULL, 0.0, 0, "Dipole phase"},
    {"PHASE2", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.phase[1]), NULL, 0.0, 0, "Quadrupole phase"},
    {"PHASE3", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.phase[2]), NULL, 0.0, 0, "Sextupole phase"},
    {"PHASE4", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.phase[3]), NULL, 0.0, 0, "Octupole phase"},
    {"PHASE5", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&mrfdf_example.phase[4]), NULL, 0.0, 0, "Decapole phase"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&mrfdf_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    } ;

EHCOR ehcor_example;
/* horizontal corrector physical parameters */
PARAMETER ehcor_param[N_EHCOR_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehcor_example.length), NULL, 0.0, 0, "length"},
    {"KICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehcor_example.kick), NULL, 0.0, 0, "kick angle"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehcor_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehcor_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehcor_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehcor_example.dz), NULL, 0.0, 0, "misalignment"},
    {"CALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehcor_example.calibration), NULL, 1.0, 0, "factor applied to obtain kick"},
    {"LERAD", "", IS_DOUBLE, 0, (long)((char *)&ehcor_example.lEffRad), NULL, 0.0, 0, "if L=0, use this length for radiation computations"},
    {"STEERING", "", IS_SHORT, 0, (long)((char *)&ehcor_example.steering), NULL, 0.0, 1, "use for steering?"},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&ehcor_example.synchRad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&ehcor_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"STEERING_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ehcor_example.steeringMultipoles), NULL, 0.0, 0, "input file for systematic multipole content of steering kicks"},
    {"RANDOM_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ehcor_example.randomMultipoles), NULL, 0.0, 0, "input file for random multipoles content of steering kicks"},
    {"RANDOM_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehcor_example.randomMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply random multipoles"},
    {"STEERING_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehcor_example.steeringMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply steering multipoles"},
    };

EVCOR evcor_example;
/* vertical corrector physical parameters */
PARAMETER evcor_param[N_EVCOR_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&evcor_example.length), NULL, 0.0, 0, "length"},
    {"KICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&evcor_example.kick), NULL, 0.0, 0, "kick angle"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&evcor_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&evcor_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&evcor_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&evcor_example.dz), NULL, 0.0, 0, "misalignment"},
    {"CALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&evcor_example.calibration), NULL, 1.0, 0, "factor applied to obtain kick"},
    {"LERAD", "", IS_DOUBLE, 0, (long)((char *)&evcor_example.lEffRad), NULL, 0.0, 0, "if L=0, use this length for radiation computations"},
    {"STEERING", "", IS_SHORT, 0, (long)((char *)&evcor_example.steering), NULL, 0.0, 1, "use for steering?"},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&evcor_example.synchRad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&evcor_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"STEERING_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&evcor_example.steeringMultipoles), NULL, 0.0, 0, "input file for systematic multipole content of steering kicks"},
    {"RANDOM_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&evcor_example.randomMultipoles), NULL, 0.0, 0, "input file for random multipoles content of steering kicks"},
    {"RANDOM_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&evcor_example.randomMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply random multipoles"},
    {"STEERING_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&evcor_example.steeringMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply steering multipoles"},
    };

EHVCOR ehvcor_example;
/* horizontal/vertical corrector physical parameters */
PARAMETER ehvcor_param[N_EHVCOR_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehvcor_example.length), NULL, 0.0, 0, "length"},
    {"HKICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehvcor_example.xkick), NULL, 0.0, 0, "horizontal kick angle"},
    {"VKICK", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehvcor_example.ykick), NULL, 0.0, 0, "vertical kick angle"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehvcor_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehvcor_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehvcor_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehvcor_example.dz), NULL, 0.0, 0, "misalignment"},
    {"HCALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehvcor_example.xcalibration), NULL, 1.0, 0, "factor applied to obtain horizontal kick"},
    {"VCALIBRATION", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehvcor_example.ycalibration), NULL, 1.0, 0, "factor applied to obtain vertical kick"},
    {"LERAD", "", IS_DOUBLE, 0, (long)((char *)&ehvcor_example.lEffRad), NULL, 0.0, 0, "if L=0, use this length for radiation computations"},
    {"STEERING", "", IS_SHORT, 0, (long)((char *)&ehvcor_example.steering), NULL, 0.0, 1, "use for steering?"},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&ehvcor_example.synchRad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&ehvcor_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"STEERING_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ehvcor_example.steeringMultipoles), NULL, 0.0, 0, "input file for systematic multipole content of steering kicks"},
    {"RANDOM_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ehvcor_example.randomMultipoles), NULL, 0.0, 0, "input file for random multipoles content of steering kicks"},
    {"RANDOM_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehvcor_example.randomMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply random multipoles"},
    {"STEERING_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ehvcor_example.steeringMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply steering multipoles"},
    };

BMAPXYZ bmapxyz_example;
PARAMETER bmapxyz_param[N_BMAPXYZ_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bmapxyz_example.length), NULL, 0.0, 0, "insertion length"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bmapxyz_example.dxError), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bmapxyz_example.dyError), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bmapxyz_example.dzError), NULL, 0.0, 0, "misalignment"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bmapxyz_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"LFIELD", "M", IS_DOUBLE, 0, (long)((char *)&bmapxyz_example.fieldLength), NULL, -1.0, 0, "expected length of the field map. If negative, determined from field data."},
    {"STRENGTH", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bmapxyz_example.strength), NULL, 1.0, 0, "factor by which to multiply field"},
    {"ACCURACY", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bmapxyz_example.accuracy), NULL, 0.0, 0, "integration accuracy"},
    {"METHOD", NULL, IS_STRING, PARAM_CHANGES_MATRIX, (long)((char*)&bmapxyz_example.method), NULL, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, modified-midpoint, two-pass modified-midpoint, leap-frog, non-adaptive runge-kutta"},
    {"FILENAME", NULL, IS_STRING, PARAM_CHANGES_MATRIX, (long)((char*)&bmapxyz_example.filename), NULL, 0.0, 0, "name of file containing columns (x, y, z) and either (Bx, By, Bz) or (Fx, Fy, Fz)"},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&bmapxyz_example.synchRad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"CHECK_FIELDS", "", IS_SHORT, 0, (long)((char *)&bmapxyz_example.checkFields), NULL, 0.0, 0, "check fields by computing divB and curlB errors?"},
    {"INJECT_AT_Z0", "", IS_SHORT, 0, (long)((char *)&bmapxyz_example.injectAtZero), NULL, 0.0, 0, "By default, particles are placed at the entrance to the field map regardless of the z coordinate values. If nonzero, particles start at z=0."},
    {"DRIFT_MATRIX", "", IS_SHORT, 0, (long)((char *)&bmapxyz_example.driftMatrix), NULL, 0.0, 0, "If non-zero, instead of tracking to determine the matrix, just assume a drift-space matrix."},
    {"PARTICLE_OUTPUT_FILE", NULL, IS_STRING, 0, (long)((char*)&bmapxyz_example.particleOutputFile), NULL, 0.0, 0, "name of file for phase-space output inside element. Use for debugging only in serial version."},
};  

BRAT brat_example;
PARAMETER brat_param[N_BRAT_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.length), NULL, 0.0, 0, "length"},
    {"ANGLE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.angle), NULL, 0.0, 0, "Nominal bending angle. Will be refined to match geometry specified by input/output and vertex coordinates"},
    {"FSE", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"ACCURACY", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.accuracy), NULL, 0.0, 0, "integration accuracy"},
    {"METHOD", NULL, IS_STRING, PARAM_CHANGES_MATRIX, (long)((char*)&brat_example.method), NULL, 0.0, 0, "Ignored. Method defaults to Bulirsch-Stoer."},
    {"FILENAME", NULL, IS_STRING, PARAM_CHANGES_MATRIX, (long)((char*)&brat_example.filename), NULL, 0.0, 0, "name of file containing columns (x, y, z, Bx, By, Bz)"},
    {"XVERTEX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.xVertex), NULL, 0.0, 0, "x coordinate of vertex point"},
    {"ZVERTEX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.zVertex), NULL, 0.0, 0, "z coordinate of vertex point"},
    {"XENTRY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.xEntry), NULL, 0.0, 0, "x coordinate of nominal entry point"},
    {"ZENTRY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.zEntry), NULL, 0.0, 0, "z coordinate of nominal entry point"},
    {"XEXIT", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.xExit), NULL, 0.0, 0, "x coordinate of nominal exit point"},
    {"ZEXIT", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.zExit), NULL, 0.0, 0, "z coordinate of nominal exit point"},
    {"DXMAP", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.dxMap), NULL, 0.0, 0, "x displacement of map"},
    {"DZMAP", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.dzMap), NULL, 0.0, 0, "z displacement of map"},
    {"YAWMAP", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.yawMap), NULL, 0.0, 0, "yaw of map about x=z=0"},
    {"FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.fieldFactor), NULL, 1.0, 0, "factor by which to multiply fields"},
    {"USE_FTABLE", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&brat_example.useFTABLE), NULL, 0.0, 0, "If nonzero, use FTABLE method for integration. Value gives the number of kicks."},
};  


BGGEXP bggexp_example;
PARAMETER bggexp_param[N_BGGEXP_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.length), NULL, 0.0, 0, "insertion length"},
    {"LFIELD", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.fieldLength), NULL, -1.0, 0, "expected length of the field map. If negative, use L."},
    {"FILENAME", NULL, IS_STRING, PARAM_CHANGES_MATRIX, (long)((char*)&bggexp_example.filename), NULL, 0.0, 0, "name of file generalized gradient data"},
    {"STRENGTH", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.strength), NULL, 1.0, 0, "factor by which to multiply field"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.dz), NULL, 0.0, 0, "misalignment"},
    {"BX", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.Bx), NULL, 0.0, 0, "add BX*STRENGTH to Bx field"},
    {"BY", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.By), NULL, 0.0, 0, "add BY*STRENGTH to By field"},
    {"MAXIMUM_M", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.mMaximum), NULL, 0.0, -1, "data with m greater than this is ignored"},
    {"MAXIMUM_2N", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.maximum2n), NULL, 0.0, -1, "data with 2*n greater than this is ignored"},
    {"Z_INTERVAL", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.zInterval), NULL, 0.0, 1, "input z data is sampled at this interval"},
    {"SYMPLECTIC", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.symplectic), NULL, 0.0, 0, "if nonzero, use implicit symplectic integrator. At minimum, should always be used to validate the sufficiency of the non-symplectic integrator."},
    {"SYNCH_RAD", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.synchRad), NULL, 0.0, 0, "if nonzero, include classical, single-particle synchrotron radiation"},
    {"ISR", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.isr), NULL, 0.0, 0, "if nonzero, include incoherent synchrotron radiation (quantum excitation)"},
    {"PARTICLE_OUTPUT_FILE", "", IS_STRING, 0, (long)((char*)&bggexp_example.particleOutputFile), NULL, 0.0, 0, "name of file for phase-space and field output. Use for debugging only!"},
    {"IS_BEND", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.isBend), NULL, 0.0, 0, "if nonzero, magnet is a bending magnet; vertex, entry, and exit points should be defined."},
    {"XVERTEX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.xVertex), NULL, 0.0, 0, "For dipoles: x position of vertex in coordinate system of the fields."},
    {"ZVERTEX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.zVertex), NULL, 0.0, 0, "For dipoles: z position of vertex in coordinate system of the fields."},
    {"XENTRY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.xEntry), NULL, 0.0, 0, "For dipoles: x position of reference entry point in coordinate system of the fields."},
    {"ZENTRY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.zEntry), NULL, 0.0, 0, "For dipoles: z position of reference entry point in coordinate system of the fields."},
    {"XEXIT", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.xExit), NULL, 0.0, 0, "For dipoles: x position of reference exit point in coordinate system of the fields."},
    {"ZEXIT", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.zExit), NULL, 0.0, 0, "For dipoles: z position of reference exit point in coordinate system of the fields."},
    {"DXEXPANSION", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&bggexp_example.dxExpansion), NULL, 0.0, 0, "x position of expansion relative to coordinate system of the fields."},
};  

IONEFFECTS ionEffects_example;
PARAMETER ionEffects_param[N_IONEFFECTS_PARAMS] = {
    {"DISABLE", "", IS_LONG, 0, (long)((char *)&ionEffects_example.disable), NULL, 0.0, 0, "If non-zero, turn off ion effects in the region covered by this element."},
    {"MACRO_IONS", "", IS_LONG, 0, (long)((char *)&ionEffects_example.macroIons), NULL, 0.0, 0, "If positive, overrides the default value set in the ion_effects command, giving the number of macro ions generated per bunch passage."},
    {"GENERATION_INTERVAL", "", IS_LONG, 0, (long)((char *)&ionEffects_example.generationInterval), NULL, 0.0, 0, "If positive, overrides the default value set in the ion_effects command, giving the number of macro ions generated per bunch passage."},
    {"X_SPAN", "", IS_DOUBLE, 0, (long)((char *)&ionEffects_example.span[0]), NULL, 0.0, 0, "If positive, gives the region over which ions are kept."},
    {"Y_SPAN", "", IS_DOUBLE, 0, (long)((char *)&ionEffects_example.span[1]), NULL, 0.0, 0, "If positive, gives the region over which ions are kept."},
    {"X_BIN_DIVISOR", "", IS_DOUBLE, 0, (long)((char *)&ionEffects_example.binDivisor[0]), NULL, 0.0, 0, "If positive, gives the ratio of electron beam sigma to bin size for ion field calculation."},
    {"Y_BIN_DIVISOR", "", IS_DOUBLE, 0, (long)((char *)&ionEffects_example.binDivisor[1]), NULL, 0.0, 0, "If positive, gives the ratio of electron beam sigma to bin size for ion field calculation."},
    {"X_RANGE_MULTIPLIER", "", IS_DOUBLE, 0, (long)((char *)&ionEffects_example.rangeMultiplier[0]), NULL, 0.0, 0, "If positive, gives the ratio of ion binning region size to ion 80% x range."},
    {"Y_RANGE_MULTIPLIER", "", IS_DOUBLE, 0, (long)((char *)&ionEffects_example.rangeMultiplier[1]), NULL, 0.0, 0, "If positive, gives the ratio of ion binning region size to ion 80% y range."},
    {"X_SIGMA_LIMIT_MULTIPLIER", "", IS_DOUBLE, 0, (long)((char *)&ionEffects_example.sigmaLimitMultiplier[0]), NULL, 0.0, 0, "If positive, gives lower limit on bi-gaussian fit sigma values in units of the ion bin size."},
    {"Y_SIGMA_LIMIT_MULTIPLIER", "", IS_DOUBLE, 0, (long)((char *)&ionEffects_example.sigmaLimitMultiplier[1]), NULL, 0.0, 0, "If positive, gives lower limit on bi-gaussian fit sigma values in units of the ion bin size."},
    {"STARTPASS", "", IS_LONG, 0, (long)((char *)&ionEffects_example.startPass), NULL, 0.0, 0, "If positive, gives the pass on which ion effects start."},
    {"ENDPASS", "", IS_LONG, 0, (long)((char *)&ionEffects_example.endPass), NULL, 0.0, -1, "If positive, gives the pass on which ion effects end."},
    {"PASSINTERVAL", "", IS_LONG, 0, (long)((char *)&ionEffects_example.passInterval), NULL, 0.0, 1, "Interval between ion effects modeling."},

};  

SLICE_POINT slice_point_example;
/* names for slice point */
PARAMETER slice_point_param[N_SLICE_POINT_PARAMS] = {
    {"N_SLICES", "", IS_LONG, 0, (long)((char *)&slice_point_example.nSlices), NULL, 0.0, 10, "number of slices"},
    {"START_PID", "", IS_LONG, 0, (long)((char *)&slice_point_example.startPID), NULL, 0.0, -1, "starting particleID for particles to dump"},
    {"END_PID", "", IS_LONG, 0, (long)((char *)&slice_point_example.endPID), NULL, 0.0, -1, "ending particleID for particles to dump"},
    {"INTERVAL", "", IS_LONG, 0, (long)((char *)&slice_point_example.interval), NULL, 0.0, 1, "interval for data output (in turns)"},
    {"START_PASS", "", IS_LONG, 0, (long)((char*)&slice_point_example.start_pass), NULL, 0.0, 0, "pass on which to start"},
    {"END_PASS", "", IS_LONG, 0, (long)((char*)&slice_point_example.end_pass), NULL, 0.0, -1, "pass on which to end (inclusive).  Ignored if negative."},
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&slice_point_example.filename), "", 0.0, 0, "output filename, possibly incomplete (see below)"},
    {"LABEL", "", IS_STRING, 0, (long)((char *)&slice_point_example.label), "", 0.0, 0, "output label"},
    {"INDEX_OFFSET", "", IS_LONG, 0, (long)((char *)&slice_point_example.indexOffset), NULL, 0.0, 0, "Offset for file indices for sequential file naming."},
    {"REFERENCE_FREQUENCY", "", IS_DOUBLE, 0, (long)((char *)&slice_point_example.referenceFrequency), NULL, -1.0, -1, "If non-zero, the indicated frequency is used to define the bucket center for purposes of computing time offsets."},
    {"DISABLE", "", IS_SHORT, 0, (long)((char *)&slice_point_example.disable), NULL, 0.0, 0, "If nonzero, no output will be generated."},    
    {"USE_DISCONNECT", "", IS_SHORT, 0, (long)((char *)&slice_point_example.useDisconnect), NULL, 0.0, 0, "If nonzero, files are disconnected between each write operation. May be useful for parallel operation.  Ignored otherwise."},
    } ;

SPEEDBUMP speedbump_example;
/* speed bump physical parameters */
PARAMETER speedbump_param[N_SPEEDBUMP_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&speedbump_example.length), NULL, 0.0, 0, "insertion length"},
    {"CHORD", "M", IS_DOUBLE, 0, (long)((char *)&speedbump_example.chord), NULL, 0.0, 0, "z length of speed bump"},
    {"DZCENTER", "M", IS_DOUBLE, 0, (long)((char *)&speedbump_example.dzCenter), NULL, 0.0, 0, "z center displacement of speed bump relative to middle of object"},
    {"HEIGHT", "M", IS_DOUBLE, 0, (long)((char *)&speedbump_example.height), NULL, 0.0, 0, "height above the surrounding chamber"},
    {"POSITION", "M", IS_DOUBLE, 0, (long)((char *)&speedbump_example.position), NULL, 0.0, 0, "position of peak relative to ideal trajectory"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&speedbump_example.dx), NULL, 0.0, 0, "horizontal misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&speedbump_example.dy), NULL, 0.0, 0, "vertical misalignment"},
    {"INSERT_FROM", "", IS_STRING, 0, (long)((char *)&speedbump_example.insertFrom), NULL, 0.0, 0, "direction from which inserted (x, +x, -x, y, +y, -y)"},
    };


CCBEND ccbend_example;
/* canonically-integrated rectangular bending magnet physical parameters */
PARAMETER ccbend_param[N_CCBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&ccbend_example.length), NULL, 0.0, 0, "arc length (not chord length!)"},
    {"ANGLE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX|PARAM_DIVISION_RELATED, (long)((char *)&ccbend_example.angle), NULL, 0.0, 0, "bend angle"},
    {"K1", "1/M$a2$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.K1), NULL, 0.0, 0, "geometric quadrupole strength"},
    {"K2", "1/M$a3$n", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.K2), NULL, 0.0, 0, "geometric sextupole strength"},
    {"K3", "1/M$a4$n", IS_DOUBLE, 0, (long)((char *)&ccbend_example.K3), NULL, 0.0, 0, "geometric octupole strength"},
    {"K4", "1/M$a5$n", IS_DOUBLE, 0, (long)((char *)&ccbend_example.K4), NULL, 0.0, 0, "geometric decapole strength"},
    {"K5", "1/M$a6$n", IS_DOUBLE, 0, (long)((char *)&ccbend_example.K5), NULL, 0.0, 0, "geometric 12-pole strength"},
    {"K6", "1/M$a7$n", IS_DOUBLE, 0, (long)((char *)&ccbend_example.K6), NULL, 0.0, 0, "geometric 14-pole strength"},
    {"K7", "1/M$a8$n", IS_DOUBLE, 0, (long)((char *)&ccbend_example.K7), NULL, 0.0, 0, "geometric 16-pole strength"},
    {"K8", "1/M$a9$n", IS_DOUBLE, 0, (long)((char *)&ccbend_example.K8), NULL, 0.0, 0, "geometric 18-pole strength"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"YAW", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.yaw), NULL, 0.0, 0, "rotation about vertical axis through entrance point"},
    {"HGAP", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.hgap), NULL, 0.0, 0, "half-gap between poles"},
    {"FINT1", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.fint1), NULL, 0.0, 0, "edge integral for entrance"},
    {"FINT2", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.fint2), NULL, 0.0, 0, "edge integral for exit"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"FSE_DIPOLE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.fseDipole), NULL, 0.0, 0, "fractional strength error of dipole component"},
    {"FSE_QUADRUPOLE", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.fseQuadrupole), NULL, 0.0, 0, "fractional strength error of quadrupole component"},
    {"ETILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.etilt), NULL, 0.0, 0, "error rotation about incoming longitudinal axis"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&ccbend_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks"},
    {"INTEGRATION_ORDER", "", IS_SHORT, 0, (long)((char *)&ccbend_example.integration_order), NULL, 0.0, 4, "integration order (2 or 4)"},
    {"SYSTEMATIC_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ccbend_example.systematic_multipoles), NULL, 0.0, 0, "input file for systematic multipoles"},
    {"EDGE_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ccbend_example.edge_multipoles), NULL, 0.0, 0, "input file for systematic entrance/exit edge multipoles"},
    {"EDGE1_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ccbend_example.edge1_multipoles), NULL, 0.0, 0, "input file for systematic entrance edge multipoles. Overrides EDGE_MULTIPOLES."},
    {"EDGE2_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ccbend_example.edge2_multipoles), NULL, 0.0, 0, "input file for systematic exit edge multipoles. Overrides EDGE_MULTIPOLES."},
    {"RANDOM_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ccbend_example.random_multipoles), NULL, 0.0, 0, "input file for random multipoles"},
    {"SYSTEMATIC_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.systematicMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply systematic and edge multipoles"},
    {"RANDOM_MULTIPOLE_FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&ccbend_example.randomMultipoleFactor), NULL, 1.0, 0, "Factor by which to multiply random multipoles"},
    {"REFERENCE_ORDER", "", IS_SHORT, 0, (long)((char *)&ccbend_example.referenceOrder), NULL, 0.0, 0, "Reference order for multipole errors. Overridden by value in multipole files, if those are given."},
    {"SYNCH_RAD", "", IS_SHORT, 0, (long)((char *)&ccbend_example.synch_rad), NULL, 0.0, 0, "include classical, single-particle synchrotron radiation?"},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&ccbend_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (quantum excitation)?"},
    {"ISR1PART", "", IS_SHORT, 0, (long)((char *)&ccbend_example.isr1Particle), NULL, 0.0, 1, "Include ISR for single-particle beam only if ISR=1 and ISR1PART=1"},
    {"USE_RAD_DIST", "", IS_SHORT, 0, (long)((char *)&ccbend_example.distributionBasedRadiation), NULL, 0.0, 0, "If nonzero, overrides SYNCH_RAD and ISR, causing simulation of radiation from distributions, optionally including opening angle."},
    {"ADD_OPENING_ANGLE", "", IS_SHORT, 0, (long)((char *)&ccbend_example.includeOpeningAngle), NULL, 0.0, 1, "If nonzero, radiation opening angle effects are added if USE_RAD_DIST is nonzero."},
    {"OPTIMIZE_FSE", "", IS_SHORT, 0, (long)((char *)&ccbend_example.optimizeFse), NULL, 0.0, 1, "Optimize strength (FSE) to obtain the ideal deflection angle."},
    {"OPTIMIZE_DX", "", IS_SHORT, 0, (long)((char *)&ccbend_example.optimizeDx), NULL, 0.0, 1, "Optimize x offset to obtain centered trajectory."},
    {"OPTIMIZE_FSE_ONCE", "", IS_SHORT, 0, (long)((char *)&ccbend_example.optimizeFseOnce), NULL, 0.0, 0, "If nonzero, the FSE offset is optimized only once, even if relevant parameters are changed."},
    {"OPTIMIZE_DX_ONCE", "", IS_SHORT, 0, (long)((char *)&ccbend_example.optimizeDxOnce), NULL, 0.0, 0, "If nonzero, the x offset is optimized only once, even if relevant parameters are changed."},
    {"COMPENSATE_KN", "", IS_SHORT, 0, (long)((char *)&ccbend_example.compensateKn), NULL, 0.0, 0, "If nonzero, K1 and K2 strengths are adjusted to compensate for the changes in FSE needed to center the trajectory."},
    {"EDGE_ORDER", "", IS_SHORT, 0, (long)((char *)&ccbend_example.edgeOrder), NULL, 0.0, 3, "Gives order of edge effects. Does not affect edge multipoles."},
    {"VERBOSE", "", IS_SHORT, 0, (long)((char *)&ccbend_example.verbose), NULL, 0.0, 0, "If nonzero, print messages showing optimized FSE and x offset."},
    };

BOFFAXE boffaxe_example;
PARAMETER boffaxe_param[N_BOFFAXE_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&boffaxe_example.length), NULL, 0.0, 0, "insertion length"},
    {"LFIELD", "M", IS_DOUBLE, 0, (long)((char *)&boffaxe_example.fieldLength), NULL, -1.0, 0, "expected length of the field map for verification purposes only."},
    {"FILENAME", NULL, IS_STRING, 0, (long)((char*)&boffaxe_example.filename), NULL, 0.0, 0, "name of file containing derivative data"},
    {"Z_COLUMN", NULL, IS_STRING, 0, (long)((char*)&boffaxe_example.zColumn), "z", 0.0, 0, "name of longitunidal coordinate column in the data file"},
    {"FIELD_COLUMN", NULL, IS_STRING, 0, (long)((char*)&boffaxe_example.fieldColumn), NULL, 0.0, 0, "name of derivative column in the data file"},
    {"ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&boffaxe_example.order), NULL, 0.0, 1, "order of transverse derivative"},
    {"EXPANSION_ORDER", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&boffaxe_example.expansionOrder), NULL, 0.0, 0, "order of expansion in x and y. If zero, determined by data in file."},
    {"STRENGTH", NULL, IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&boffaxe_example.strength), NULL, 1.0, 0, "factor by which to multiply field"},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&boffaxe_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&boffaxe_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&boffaxe_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&boffaxe_example.dz), NULL, 0.0, 0, "misalignment"},
    {"BX", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&boffaxe_example.Bx), NULL, 0.0, 0, "add BX*STRENGTH to Bx field"},
    {"BY", "T", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&boffaxe_example.By), NULL, 0.0, 0, "add BY*STRENGTH to By field"},
    {"Z_INTERVAL", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&boffaxe_example.zInterval), NULL, 0.0, 1, "input z data is sampled at this interval"},
    {"Z_SUBDIVISIONS", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&boffaxe_example.zSubdivisions), NULL, 0.0, 1, "Number of subdivisions of z interval to use in integration"},
    {"SYNCH_RAD", "", IS_SHORT, PARAM_CHANGES_MATRIX, (long)((char *)&boffaxe_example.synchRad), NULL, 0.0, 0, "if nonzero, include classical, single-particle synchrotron radiation"},
    {"ISR", "", IS_SHORT, 0, (long)((char *)&boffaxe_example.isr), NULL, 0.0, 0, "if nonzero, include incoherent synchrotron radiation (quantum excitation)"},
    {"PARTICLE_OUTPUT_FILE", "", IS_STRING, 0, (long)((char*)&boffaxe_example.particleOutputFile), NULL, 0.0, 0, "name of file for phase-space and field output. Use for debugging only!"},
};  

APCONTOUR apcontour_example;
/* aperture-contour physical parameters */
PARAMETER apcontour_param[N_APCONTOUR_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&apcontour_example.length), NULL, 0.0, 0, "length"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&apcontour_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&apcontour_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&apcontour_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 0, (long)((char *)&apcontour_example.dz), NULL, 0.0, 0, "misalignment"},
    {"RESOLUTION", "M", IS_DOUBLE, 0, (long)((char *)&apcontour_example.resolution), NULL, 0.0, 0, "z resolution of finding intersection"},
    {"INVERT", "", IS_LONG, 0, (long)((char *)&apcontour_example.invert), NULL, 0.0, 0, "if non-zero, contour defines an obstruction rather than an aperture"},
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&apcontour_example.filename), NULL, 0.0, 0, "name of file containing contour data"},
    {"XCOLUMN", "", IS_STRING, 0, (long)((char *)&apcontour_example.xColumn), NULL, 0.0, 0, "name of column containing x data"},
    {"YCOLUMN", "", IS_STRING, 0, (long)((char *)&apcontour_example.yColumn), NULL, 0.0, 0, "name of containing y data"},
    };

TAPERAPC taperapc_example;
/* tapered circular aperture physical parameters */
PARAMETER taperapc_param[N_TAPERAPC_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&taperapc_example.length), NULL, 0.0, 0, "length"},
    {"RSTART", "M", IS_DOUBLE, 0, (long)((char *)&taperapc_example.r[0]), NULL, 0.0, 0, "radius at the start"},
    {"REND", "M", IS_DOUBLE, 0, (long)((char *)&taperapc_example.r[1]), NULL, 0.0, 0, "radius at the end"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&taperapc_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&taperapc_example.dy), NULL, 0.0, 0, "misalignment"},
    {"STICKY", NULL, IS_SHORT, 0, (long)((char *)&taperapc_example.sticky), NULL, 0.0, 0, "final aperture holds downstream until next TAPERAPC, TAPERAPE, TAPERAPR, or MAXAMP"},
    };

TAPERAPE taperape_example;
/* tapered elliptical aperture physical parameters */
PARAMETER taperape_param[N_TAPERAPE_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&taperape_example.length), NULL, 0.0, 0, "length"},
    {"ASTART", "M", IS_DOUBLE, 0, (long)((char *)&taperape_example.a[0]), NULL, 0.0, 0, "horizontal semi-axis at the start"},
    {"AEND", "M", IS_DOUBLE, 0, (long)((char *)&taperape_example.a[1]), NULL, 0.0, 0, "horizontal semi-axis at the end"},
    {"BSTART", "M", IS_DOUBLE, 0, (long)((char *)&taperape_example.b[0]), NULL, 0.0, 0, "vertical semi-axis at the start"},
    {"BEND", "M", IS_DOUBLE, 0, (long)((char *)&taperape_example.b[1]), NULL, 0.0, 0, "vertical semi-axis at the end"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&taperape_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&taperape_example.dy), NULL, 0.0, 0, "misalignment"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&taperape_example.tilt), NULL, 0.0, 0, "misalignment"},
    {"RESOLUTION", "M", IS_DOUBLE, 0, (long)((char *)&taperape_example.resolution), NULL, 1e-6, 0, "z resolution of finding intersection"},
    {"XEXPONENT", NULL, IS_SHORT, 0, (long)((char *)&taperape_example.xExponent), NULL, 0, 2, "super-elliptical exponent (even number)"},
    {"YEXPONENT", NULL, IS_SHORT, 0, (long)((char *)&taperape_example.yExponent), NULL, 0, 2, "super-elliptical exponent (even number)"},
    {"STICKY", NULL, IS_SHORT, 0, (long)((char *)&taperape_example.sticky), NULL, 0.0, 0, "final aperture holds downstream until next TAPERAPC, TAPERAPE, TAPERAPR, or MAXAMP"},
    };

TAPERAPR taperapr_example;
/* tapered rectangular aperture physical parameters */
PARAMETER taperapr_param[N_TAPERAPR_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&taperapr_example.length), NULL, 0.0, 0, "length"},
    {"XSTART", "M", IS_DOUBLE, 0, (long)((char *)&taperapr_example.xmax[0]), NULL, 0.0, 0, "horizontal half-aperture at the start"},
    {"XEND", "M", IS_DOUBLE, 0, (long)((char *)&taperapr_example.xmax[1]), NULL, 0.0, 0, "horizontal half-aperture at the end"},
    {"YSTART", "M", IS_DOUBLE, 0, (long)((char *)&taperapr_example.ymax[0]), NULL, 0.0, 0, "vertical half-aperture at the start"},
    {"YEND", "M", IS_DOUBLE, 0, (long)((char *)&taperapr_example.ymax[1]), NULL, 0.0, 0, "vertical half-aperture at the end"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&taperapr_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&taperapr_example.dy), NULL, 0.0, 0, "misalignment"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&taperapr_example.tilt), NULL, 0.0, 0, "misalignment"},
    {"STICKY", NULL, IS_SHORT, 0, (long)((char *)&taperapr_example.sticky), NULL, 0.0, 0, "final aperture holds downstream until next TAPERAPC, TAPERAPE, TAPERAPR, or MAXAMP"},
    };

SHRFDF shrfdf_example;
/* names for space harmonic rf deflector parameters */
PARAMETER shrfdf_param[N_SHRFDF_PARAMS] = {
    {"FACTOR", "", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.factor), NULL, 1.0, 0, "A factor by which to multiply all components."},
    {"TILT", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"PERIOD_LENGTH", "M", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.period_length), NULL, 0.0, 0, "cavity period length, or cell length"},
    {"PERIOD_PHASE", "RAD", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.period_phase), NULL, 0.0, 0, "cavity period phase advance, or so-called working mode"},
    {"V0", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.v[0]), NULL, 0.0, 0, "effective voltage of space harmonic n=0"},
    {"V1", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.v[1]), NULL, 0.0, 0, "effective voltage of space harmonic n=1"},
    {"V2", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.v[2]), NULL, 0.0, 0, "effective voltage of space harmonic n=2"},
    {"V3", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.v[3]), NULL, 0.0, 0, "effective voltage of space harmonic n=3"},
    {"V4", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.v[4]), NULL, 0.0, 0, "effective voltage of space harmonic n=4"},
    {"V5", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.v[5]), NULL, 0.0, 0, "effective voltage of space harmonic n=5"},
    {"V6", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.v[6]), NULL, 0.0, 0, "effective voltage of space harmonic n=6"},
    {"V7", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.v[7]), NULL, 0.0, 0, "effective voltage of space harmonic n=7"},
    {"V8", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.v[8]), NULL, 0.0, 0, "effective voltage of space harmonic n=8"},
    {"V9", "V", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.v[9]), NULL, 0.0, 0, "effective voltage of space harmonic n=9"},
    {"PHASE0", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.phase[0]), NULL, 0.0, 0, "Phase of space harmonic n=0"},
    {"PHASE1", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.phase[1]), NULL, 0.0, 0, "Phase of space harmonic n=1"},
    {"PHASE2", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.phase[2]), NULL, 0.0, 0, "Phase of space harmonic n=2"},
    {"PHASE3", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.phase[3]), NULL, 0.0, 0, "Phase of space harmonic n=3"},
    {"PHASE4", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.phase[4]), NULL, 0.0, 0, "Phase of space harmonic n=4"},
    {"PHASE5", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.phase[5]), NULL, 0.0, 0, "Phase of space harmonic n=5"},
    {"PHASE6", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.phase[6]), NULL, 0.0, 0, "Phase of space harmonic n=6"},
    {"PHASE7", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.phase[7]), NULL, 0.0, 0, "Phase of space harmonic n=7"},
    {"PHASE8", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.phase[8]), NULL, 0.0, 0, "Phase of space harmonic n=8"},
    {"PHASE9", "HZ", IS_DOUBLE, PARAM_CHANGES_MATRIX, (long)((char *)&shrfdf_example.phase[9]), NULL, 0.0, 0, "Phase of space harmonic n=9"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&shrfdf_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    } ;

/* END OF ELEMENT DICTIONARY ARRAYS */

/* array of parameter structures */

#define MAT_LEN     HAS_MATRIX|HAS_LENGTH
#define MAT_LEN_NCAT HAS_MATRIX|HAS_LENGTH|DONT_CONCAT

ELEMENT_DESCRIPTION entity_description[N_TYPES] = {
    {                0,           0,                  0,    NULL           },
    {    N_QUAD_PARAMS,     MAT_LEN|BACKTRACK|DIVIDE_OK|IS_MAGNET|MATRIX_TRACKING|GPU_SUPPORT,       sizeof(QUAD),    quad_param     },
    {    N_BEND_PARAMS,     MAT_LEN|BACKTRACK|DIVIDE_OK|IS_MAGNET|MATRIX_TRACKING|GPU_SUPPORT,       sizeof(BEND),    sbend_param     },
    {    N_BEND_PARAMS,     MAT_LEN|BACKTRACK|DIVIDE_OK|IS_MAGNET|MATRIX_TRACKING|GPU_SUPPORT,       sizeof(BEND),    rbend_param     },
    {    N_DRIFT_PARAMS,     MAT_LEN|BACKTRACK|DIVIDE_OK|MATRIX_TRACKING|GPU_SUPPORT,      sizeof(DRIFT),    drift_param    }, 
    {    N_SEXT_PARAMS,     MAT_LEN|BACKTRACK|DIVIDE_OK|IS_MAGNET|MATRIX_TRACKING|GPU_SUPPORT,       sizeof(SEXT),    sext_param     },
    {    N_OCTU_PARAMS,     MAT_LEN|DIVIDE_OK|IS_MAGNET|MATRIX_TRACKING|GPU_SUPPORT,       sizeof(OCTU),    octu_param     },
    {    N_MULT_PARAMS,  MAT_LEN_NCAT|IS_MAGNET,       sizeof(MULT),    mult_param     }, 
    {    N_SOLE_PARAMS,     MAT_LEN|IS_MAGNET|MAT_CHW_ENERGY|DIVIDE_OK|GPU_SUPPORT|BACKTRACK,
           sizeof(SOLE),    sole_param     }, 
    {    N_HCOR_PARAMS,     MAT_LEN|IS_MAGNET|GPU_SUPPORT,        sizeof(HCOR),    hcor_param     }, 
    {    N_VCOR_PARAMS,     MAT_LEN|IS_MAGNET|GPU_SUPPORT,        sizeof(VCOR),    vcor_param     }, 
    {    N_RFCA_PARAMS,     MAT_LEN_NCAT|BACKTRACK|HAS_RF_MATRIX|MAY_CHANGE_ENERGY|MPALGORITHM|DIVIDE_OK|GPU_SUPPORT,
                     sizeof(RFCA),    rfca_param     }, 
    {                0,           NO_DICT_OUTPUT,                  0,    NULL           },
    {    N_HMON_PARAMS,     MAT_LEN_NCAT|BACKTRACK|MATRIX_TRACKING|GPU_SUPPORT|NO_APERTURE,       sizeof(HMON),    hmon_param     }, 
    {    N_VMON_PARAMS,     MAT_LEN_NCAT|BACKTRACK|MATRIX_TRACKING|GPU_SUPPORT|NO_APERTURE,       sizeof(VMON),    vmon_param     }, 
    {    N_MONI_PARAMS,     MAT_LEN_NCAT|BACKTRACK|MATRIX_TRACKING|GPU_SUPPORT|NO_APERTURE,       sizeof(MONI),    moni_param     }, 
    {    N_RCOL_PARAMS,  MAT_LEN_NCAT|BACKTRACK|GPU_SUPPORT,       sizeof(RCOL),    rcol_param     }, 
    {    N_ECOL_PARAMS,  MAT_LEN_NCAT|BACKTRACK|GPU_SUPPORT,       sizeof(ECOL),    ecol_param     }, 
    {    N_MARK_PARAMS,  0|MPALGORITHM|BACKTRACK|NO_APERTURE,       sizeof(MARK),    mark_param     }, 
    {    N_MATR_PARAMS,  MAT_LEN|HAS_RF_MATRIX|GPU_SUPPORT,  sizeof(MATR),    matr_param     }, 
    {    N_ALPH_PARAMS,  HAS_MATRIX|IS_MAGNET|MAT_CHW_ENERGY|GPU_SUPPORT,  sizeof(ALPH),    alph_param     }, 
    {    N_RFDF_PARAMS,  MAT_LEN_NCAT|HAS_RF_MATRIX|MPALGORITHM,       sizeof(RFDF),    rfdf_param     }, 
    {    N_RFTMEZ0_PARAMS,  MAT_LEN_NCAT|MAY_CHANGE_ENERGY|MPALGORITHM,    sizeof(RFTMEZ0),    rftmez0_param     }, 
    {    N_RMDF_PARAMS,  MAT_LEN_NCAT|UNIPROCESSOR,       sizeof(RMDF),    rmdf_param     }, 
    {    N_TMCF_PARAMS,  MAT_LEN_NCAT|MPALGORITHM,  sizeof(TMCF_MODE),    tmcf_param     }, 
    {    N_CEPL_PARAMS,  MAT_LEN_NCAT|MPALGORITHM,  sizeof(CE_PLATES),    cepl_param     }, 
    {   N_WATCH_PARAMS,  MPALGORITHM|BACKTRACK|RUN_ZERO_PARTICLES|NO_APERTURE,      sizeof(WATCH),    watch_param    }, 
    {    N_TWPL_PARAMS,  MAT_LEN_NCAT|MPALGORITHM,  sizeof(TW_PLATES),    twpl_param     }, 
    {  N_MALIGN_PARAMS,  HAS_MATRIX|BACKTRACK|DONT_CONCAT|GPU_SUPPORT,
                                         sizeof(MALIGN),    malign_param   },
    {  N_TWLA_PARAMS,  MAT_LEN_NCAT|MAY_CHANGE_ENERGY|HAS_RF_MATRIX|MPALGORITHM,   sizeof(TW_LINAC),    twla_param     },

    {  N_PEPPOT_PARAMS,  MAT_LEN_NCAT,     sizeof(PEPPOT),    peppot_param   },
    {  N_ENERGY_PARAMS,          MPALGORITHM|GPU_SUPPORT,     sizeof(ENERGY),    energy_param   },
    {  N_MAXAMP_PARAMS,           GPU_SUPPORT|BACKTRACK,     sizeof(MAXAMP),    maxamp_param   },
    {  N_ROTATE_PARAMS,  HAS_MATRIX|BACKTRACK|MATRIX_TRACKING|GPU_SUPPORT,     sizeof(ROTATE),    rotate_param   },
    { N_TRCOUNT_PARAMS,           UNIDIAGNOSTIC|NO_APERTURE,    sizeof(TRCOUNT),    trcount_param  },
    {  N_RECIRC_PARAMS,           0|NO_APERTURE,     sizeof(RECIRC),    recirc_param   },
    {  N_QFRING_PARAMS,     MAT_LEN|MATRIX_TRACKING,     sizeof(QFRING),    qfring_param   },
    { N_SCRAPER_PARAMS,  MAT_LEN_NCAT|BACKTRACK|GPU_SUPPORT,    sizeof(SCRAPER),    scraper_param  },
    {  N_CENTER_PARAMS,           0|GPU_SUPPORT,     sizeof(CENTER),    center_param   },
    {  N_KICKER_PARAMS,  MAT_LEN_NCAT|IS_MAGNET,     sizeof(KICKER),    kicker_param   },
    {   N_KSEXT_PARAMS, MAT_LEN_NCAT|IS_MAGNET|MAT_CHW_ENERGY|DIVIDE_OK|GPU_SUPPORT,      
                                          sizeof(KSEXT),    ksext_param    },
    {  N_KSBEND_PARAMS, MAT_LEN_NCAT|IS_MAGNET|NO_DICT_OUTPUT,
                                         sizeof(KSBEND),    ksbend_param   },
    {   N_KQUAD_PARAMS, MAT_LEN_NCAT|IS_MAGNET|MAT_CHW_ENERGY|DIVIDE_OK|GPU_SUPPORT, 
                                          sizeof(KQUAD),    kquad_param    },
    { N_MAGNIFY_PARAMS, HAS_MATRIX|MATRIX_TRACKING,     sizeof(MAGNIFY),    magnify_param  },
    {  N_SAMPLE_PARAMS,          0|NO_APERTURE,      sizeof(SAMPLE),    sample_param   },
    {   N_HVCOR_PARAMS,    MAT_LEN|IS_MAGNET,        sizeof(HVCOR),    hvcor_param    }, 
    { N_SCATTER_PARAMS,          0,     sizeof(SCATTER),    scatter_param  },
    {  N_NIBEND_PARAMS, MAT_LEN_NCAT|IS_MAGNET,
                                         sizeof(NIBEND),    nibend_param   },
    {   N_KPOLY_PARAMS,          0,       sizeof(KPOLY),    kpoly_param    }, 
    {  N_NISEPT_PARAMS, MAT_LEN_NCAT|IS_MAGNET,
                                         sizeof(NISEPT),    nisept_param   },
    {  N_RAMPRF_PARAMS, MAT_LEN_NCAT|HAS_RF_MATRIX|MAY_CHANGE_ENERGY|MPALGORITHM,    sizeof(RAMPRF),    ramprf_param   },
    {   N_RAMPP_PARAMS,  MAY_CHANGE_ENERGY|MPALGORITHM,       sizeof(RAMPP),    rampp_param    },
    {   N_STRAY_PARAMS,    MAT_LEN|MAT_CHW_ENERGY,
                                          sizeof(STRAY),    stray_param    },
    {  N_CSBEND_PARAMS, MAT_LEN_NCAT|IS_MAGNET|DIVIDE_OK|GPU_SUPPORT, 
                                        sizeof(CSBEND),    csbend_param   },
    {   N_TWMTA_PARAMS, MAT_LEN_NCAT|MPALGORITHM,     sizeof(TWMTA),    twmta_param    },
    {  N_MATTER_PARAMS,    MAT_LEN|GPU_SUPPORT,      sizeof(MATTER),   matter_param    },
    {  N_RFMODE_PARAMS,          MPALGORITHM,      sizeof(RFMODE),   rfmode_param    },
    { N_TRFMODE_PARAMS,          MPALGORITHM,     sizeof(TRFMODE),  trfmode_param    },
    { N_ZLONGIT_PARAMS,          0,     sizeof(ZLONGIT),  zlongit_param    },
    { N_SREFFECTS_PARAMS,        HAS_MATRIX|DONT_CONCAT,   sizeof(SREFFECTS),  sreffects_param  },
    { N_MODRF_PARAMS, MAT_LEN_NCAT|HAS_RF_MATRIX|MAY_CHANGE_ENERGY|MPALGORITHM,       sizeof(MODRF),    modrf_param     }, 
    { N_BMAPXY_PARAMS,     MAT_LEN_NCAT|IS_MAGNET,   sizeof(BMAPXY),  bmapxy_param      },
    { N_ZTRANSVERSE_PARAMS,      0,     sizeof(ZTRANSVERSE),  ztransverse_param    },
    { N_IBSCATTER_PARAMS,        MPALGORITHM,   sizeof(IBSCATTER),  ibscatter_param  },
    { N_FMULT_PARAMS,  MAT_LEN_NCAT|IS_MAGNET,       sizeof(FMULT),    fmult_param     }, 
    { N_WAKE_PARAMS, MAY_CHANGE_ENERGY|BACKTRACK|MPALGORITHM|GPU_SUPPORT, sizeof(WAKE), wake_param},
    { N_TRWAKE_PARAMS, 0|BACKTRACK|MPALGORITHM|GPU_SUPPORT, sizeof(TRWAKE), trwake_param},
    { N_TUBEND_PARAMS, 0, sizeof(TUBEND), tubend_param},
    { N_CHARGE_PARAMS, BACKTRACK|MPALGORITHM|NO_APERTURE, sizeof(CHARGE), charge_param},
    { N_PFILTER_PARAMS, MPALGORITHM|NO_APERTURE, sizeof(PFILTER), pfilter_param},
    { N_HISTOGRAM_PARAMS, RUN_ZERO_PARTICLES|BACKTRACK|MPALGORITHM|NO_APERTURE, sizeof(HISTOGRAM), histogram_param},
    {  N_CSRCSBEND_PARAMS, MAT_LEN_NCAT|IS_MAGNET|MPALGORITHM|GPU_SUPPORT,
       sizeof(CSRCSBEND),    csrcsbend_param   },
    {  N_CSRDRIFT_PARAMS, MAT_LEN_NCAT|MPALGORITHM|DIVIDE_OK|GPU_SUPPORT,
       sizeof(CSRDRIFT),    csrdrift_param   },
    {  N_RFCW_PARAMS,     MAT_LEN_NCAT|HAS_RF_MATRIX|MAY_CHANGE_ENERGY|MPALGORITHM|GPU_SUPPORT,
       sizeof(RFCW),    rfcw_param     }, 
    { N_REMCOR_PARAMS,           UNIPROCESSOR,     sizeof(REMCOR),    remcor_param   },
    { N_MAPSOLENOID_PARAMS,  MAT_LEN_NCAT,    sizeof(MAP_SOLENOID),    mapSolenoid_param    }, 
    { N_REFLECT_PARAMS,      HAS_MATRIX|MATRIX_TRACKING,    sizeof(REFLECT),    reflect_param  },
    { N_CLEAN_PARAMS,  0|NO_APERTURE, sizeof(CLEAN), clean_param },
    { N_TWISSELEMENT_PARAMS, HAS_MATRIX|DONT_CONCAT|MPALGORITHM,  sizeof(TWISSELEMENT), twissElement_param},
    { N_WIGGLER_PARAMS, MAT_LEN|MATRIX_TRACKING, sizeof(WIGGLER), wiggler_param},
    { N_SCRIPT_PARAMS,  MAT_LEN|DONT_CONCAT|MPALGORITHM, sizeof(SCRIPT),    script_param     }, 
    { N_FLOORELEMENT_PARAMS,  0|NO_APERTURE, sizeof(FLOORELEMENT),    floor_param     }, 
    { N_LTHINLENS_PARAMS,  HAS_MATRIX|IS_MAGNET|MATRIX_TRACKING,       sizeof(LTHINLENS),    lthinlens_param     },
    {  N_LMIRROR_PARAMS,  HAS_MATRIX|IS_MAGNET|MATRIX_TRACKING,       sizeof(LMIRROR),    lmirror_param     },
    {  N_EMATRIX_PARAMS,  HAS_MATRIX|HAS_RF_MATRIX|HAS_LENGTH|MAY_CHANGE_ENERGY|MPALGORITHM,  sizeof(EMATRIX),    ematrix_param     }, 
    {  N_FRFMODE_PARAMS,         MPALGORITHM,      sizeof(FRFMODE),   frfmode_param    },
    { N_FTRFMODE_PARAMS,         MPALGORITHM,     sizeof(FTRFMODE),  ftrfmode_param    },
    { N_TFBPICKUP_PARAMS,         MPALGORITHM|NO_APERTURE,      sizeof(TFBPICKUP),  tfbpickup_param    },
    { N_TFBDRIVER_PARAMS, MPALGORITHM|RUN_ZERO_PARTICLES,     sizeof(TFBDRIVER),  tfbdriver_param    },
    { N_LSCDRIFT_PARAMS, MAT_LEN_NCAT|MPALGORITHM|GPU_SUPPORT,     sizeof(LSCDRIFT),  lscdrift_param    },
    { N_DSCATTER_PARAMS,          0,     sizeof(DSCATTER),    dscatter_param  },
    { N_LSRMDLTR_PARAMS,    MAT_LEN_NCAT|MPALGORITHM, sizeof(LSRMDLTR), lsrMdltr_param },
    { N_TAYLORSERIES_PARAMS, MAT_LEN_NCAT|IS_MAGNET|NO_DICT_OUTPUT,    sizeof(TAYLORSERIES),  taylorSeries_param  },
    {    N_RFTM110_PARAMS,  0|MPALGORITHM,       sizeof(RFTM110),    rftm110_param     }, 
    {   N_CWIGGLER_PARAMS,  MAT_LEN_NCAT|IS_MAGNET, sizeof(CWIGGLER),    cwiggler_param     }, 
    {   N_EDRIFT_PARAMS, MAT_LEN|DIVIDE_OK|GPU_SUPPORT, sizeof(EDRIFT),    edrift_param   },
    {   N_SCMULT_PARAMS,    0,       sizeof(SCMULT),    scmult_param     },   
    {  N_ILMATRIX_PARAMS,  HAS_RF_MATRIX|MAT_LEN_NCAT,  sizeof(ILMATRIX),    ilmatrix_param     }, 
    {   N_TSCATTER_PARAMS,  0,       sizeof(TSCATTER),  tscatter_param     },   
    {   N_KQUSE_PARAMS, MAT_LEN_NCAT|IS_MAGNET|MAT_CHW_ENERGY|DIVIDE_OK,      
                                          sizeof(KQUSE),    kquse_param    },
    {   N_UKICKMAP_PARAMS, MAT_LEN_NCAT|IS_MAGNET|MPALGORITHM, sizeof(UKICKMAP),    ukickmap_param    },
    {  N_MKICKER_PARAMS,  MAT_LEN_NCAT|IS_MAGNET,     sizeof(MKICKER),    mkicker_param   },
    {  N_EMITTANCEELEMENT_PARAMS,  MPALGORITHM,    sizeof(EMITTANCEELEMENT),    emittanceElement_param   },
    { N_MHISTOGRAM_PARAMS, UNIPROCESSOR, sizeof(MHISTOGRAM), mhistogram_param},
    { N_FTABLE_PARAMS, MAT_LEN_NCAT|IS_MAGNET, sizeof(FTABLE), ftable_param},
    {   N_KOCT_PARAMS, MAT_LEN_NCAT|IS_MAGNET|MAT_CHW_ENERGY|DIVIDE_OK,      
                                          sizeof(KOCT),    koct_param    },
    {N_MRADITEGRALS_PARAMS, 0, sizeof(MRADINTEGRALS), mRadIntegrals_param },
    { N_APPLE_PARAMS,  MAT_LEN_NCAT|IS_MAGNET, sizeof(APPLE),    apple_param}, 
    { N_MRFDF_PARAMS,  MPALGORITHM,   sizeof(MRFDF),    mrfdf_param     }, 
    { N_CORGPIPE_PARAMS, MAY_CHANGE_ENERGY|MPALGORITHM|MAT_LEN_NCAT, sizeof(CORGPIPE), corgpipe_param},
    { N_LRWAKE_PARAMS, MPALGORITHM, sizeof(LRWAKE), lrwake_param},
    {    N_EHCOR_PARAMS,     MAT_LEN_NCAT|IS_MAGNET,        sizeof(EHCOR),    ehcor_param     }, 
    {    N_EVCOR_PARAMS,     MAT_LEN_NCAT|IS_MAGNET,        sizeof(EVCOR),    evcor_param     }, 
    {    N_EHVCOR_PARAMS,    MAT_LEN_NCAT|IS_MAGNET,        sizeof(EHVCOR),   ehvcor_param     }, 
    { N_BMAPXYZ_PARAMS,    MAT_LEN_NCAT|IS_MAGNET,   sizeof(BMAPXYZ),  bmapxyz_param      },
    { N_BRAT_PARAMS,     MAT_LEN_NCAT|IS_MAGNET,   sizeof(BRAT),  brat_param      },
    { N_BGGEXP_PARAMS,   MAT_LEN_NCAT|IS_MAGNET,   sizeof(BGGEXP),  bggexp_param      },
    { N_BRANCH_PARAMS,   0, sizeof(BRANCH),  branch_param      },
    { N_IONEFFECTS_PARAMS, MPALGORITHM,     sizeof(IONEFFECTS),  ionEffects_param    },
    { N_SLICE_POINT_PARAMS, MPALGORITHM|RUN_ZERO_PARTICLES|NO_APERTURE, sizeof(SLICE_POINT),    slice_point_param   }, 
    { N_SPEEDBUMP_PARAMS, MAT_LEN_NCAT|BACKTRACK, sizeof(SPEEDBUMP),    speedbump_param   }, 
    { N_CCBEND_PARAMS, MAT_LEN_NCAT, sizeof(CCBEND),    ccbend_param   }, 
    { N_HKPOLY_PARAMS, MAT_LEN_NCAT, sizeof(HKPOLY),    hkpoly_param   },
    { N_BOFFAXE_PARAMS,  MAT_LEN_NCAT|IS_MAGNET,   sizeof(BOFFAXE),  boffaxe_param  },
    { N_APCONTOUR_PARAMS, MAT_LEN_NCAT, sizeof(APCONTOUR), apcontour_param},
    { N_TAPERAPC_PARAMS, MAT_LEN_NCAT, sizeof(TAPERAPC), taperapc_param},
    { N_TAPERAPE_PARAMS, MAT_LEN_NCAT, sizeof(TAPERAPE), taperape_param},
    { N_TAPERAPR_PARAMS, MAT_LEN_NCAT, sizeof(TAPERAPR), taperapr_param},
    { N_SHRFDF_PARAMS,  MPALGORITHM,   sizeof(SHRFDF),    shrfdf_param     },
} ;

void compute_offsets()
{
  long i, j, typeSize=0;
  for (i=0; i<N_TYPES; i++) {
    for (j=entity_description[i].n_params-1; j>=0; j--)
      entity_description[i].parameter[j].offset -= entity_description[i].parameter[0].offset;
    if ((j=entity_description[i].n_params-1)>=0) {
      switch (entity_description[i].parameter[j].type) {
      case IS_DOUBLE:
	typeSize = sizeof(double);
	break;
      case IS_LONG:
	typeSize = sizeof(long);
	break;
      case IS_SHORT:
	typeSize = sizeof(short);
	break;
      case IS_STRING:
	typeSize = sizeof(char*);
	break;
      default:
	fprintf(stderr, "Error: invalid item type code %ld for %s parameter of %s\n",
		entity_description[i].parameter[j].type, 
		entity_description[i].parameter[j].name,
		entity_name[i]);
	exit(1);
	break;
      }
      entity_description[i].user_structure_size = entity_description[i].parameter[j].offset + typeSize;
    }
  }

  for (i=0; i<N_TYPES; i++) {
    size_t difference;
    if (entity_description[i].n_params==0) continue;
    if (entity_description[i].parameter[0].offset<0) {
      printf("error: bad initial parameter offset for element type %s\n", entity_name[i]);
      fflush(stdout);
      exitElegant(1);
    }
    for (j=1; j<entity_description[i].n_params; j++) {
      /* difference of offsets must be positive and less than size of double */
      if ((difference=(entity_description[i].parameter[j].offset-entity_description[i].parameter[j-1].offset))<=0) {
        printf("error: bad parameter offset (retrograde) for element type %s, parameter %s\n",
               entity_name[i], entity_description[i].parameter[j].name?entity_description[i].parameter[j].name:"NULL");
        fflush(stdout);
        exitElegant(1);
      }
      if (difference>sizeof(double) && i!=T_TWISSELEMENT && i!=T_EMATRIX) {
        printf("error: bad parameter offset (too large %ld) for element type %s, parameter %s\n",
               difference, entity_name[i], entity_description[i].parameter[j].name?entity_description[i].parameter[j].name:"NULL");
        fflush(stdout);
        exitElegant(1);
      }
    }
    entity_description[i].flags |= OFFSETS_CHECKED;
  }
}


/* The sigma matrix s[i][j] is stored in a 21-element array.  These indices give the i and j values 
 * corresponding to an element the array.  We have i<=j (upper triangular).  Values are filled in
 * by setSigmaIndices, which is called on start-up.
 */
long sigmaIndex1[21], sigmaIndex2[21];

/* This array gives the index in the 21-element array for given i and j.  Values are filled in
 * by setSigmaIndices, which is called on start-up.
 */
long sigmaIndex3[6][6];

#if USE_MPI
char *mpiAbortDescription[N_MPI_ABORT_TYPES] = {
  "-",
  "Bunch too long in ZLONGIT element",
  "Bunch too long in ZTRANSVERSE element",
  "Bunch too long in RFMODE element",
  "Bucket assignment error",
  "Error in particle pointer order",
  "Non-positive particle ID value",
  "RF cavity fiducialization failed",
  "SREFFECTS setup error: EYREF and COUPLING both nonzero",
};
#endif

char *chamberShapeChoice[N_CHAMBER_SHAPES] = {
  "round", "rect.", "ellip.", "sup.ellip.", "?",
};

