/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: track_data4.c
 * contents: data arrays for elegant 
 *
 * Michael Borland, 1989-1993
 */
#include "mdb.h"
#include "track.h"

/* flag word for trace mode */
long trace_mode = 0;

long particleID = 1;

char *entity_name[N_TYPES] = {
    "LINE", "QUAD", "SBEN", "RBEN", "DRIF", "SEXT", "OCTU", "MULT", "SOLE", 
    "HKICK", "VKICK", "RFCA", "ELSE", "HMON", "VMON", "MONI", "RCOL", "ECOL", 
    "MARK", "MATR", "ALPH", "RFDF", "RFTM", "RMDF", "TMCF", "CEPL", "WATCH",
    "TWPL", "MALIGN", "TWLA", "PEPPOT", "ENERGY", "MAXAMP", "ROTATE",
    "TRCOUNT", "RECIRC", "QUFRINGE", "SCRAPER", "CENTER", "BUMPER",
    "KSEXT", "KSBEND", "KQUAD", "MAGNIFY", "SAMPLE", "KICKER", "SCATTER",
    "NIBEND", "KPOLY", "NISEPT", "RAMPRF", "RAMPP", "STRAY", "CSBEND",
    "TWMTA", "MATTER", "RFMODE", "TRFMODE", "ZLONGIT", "SREFFECTS",
    "MODRF", "BMAPXY", "ZTRANSVERSE", "IBSCATTER", "FMULT",
    "WAKE", "TRWAKE", "TUBEND",
    };

char *madcom_name[N_MADCOMS] = {
    "", "USE", "TITLE", "RETURN", "RENAME"
        } ;

char *entity_text[N_TYPES] = {
    NULL, 
    "A quadrupole implemented as a matrix, up to 2nd order.",
    "A sector dipole implemented as a matrix, up to 2nd order.",
    "A rectangular dipole, implemented as a SBEND with edge angles.",
    "A drift space implemented as a matrix, up to 2nd order",
    "A sextupole implemented as a matrix, up to 2nd order",
    "Not implemented--use the MULT element.",
    "A canonical kick multipole.",
    "A solenoid implemented as a matrix, up to 2nd order.",
    "A horizontal steering dipole implemented as a matrix, up to 2nd order.",
    "A vertical steering dipole implemented as a matrix, up to 2nd order.",
    "A first-order matrix RF cavity with exact phase dependence.",
    "Not implemented.",
    "A horizontal position monitor, accepting a rpn equation for the readout as a\n\
function of the actual position (x).",
    "A vertical position monitor, accepting a rpn equation for the readout as a\n\
function of the actual position (y).",
    "A two-plane position monitor, accepting two rpn equations for the readouts\n\
as a function of the actual positions (x and y).",
    "A rectangular collimator.",
    "An elliptical collimator.",
    "A marker, equivalent to a zero-length drift space.",
    "Explicit matrix input from a text file, in the format written by the print_matrix\n\
command.",
    "An alpha magnet implemented as a matrix, up to 3rd order.  PART is used to split\n\
the magnet into halves.  XSN and DPN allow momentum filtration at the midpoint.",
    "A deflecting TM RF cavity, using an approximate analytical solution.",
    "Not implemented.",
    "A linearly-ramped electric field deflector, using an approximate analytical solution.",
    "A numerically-integrated accelerating TM RF cavity with spatially-constant fields.",
    "A numerically-integrated linearly-ramped electric field deflector.",
    "A beam property/motion monitor--allowed modes are centroid, coordinate, and fft.\n\
Output is in SDDS format.",
    "A numerically-integrated traveling-wave stripline deflector.",
    "A misalignment of the beam, implemented as a zero-order matrix.",
    "A numerically-integrated first-space-harmonic traveling-wave linear accelerator.",
    "A pepper-pot plate.",
    "An element that matches the central momentum to the beam momentum, or changes\n\
the central momentum or energy to a specified value.",
    "A collimating element that sets the maximum transmitted particle amplitudes for\n\
all following elements, until the next MAXAMP.",
    "An element that rotates the beam coordinates about the longitudinal axis.",
    "An element that defines the point from which transmission calculations are made.",
    "An element that defines the point to which particles recirculate in multi-pass\n\
tracking",
    "An element consisting of a linearly increasing or decreasing quadrupole field.",
    "A collimating element that sticks into the beam from one side only.  The\n\
directions 0, 1, 2, and 3 are from +x, +y, -x, and -y, respectively.",
    "An element that centers the beam transversely on the ideal trajectory.",
    "A time-dependent uniform-field rectangular kicker magnet with no fringe effects.\n\
The waveform is in mpl format, with time in seconds and amplitude normalized to 1.",
    "A canonical kick sextupole, which differs from the MULT element with ORDER=2 in\n\
that it can be used for chromaticity correction.",
    "A kick bending magnet which is NOT canonical, but is better than a 2nd order\n\
matrix implementation.",
    "A canonical kick quadrupole, which differs from the MULT element with ORDER=1 in\n\
that it can be used for tune correction.",
    "An element that allows multiplication of phase-space coordinates of all particles\n\
by constants.",
    "An element that reduces the number of particles in the beam by interval-based or\n\
random sampling.",
    "A combined horizontal-vertical steering magnet implemented as a matrix, up to\n\
2nd order.",
    "A scattering element to add gaussian random numbers to particle coordinates.",
    "A numerically-integrated dipole magnet with various extended-fringe-field models.",
    "A thin kick element with polynomial dependence on the coordinates in one plane.",
    "A numerically-integrated dipole magnet with a Cartesian gradient.",
    "A voltage-ramped RF cavity, implemented like RFCA.  The voltage ramp pattern is\n\
given by a mpl-format file of the voltage factor vs time in seconds.",
    "A momentum-ramping element that changes the central momentum according to a mpl\n\
format file of the momentum factor vs time in seconds.",
    "A stray field element with local and global components.  Global components are\n\
defined relative to the initial beamline direction.",
    "A canonical kick sector dipole magnet.  INTEGRATION_ORDER may be 2 or 4.",
    "A numerically-integrated traveling-wave muffin-tin accelerator.",
    "A Coulomb-scattering and energy-absorbing element simulating material in the\n\
beam path.",
    "A simulation of a beam-driven TM monopole mode of a RF cavity.",
    "A simulation of a beam-driven TM dipole mode of a RF cavity.",
    "A simulation of a single-pass broad-band or functionally specified longitudinal\n\
impedance.",
    "Simulation of synchrotron radiation effects (damping and quantum excitation).",
    "A first-order matrix RF cavity with exact phase dependence, plus optional amplitude\n\
and phase modulation.",
    "A map of Bx and By vs x and y.",
    "A simulation of a single-pass broad-band or functionally-specified transverse dipole impedance.",
    "A simulation of intra-beam scattering.",
    "Multipole kick element with coefficient input from an SDDS file.",
    "Longitudinal wake specified as a function of time lag behind the particle.",
    "Transverse wake specified as a function of time lag behind the particle.",
    "A special rectangular bend element for top-up backtracking.",
    } ;

QUAD quad_example;
/* quadrupole physical parameters */
PARAMETER quad_param[N_QUAD_PARAMS]={
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&quad_example.length), NULL, 0.0, 0},
    {"K1", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&quad_example.k1), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&quad_example.tilt), NULL, 0.0, 0},
    {"FFRINGE", "", IS_DOUBLE, 1, (long)((char *)&quad_example.ffringe), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&quad_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&quad_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&quad_example.dz), NULL, 0.0, 0},
    {"FSE", "M", IS_DOUBLE, 1, (long)((char *)&quad_example.fse), NULL, 0.0, 0},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&quad_example.order), NULL, 0.0, 0}
    };

BEND bend_example;
/* bending magnet physical parameters */
PARAMETER bend_param[N_BEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&bend_example.length), NULL, 0.0, 0},
    {"ANGLE", "RAD", IS_DOUBLE, 1, (long)((char *)&bend_example.angle), NULL, 0.0, 0},
    {"K1", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&bend_example.k1), NULL, 0.0, 0},
    {"E1", "RAD", IS_DOUBLE, 1, (long)((char *)&bend_example.e1), NULL, 0.0, 0},
    {"E2", "RAD", IS_DOUBLE, 1, (long)((char *)&bend_example.e2), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&bend_example.tilt), NULL, 0.0, 0},
    {"K2", "1/M$a3$n", IS_DOUBLE, 1, (long)((char *)&bend_example.k2), NULL, 0.0, 0},
    {"H1", "1/M", IS_DOUBLE, 1, (long)((char *)&bend_example.h1), NULL, 0.0, 0},
    {"H2", "1/M", IS_DOUBLE, 1, (long)((char *)&bend_example.h2), NULL, 0.0, 0},
    {"HGAP", "M", IS_DOUBLE, 1, (long)((char *)&bend_example.hgap), NULL, 0.0, 0},
    {"FINT", "", IS_DOUBLE, 1, (long)((char *)&bend_example.fint), NULL, DEFAULT_FINT, 0},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&bend_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&bend_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&bend_example.dz), NULL, 0.0, 0},
    {"FSE", "", IS_DOUBLE, 1, (long)((char *)&bend_example.fse), NULL, 0.0, 0},
    {"ETILT", "", IS_DOUBLE, 1, (long)((char *)&bend_example.etilt), NULL, 0.0, 0},
    {"EDGE1_EFFECTS", "", IS_LONG, 1, (long)((char *)&bend_example.edge1_effects), NULL, 0.0, 1},
    {"EDGE2_EFFECTS", "", IS_LONG, 1, (long)((char *)&bend_example.edge2_effects), NULL, 0.0, 1},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&bend_example.order), NULL, 0.0, 0},
    {"EDGE_ORDER", "", IS_LONG, 1, (long)((char *)&bend_example.edge_order), NULL, 0.0, 0},
    {"TRANSPORT", "", IS_LONG, 1, (long)((char *)&bend_example.TRANSPORT), NULL, 0.0, 0}
    };

DRIFT drift_example;
/* drift length physical parameters */
PARAMETER drift_param[N_DRIFT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&drift_example.length), NULL, 0.0, 0},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&drift_example.order), NULL, 0.0, 0}
    };

SEXT sext_example;
/* sextupole physical parameters */
PARAMETER sext_param[N_SEXT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&sext_example.length), NULL, 0.0, 0},
    {"K2", "1/M$a3$n", IS_DOUBLE, 1, (long)((char *)&sext_example.k2), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&sext_example.tilt), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&sext_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&sext_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&sext_example.dz), NULL, 0.0, 0},
    {"FSE", "M", IS_DOUBLE, 1, (long)((char *)&sext_example.fse), NULL, 0.0, 0},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&sext_example.order), NULL, 0.0, 0}
    };
   
MULT mult_example;
/* multipole physical parameters */
PARAMETER mult_param[N_MULT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.length), NULL, 0.0, 0},
    {"KNL", "M$a(1-ORDER)$n", IS_DOUBLE, 0, (long)((char *)&mult_example.KnL), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&mult_example.tilt), NULL, 0.0, 0},
    {"BORE", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.bore), NULL, 0.0, 0},
    {"BNL", "T M", IS_DOUBLE, 0, (long)((char *)&mult_example.BnL), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.dz), NULL, 0.0, 0},
    {"FACTOR", "", IS_DOUBLE, 0, (long)((char *)&mult_example.factor), NULL, 1.0, 0},
    {"ORDER", "", IS_LONG, 0, (long)((char *)&mult_example.order), NULL, 0.0, 1},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&mult_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&mult_example.synch_rad), NULL, 0.0, 0},
    };

FMULT fmult_example;
/* multipole physical parameters */
PARAMETER fmult_param[N_FMULT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&fmult_example.length), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&fmult_example.tilt), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&fmult_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&fmult_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 0, (long)((char *)&fmult_example.dz), NULL, 0.0, 0},
    {"FSE", "", IS_DOUBLE, 0, (long)((char *)&fmult_example.fse), NULL, 0.0, 0},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&fmult_example.n_kicks), NULL, 0.0, 1},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&fmult_example.synch_rad), NULL, 0.0, 0},
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&fmult_example.filename), NULL, 0.0, 0},
    };

SOLE sole_example;
/* solenoid physical parameters */
PARAMETER sole_param[N_SOLE_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&sole_example.length), NULL, 0.0, 0},
    {"KS", "RAD/M", IS_DOUBLE, 1, (long)((char *)&sole_example.ks), NULL, 0.0, 0},
    {"B", "T", IS_DOUBLE, 1, (long)((char *)&sole_example.B), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&sole_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&sole_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&sole_example.dz), NULL, 0.0, 0},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&sole_example.order), NULL, 0.0, 0},
    };
   
HCOR hcor_example;
/* horizontal corrector physical parameters */
PARAMETER hcor_param[N_HCOR_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&hcor_example.length), NULL, 0.0, 0},
    {"KICK", "RAD", IS_DOUBLE, 1, (long)((char *)&hcor_example.kick), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&hcor_example.tilt), NULL, 0.0, 0},
    {"B2", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&hcor_example.b2), NULL, 0.0, 0},
    {"CALIBRATION", "", IS_DOUBLE, 1, (long)((char *)&hcor_example.calibration), NULL, 1.0, 0},
    {"EDGE_EFFECTS", "", IS_LONG, 0, (long)((char *)&hcor_example.edge_effects), NULL, 0.0, 0},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&hcor_example.order), NULL, 0.0, 0},
    {"STEERING", "", IS_LONG, 0, (long)((char *)&hcor_example.steering), NULL, 0.0, 1},
    };

VCOR vcor_example;
/* vertical corrector physical parameters */
PARAMETER vcor_param[N_VCOR_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&vcor_example.length), NULL, 0.0, 0},
    {"KICK", "RAD", IS_DOUBLE, 1, (long)((char *)&vcor_example.kick), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&vcor_example.tilt), NULL, 0.0, 0},
    {"B2", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&vcor_example.b2), NULL, 0.0, 0},
    {"CALIBRATION", "", IS_DOUBLE, 1, (long)((char *)&vcor_example.calibration), NULL, 1.0, 0},
    {"EDGE_EFFECTS", "", IS_LONG, 0, (long)((char *)&vcor_example.edge_effects), NULL, 0.0, 0},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&vcor_example.order), NULL, 0.0, 0},
    {"STEERING", "", IS_LONG, 0, (long)((char *)&vcor_example.steering), NULL, 0.0, 1},
    };

RFCA rfca_example;
/* rf cavity physical parameters */
PARAMETER rfca_param[N_RFCA_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&rfca_example.length), NULL, 0.0, 0},
    {"VOLT", "V", IS_DOUBLE, 0, (long)((char *)&rfca_example.volt), NULL, 0.0, 0},
    {"PHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&rfca_example.phase), NULL, 0.0, 0},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&rfca_example.freq), NULL, 500.0e6, 0},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&rfca_example.Q), NULL, 0.0, 0},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&rfca_example.phase_reference), NULL, 0.0, 0},
    {"CHANGE_P0", "", IS_LONG, 0, (long)((char *)&rfca_example.change_p0), NULL, 0.0, 0}, 
    {"CHANGE_T", "", IS_LONG, 0, (long)((char *)&rfca_example.change_t), NULL, 0.0, 0}, 
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&rfca_example.fiducial), NULL, 0.0, 0}
    };
   
HMON hmon_example;
/* name for horizontal monitor physical parameters */
PARAMETER hmon_param[N_HMON_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&hmon_example.length), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&hmon_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&hmon_example.dy), NULL, 0.0, 0},
    {"WEIGHT", "", IS_DOUBLE, 0, (long)((char *)&hmon_example.weight), NULL, 1.0, 0},
    {"TILT", "", IS_DOUBLE, 0, (long)((char *)&hmon_example.tilt), NULL, 0.0, 0},
    {"CALIBRATION", "", IS_DOUBLE, 0, (long)((char *)&hmon_example.calibration), NULL, 1.0, 0},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&hmon_example.order), NULL, 1.0, 0},
    {"READOUT", "", IS_STRING, 0, (long)((char *)&hmon_example.readout), NULL, 0.0, 1}
    } ;
   
VMON vmon_example;
/* name for vertical monitor physical parameters */
PARAMETER vmon_param[N_VMON_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&vmon_example.length), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&vmon_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&vmon_example.dy), NULL, 0.0, 0},
    {"WEIGHT", "", IS_DOUBLE, 0, (long)((char *)&vmon_example.weight), NULL, 1.0, 0},
    {"TILT", "", IS_DOUBLE, 0, (long)((char *)&vmon_example.tilt), NULL, 0.0, 0},
    {"CALIBRATION", "", IS_DOUBLE, 0, (long)((char *)&vmon_example.calibration), NULL, 1.0, 0},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&vmon_example.order), NULL, 1.0, 0},
    {"READOUT", "", IS_STRING, 0, (long)((char *)&vmon_example.readout), NULL, 0.0, 1}
    } ;

MONI moni_example;   
/* name for two-plane monitor physical parameters */
PARAMETER moni_param[N_MONI_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&moni_example.length), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&moni_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&moni_example.dy), NULL, 0.0, 0},
    {"WEIGHT", "", IS_DOUBLE, 0, (long)((char *)&moni_example.weight), NULL, 1.0, 0},
    {"TILT", "", IS_DOUBLE, 0, (long)((char *)&moni_example.tilt), NULL, 0.0, 0},
    {"XCALIBRATION", "", IS_DOUBLE, 0, (long)((char *)&moni_example.xcalibration), NULL, 1.0, 0},
    {"YCALIBRATION", "", IS_DOUBLE, 0, (long)((char *)&moni_example.ycalibration), NULL, 1.0, 0},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&moni_example.order), NULL, 1.0, 0},
    {"XREADOUT", "", IS_STRING, 0, (long)((char *)&moni_example.x_readout), NULL, 0.0, 1},
    {"YREADOUT", "", IS_STRING, 0, (long)((char *)&moni_example.y_readout), NULL, 0.0, 1}
    } ;

RCOL rcol_example;   
/* name for rectangular collimator physical parameters */
PARAMETER rcol_param[N_RCOL_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&rcol_example.length), NULL, 0.0, 0},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&rcol_example.x_max), NULL, 0.0, 0},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&rcol_example.y_max), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&rcol_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&rcol_example.dy), NULL, 0.0, 0},
    } ;
   
ECOL ecol_example;
/* name for elliptical collimator physical parameters */
PARAMETER ecol_param[N_ECOL_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&ecol_example.length), NULL, 0.0, 0},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&ecol_example.x_max), NULL, 0.0, 0},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&ecol_example.y_max), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&ecol_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&ecol_example.dy), NULL, 0.0, 0},
    } ;

MARK mark_example;
/* name for marker parameters */
PARAMETER mark_param[N_MARK_PARAMS] = {
    {"FITPOINT", "", IS_LONG, 0, (long)((char *)&mark_example.fitpoint), NULL, 0.0, 0},
    } ;

MATR matr_example;
/* name for matrix parameters */
PARAMETER matr_param[N_MATR_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&matr_example.length), NULL, 0.0, 0},
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&matr_example.filename), "", 0.0, 1},
    {"ORDER", "", IS_LONG, 0, (long)((char *)&matr_example.order), NULL, 0.0, 1},
    } ;

ALPH alph_example;
/* names for alpha magnet parameters */
PARAMETER alph_param[N_ALPH_PARAMS] = {
    {"XMAX", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.xmax), NULL, 0.0, 0},
    {"XS1", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.xs1), NULL, 0.0, 0},
    {"XS2", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.xs2), NULL, 0.0, 0},
    {"DP1", "", IS_DOUBLE, 1, (long)((char *)&alph_example.dp1), NULL, -1, 0},
    {"DP2", "", IS_DOUBLE, 1, (long)((char *)&alph_example.dp2), NULL, 1, 0},
    {"XPUCK", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.xPuck), NULL, -1, 0},
    {"WIDTHPUCK", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.widthPuck), NULL, 0, 0},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.dz), NULL, 0.0, 0},
    {"PART", "", IS_LONG, 0, (long)((char *)&alph_example.part), NULL, 0.0, 0},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&alph_example.order), NULL, 0.0, 0}
    } ;

RFDF rfdf_example;
/* names for rf deflector parameters */
PARAMETER rfdf_param[N_RFDF_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&rfdf_example.length), NULL, 0.0, 0},
    {"PHASE", "RAD", IS_DOUBLE, 0, (long)((char *)&rfdf_example.phase), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&rfdf_example.tilt), NULL, 0.0, 0},
    {"FREQUENCY", "HZ", IS_DOUBLE, 0, (long)((char *)&rfdf_example.frequency), NULL, DEFAULT_FREQUENCY, 0},
    {"VOLTAGE", "V", IS_DOUBLE, 0, (long)((char *)&rfdf_example.voltage), NULL, 0.0, 0},
    {"GAP", "M", IS_DOUBLE, 0, (long)((char *)&rfdf_example.gap), NULL, DEFAULT_GAP, 0},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&rfdf_example.time_offset), NULL, 0.0, 0},
    {"B_FIELD", "T", IS_DOUBLE, 0, (long)((char *)&rfdf_example.B_field), NULL, 0.0, 0},
    {"N_SECTIONS", "", IS_LONG, 0, (long)((char *)&rfdf_example.n_sections), NULL, 0.0, DEFAULT_N_SECTIONS},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&rfdf_example.phase_reference), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&rfdf_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&rfdf_example.dy), NULL, 0.0, 0},
    } ;

TM_MODE rftm_example;
/* names for tm mode cavity (from superfish and using nag integrator) 
 * parameters 
 */
PARAMETER rftm_param[N_RFTM_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&rftm_example.length), NULL, 0.0, 0},
    {"FREQUENCY", "HZ", IS_DOUBLE, 0, (long)((char *)&rftm_example.frequency), NULL, DEFAULT_FREQUENCY, 0},
    {"PHASE", "RAD", IS_DOUBLE, 0, (long)((char *)&rftm_example.phase), NULL, 0.0, 0},
    {"EZ_PEAK", "V", IS_DOUBLE, 0, (long)((char *)&rftm_example.Ez_peak), NULL, 0.0, 0},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&rftm_example.time_offset), NULL, 0.0, 0},
    {"RADIAL_OFFSET", "M", IS_DOUBLE, 0, (long)((char *)&rftm_example.radial_offset), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&rftm_example.tilt), NULL, 0.0, 0},
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&rftm_example.accuracy), NULL, DEFAULT_ACCURACY, 0},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&rftm_example.phase_reference), NULL, 0.0, 0},
    {"N_STEPS", "", IS_LONG, 0, (long)((char *)&rftm_example.n_steps), NULL, 0.0, 100},
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&rftm_example.filename), "", 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&rftm_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&rftm_example.dy), NULL, 0.0, 0},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&rftm_example.method), DEFAULT_INTEG_METHOD, 0.0, 0},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&rftm_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0}
    } ;

RMDF rmdf_example;
/* names for rf deflector parameters */
PARAMETER rmdf_param[N_RMDF_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&rmdf_example.length), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&rmdf_example.tilt), NULL, 0.0, 0},
    {"RAMP_TIME", "S", IS_DOUBLE, 0, (long)((char *)&rmdf_example.ramp_time), NULL, DEFAULT_RAMP_TIME, 0},
    {"VOLTAGE", "V", IS_DOUBLE, 0, (long)((char *)&rmdf_example.voltage), NULL, 0.0, 0},
    {"GAP", "M", IS_DOUBLE, 0, (long)((char *)&rmdf_example.gap), NULL, DEFAULT_GAP, 0},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&rmdf_example.time_offset), NULL, 0.0, 0},
    {"N_SECTIONS", "", IS_LONG, 0, (long)((char *)&rmdf_example.n_sections), NULL, 0.0, 10},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&rmdf_example.phase_reference), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&rmdf_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&rmdf_example.dy), NULL, 0.0, 0},
    } ;

TMCF_MODE tmcf_example;
/* names for tm mode cavity with spatially constant fields, using nag 
 * integrator parameters 
 */
PARAMETER tmcf_param[N_TMCF_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&tmcf_example.length), NULL, 0.0, 0},
    {"FREQUENCY", "HZ", IS_DOUBLE, 0, (long)((char *)&tmcf_example.frequency), NULL, DEFAULT_FREQUENCY, 0},
    {"PHASE", "S", IS_DOUBLE, 0, (long)((char *)&tmcf_example.phase), NULL, 0.0, 0},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&tmcf_example.time_offset), NULL, 0.0, 0},
    {"RADIAL_OFFSET", "M", IS_DOUBLE, 0, (long)((char *)&tmcf_example.radial_offset), NULL, DEFAULT_RADIAL_OFFSET, 0},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&tmcf_example.tilt), NULL, 0.0, 0},
    {"ER", "V", IS_DOUBLE, 0, (long)((char *)&tmcf_example.Er), NULL, 0.0, 0},
    {"BPHI", "T", IS_DOUBLE, 0, (long)((char *)&tmcf_example.Bphi), NULL, 0.0, 0},
    {"EZ", "V", IS_DOUBLE, 0, (long)((char *)&tmcf_example.Ez), NULL, 0.0, 0},
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&tmcf_example.accuracy), NULL, DEFAULT_ACCURACY, 0},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&tmcf_example.x_max), NULL, 0.0, 0},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&tmcf_example.y_max), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&tmcf_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&tmcf_example.dy), NULL, 0.0, 0},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&tmcf_example.phase_reference), NULL, 0.0, 0},
    {"N_STEPS", "", IS_LONG, 0, (long)((char *)&tmcf_example.n_steps), NULL, 0.0, 100},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&tmcf_example.method), DEFAULT_INTEG_METHOD, 0.0, 0},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&tmcf_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0}
    } ;

CE_PLATES cepl_example;
/* names for ramped deflector using nag integrator parameters 
 */
PARAMETER cepl_param[N_CEPL_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&cepl_example.length), NULL, 0.0, 0},
    {"RAMP_TIME", "S", IS_DOUBLE, 0, (long)((char *)&cepl_example.ramp_time), NULL, DEFAULT_RAMP_TIME, 0},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&cepl_example.time_offset), NULL, 0.0, 0},
    {"VOLTAGE", "V", IS_DOUBLE, 0, (long)((char *)&cepl_example.voltage), NULL, 0.0, 0},
    {"GAP", "M", IS_DOUBLE, 0, (long)((char *)&cepl_example.gap), NULL, DEFAULT_GAP, 0},
    {"STATIC_VOLTAGE", "V", IS_DOUBLE, 0, (long)((char *)&cepl_example.static_voltage), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&cepl_example.tilt), NULL, 0.0, 0},
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&cepl_example.accuracy), NULL, DEFAULT_ACCURACY, 0},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&cepl_example.x_max), NULL, 0.0, 0},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&cepl_example.y_max), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&cepl_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&cepl_example.dy), NULL, 0.0, 0},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&cepl_example.phase_reference), NULL, 0.0, 0},
    {"N_STEPS", "", IS_LONG, 0, (long)((char *)&cepl_example.n_steps), NULL, 0.0, 100},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&cepl_example.method), DEFAULT_INTEG_METHOD, 0.0, 0},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&cepl_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0}
    } ;

WATCH watch_example;
/* names for watch point */
PARAMETER watch_param[N_WATCH_PARAMS] = {
    {"FRACTION", "", IS_DOUBLE, 0, (long)((char *)&watch_example.fraction), NULL, 1.0, 0},
    {"INTERVAL", "", IS_LONG, 0, (long)((char *)&watch_example.interval), NULL, 0.0, 1},
    {"START_PASS", "", IS_LONG, 0, (long)((char*)&watch_example.start_pass), NULL, 0.0, 0},
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&watch_example.filename), "", 0.0, 0},
    {"LABEL", "", IS_STRING, 0, (long)((char *)&watch_example.label), "", 0.0, 0},
    {"MODE", "", IS_STRING, 0, (long)((char *)&watch_example.mode), "coordinates", 0.0, 0},
    {"X_DATA", "", IS_LONG, 0, (long)((char*)&watch_example.xData), NULL, 0.0, 1},
    {"Y_DATA", "", IS_LONG, 0, (long)((char*)&watch_example.yData), NULL, 0.0, 1},
    {"LONGIT_DATA", "", IS_LONG, 0, (long)((char*)&watch_example.longitData), NULL, 0.0, 1},
    } ;

TW_PLATES twpl_example;
/* names for tw (tem) ramped deflector using nag integrator parameters 
 */
PARAMETER twpl_param[N_TWPL_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&twpl_example.length), NULL, 0.0, 0},
    {"RAMP_TIME", "S", IS_DOUBLE, 0, (long)((char *)&twpl_example.ramp_time), NULL, DEFAULT_RAMP_TIME, 0},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&twpl_example.time_offset), NULL, 0.0, 0},
    {"VOLTAGE", "V", IS_DOUBLE, 0, (long)((char *)&twpl_example.voltage), NULL, 0.0, 0},
    {"GAP", "M", IS_DOUBLE, 0, (long)((char *)&twpl_example.gap), NULL, DEFAULT_GAP, 0},
    {"STATIC_VOLTAGE", "V", IS_DOUBLE, 0, (long)((char *)&twpl_example.static_voltage), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&twpl_example.tilt), NULL, 0.0, 0},
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&twpl_example.accuracy), NULL, DEFAULT_ACCURACY, 0},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twpl_example.x_max), NULL, 0.0, 0},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twpl_example.y_max), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&twpl_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&twpl_example.dy), NULL, 0.0, 0},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&twpl_example.phase_reference), NULL, 0.0, 0},
    {"N_STEPS", "", IS_LONG, 0, (long)((char *)&twpl_example.n_steps), NULL, 0.0, 100},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&twpl_example.method), DEFAULT_INTEG_METHOD, 0.0, 0},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&twpl_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0}
    } ;

MALIGN malign_example;
/* names for misaligment parameters */
PARAMETER malign_param[N_MALIGN_PARAMS] = {
    {"DXP", "", IS_DOUBLE, 1, (long)((char *)&malign_example.dxp), NULL, 0.0, 0},
    {"DYP", "", IS_DOUBLE, 1, (long)((char *)&malign_example.dyp), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&malign_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&malign_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&malign_example.dz), NULL, 0.0, 0},
    {"DT", "S", IS_DOUBLE, 1, (long)((char *)&malign_example.dt), NULL, 0.0, 0},
    {"DP", "", IS_DOUBLE, 1, (long)((char *)&malign_example.dp), NULL, 0.0, 0},
    {"DE", "", IS_DOUBLE, 1, (long)((char *)&malign_example.de), NULL, 0.0, 0},
    {"ON_PASS", "", IS_LONG, 0, (long)((char *)&malign_example.on_pass), NULL, 0.0, -1},
    } ;

TW_LINAC twla_example;
/* names for traveling-wave linac parameters
 */
PARAMETER twla_param[N_TWLA_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&twla_example.length), NULL, 0.0, 0},
    {"FREQUENCY", "HZ", IS_DOUBLE, 0, (long)((char *)&twla_example.frequency), NULL, DEFAULT_FREQUENCY, 0},
    {"PHASE", "RAD", IS_DOUBLE, 0, (long)((char *)&twla_example.phase), NULL, 0.0, 0},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&twla_example.time_offset), NULL, 0.0, 0},
    {"EZ", "V/M", IS_DOUBLE, 0, (long)((char *)&twla_example.Ez), NULL, 0.0, 0},
    {"B_SOLENOID", "T", IS_DOUBLE, 0, (long)((char *)&twla_example.B_solenoid), NULL, 0.0, 0},
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&twla_example.accuracy), NULL, DEFAULT_ACCURACY, 0},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twla_example.x_max), NULL, 0.0, 0},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twla_example.y_max), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&twla_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&twla_example.dy), NULL, 0.0, 0},
    {"BETA_WAVE", "", IS_DOUBLE, 0, (long)((char *)&twla_example.beta_wave), NULL, DEFAULT_BETA_WAVE, 0},
    {"ALPHA", "1/M", IS_DOUBLE, 0, (long)((char *)&twla_example.alpha), NULL, 0.0, 0},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&twla_example.phase_reference), NULL, 0.0, 0},
    {"N_STEPS", "", IS_LONG, 0, (long)((char *)&twla_example.n_steps), NULL, 0.0, 100},
    {"FOCUSSING", "", IS_LONG, 0, (long)((char *)&twla_example.focussing), NULL, 0.0, 1},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&twla_example.method), DEFAULT_INTEG_METHOD, 0.0, 0},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&twla_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0}
    } ;

PEPPOT peppot_example;
/* names for pepper-pot plate */
PARAMETER peppot_param[N_PEPPOT_PARAMS] = {
    {"L"           , "M", IS_DOUBLE, 0, (long)((char *)&peppot_example.length), NULL, 0.0, 0},
    {"RADII"       , "M", IS_DOUBLE, 0, (long)((char *)&peppot_example.radii), NULL, 0.0, 0},
    {"TRANSMISSION", "", IS_DOUBLE, 0, (long)((char *)&peppot_example.transmission), NULL, 0.0, 0},
    {"TILT"        , "RAD", IS_DOUBLE, 0, (long)((char *)&peppot_example.tilt), NULL, 0.0, 0},
    {"THETA_RMS"   , "RAD", IS_DOUBLE, 0, (long)((char *)&peppot_example.theta_rms), NULL, 0.0, 0},
    {"N_HOLES"     , "", IS_LONG  , 0, (long)((char *)&peppot_example.n_holes), NULL, 0.0, 0},
    } ;

ENERGY energy_example;
/* names for energy */
PARAMETER energy_param[N_ENERGY_PARAMS] = {
    {"CENTRAL_ENERGY", "MC$a2$n", IS_DOUBLE, 0, (long)((char *)&energy_example.central_energy), NULL, 0.0, 0},
    {"CENTRAL_MOMENTUM", "MC", IS_DOUBLE, 0, (long)((char *)&energy_example.central_momentum), NULL, 0.0, 0},
    {"MATCH_BEAMLINE", "", IS_LONG, 0, (long)((char *)&energy_example.match_beamline), NULL, 0.0, 0},
    {"MATCH_PARTICLES", "", IS_LONG, 0, (long)((char *)&energy_example.match_particles), NULL, 0.0, 0},
    } ;

MAXAMP maxamp_example;
/* names for max.amp. */
PARAMETER maxamp_param[N_MAXAMP_PARAMS] = {
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&maxamp_example.x_max), NULL, 0.0, 0},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&maxamp_example.y_max), NULL, 0.0, 0},
    {"ELLIPTICAL", "", IS_LONG, 0, (long)((char *)&maxamp_example.elliptical), NULL, 0.0, 0}
    } ;

ROTATE rotate_example;
/* names for beam rotation */
PARAMETER rotate_param[N_ROTATE_PARAMS] = {
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&rotate_example.tilt), NULL, 0.0, 0},
    } ;

/* names for transmission count */
PARAMETER trcount_param[N_TRCOUNT_PARAMS] = {
    {"DUMMY", "", IS_LONG, 0, 0, NULL, 0.0, 0},
    } ;

/* names for recirculation point */
PARAMETER recirc_param[N_RECIRC_PARAMS] = {
    {"I_RECIRC_ELEMENT", "", IS_LONG, 0, 0, NULL, 0.0, 0},
    } ;

QFRING qfring_example;
/* quadrupole fringe-field physical parameters */
PARAMETER qfring_param[N_QFRING_PARAMS]={
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&qfring_example.length), NULL, 0.0, 0},
    {"K1", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&qfring_example.k1), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&qfring_example.tilt), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&qfring_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&qfring_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&qfring_example.dz), NULL, 0.0, 0},
    {"FSE", "M", IS_DOUBLE, 1, (long)((char *)&qfring_example.fse), NULL, 0.0, 0},
    {"DIRECTION", "", IS_LONG, 1, (long)((char *)&qfring_example.direction), NULL, 0.0, 0},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&qfring_example.order), NULL, 0.0, 0}
    };

SCRAPER scraper_example;
/* scraper physical parameters */
PARAMETER scraper_param[N_SCRAPER_PARAMS]={
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&scraper_example.length), NULL, 0.0, 0},
    {"POSITION", "M", IS_DOUBLE, 0, (long)((char *)&scraper_example.position), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&scraper_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&scraper_example.dy), NULL, 0.0, 0},
    {"XO", "M", IS_DOUBLE, 0, (long)((char *)&scraper_example.Xo), NULL, 0.0, 0},
    {"INSERT_FROM", "", IS_STRING, 0, (long)((char *)&scraper_example.insert_from), NULL, 0.0, 0},
    {"ELASTIC", "", IS_LONG,  0, (long)((char *)&scraper_example.elastic), NULL, 0.0, 0},
    {"DIRECTION", "", IS_LONG,  0, (long)((char *)&scraper_example.direction), NULL, 0.0, 0},
    };

CENTER center_example;
/* beam centering physical parameters */
PARAMETER center_param[N_CENTER_PARAMS]={
    {"X" , "M", IS_LONG, 0, (long)((char *)&center_example.x), NULL, 0.0, 1},
    {"XP", "", IS_LONG, 0, (long)((char *)&center_example.xp), NULL, 0.0, 1},
    {"Y" , "M", IS_LONG, 0, (long)((char *)&center_example.y), NULL, 0.0, 1},
    {"YP", "", IS_LONG, 0, (long)((char *)&center_example.yp), NULL, 0.0, 1},
    };

KICKER kicker_example;
/* kicker physical parameters */
PARAMETER kicker_param[N_KICKER_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&kicker_example.length), NULL, 0.0, 0},
    {"ANGLE", "RAD", IS_DOUBLE, 0, (long)((char *)&kicker_example.angle), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&kicker_example.tilt), NULL, 0.0, 0},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&kicker_example.time_offset), NULL, 0.0, 0},
    {"PERIODIC", "", IS_LONG, 0, (long)((char *)&kicker_example.periodic), NULL, 0.0, 0},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&kicker_example.phase_reference), NULL, 0.0, 0},
    {"FIRE_ON_PASS", "", IS_LONG, 0, (long)((char *)&kicker_example.fire_on_pass), NULL, 0.0, 0},
    {"WAVEFORM", "", IS_STRING, 0, (long)((char *)&kicker_example.waveform), NULL, 0.0, 0},
    {"SPATIAL_DEPENDENCE", "", IS_STRING, 0, (long)((char *)&kicker_example.spatial_dependence), NULL, 0.0, 0},
    } ;

KSEXT ksext_example;
/* kick sextupole physical parameters */
PARAMETER ksext_param[N_KSEXT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&ksext_example.length), NULL, 0.0, 0},
    {"K2", "1/M$a3$n", IS_DOUBLE, 1, (long)((char *)&ksext_example.k2), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&ksext_example.tilt), NULL, 0.0, 0},
    {"BORE", "M", IS_DOUBLE, 1, (long)((char *)&ksext_example.bore), NULL, 0.0, 0},
    {"B", "T", IS_DOUBLE, 1, (long)((char *)&ksext_example.B), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&ksext_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&ksext_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&ksext_example.dz), NULL, 0.0, 0},
    {"FSE", "M", IS_DOUBLE, 1, (long)((char *)&ksext_example.fse), NULL, 0.0, 0},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&ksext_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&ksext_example.synch_rad), NULL, 0.0, 0},
    {"SYSTEMATIC_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ksext_example.systematic_multipoles), NULL, 0.0, 0},
    {"RANDOM_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ksext_example.random_multipoles), NULL, 0.0, 0},
    {"INTEGRATION_ORDER", "", IS_LONG, 0, (long)((char *)&ksext_example.integration_order), NULL, 0.0, 4},
    };

KSBEND ksbend_example;
/* symplectic sector bending magnet physical parameters */
PARAMETER ksbend_param[N_KSBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.length), NULL, 0.0, 0},
    {"ANGLE", "RAD", IS_DOUBLE, 1, (long)((char *)&ksbend_example.angle), NULL, 0.0, 0},
    {"K1", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&ksbend_example.k1), NULL, 0.0, 0},
    {"K2", "1/M$a3$n", IS_DOUBLE, 1, (long)((char *)&ksbend_example.k2), NULL, 0.0, 0},
    {"K3", "1/M$a3$n", IS_DOUBLE, 0, (long)((char *)&ksbend_example.k3), NULL, 0.0, 0},
    {"K4", "1/M$a4$n", IS_DOUBLE, 0, (long)((char *)&ksbend_example.k4), NULL, 0.0, 0},
    {"E1", "RAD", IS_DOUBLE, 1, (long)((char *)&ksbend_example.e1), NULL, 0.0, 0},
    {"E2", "RAD", IS_DOUBLE, 1, (long)((char *)&ksbend_example.e2), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&ksbend_example.tilt), NULL, 0.0, 0},
    {"H1", "1/M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.h1), NULL, 0.0, 0},
    {"H2", "1/M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.h2), NULL, 0.0, 0},
    {"HGAP", "M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.hgap), NULL, 0.0, 0},
    {"FINT", "", IS_DOUBLE, 1, (long)((char *)&ksbend_example.fint), NULL, DEFAULT_FINT, 0},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.dz), NULL, 0.0, 0},
    {"FSE", "", IS_DOUBLE, 1, (long)((char *)&ksbend_example.fse), NULL, 0.0, 0},
    {"ETILT", "", IS_DOUBLE, 1, (long)((char *)&ksbend_example.etilt), NULL, 0.0, 0},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&ksbend_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS},
    {"NONLINEAR", "", IS_LONG, 0, (long)((char *)&ksbend_example.nonlinear), NULL, 0.0, 1},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&ksbend_example.synch_rad), NULL, 0.0, 0},
    {"EDGE1_EFFECTS", "", IS_LONG, 1, (long)((char *)&ksbend_example.edge1_effects), NULL, 0.0, 1},
    {"EDGE2_EFFECTS", "", IS_LONG, 1, (long)((char *)&ksbend_example.edge2_effects), NULL, 0.0, 1},
    {"EDGE_ORDER", "", IS_LONG, 1, (long)((char *)&ksbend_example.edge_order), NULL, 0.0, 1},
    {"PARAXIAL", "", IS_LONG, 0, (long)((char *)&ksbend_example.paraxial), NULL, 0.0, 0},
    {"TRANSPORT", "", IS_LONG, 1, (long)((char *)&ksbend_example.TRANSPORT), NULL, 0.0, 0},
    {"METHOD", "", IS_STRING, 0, (long)((char *)&ksbend_example.method), "modified-midpoint", 0.0, 0}
    };

KQUAD kquad_example;
/* kick quadrupole physical parameters */
PARAMETER kquad_param[N_KQUAD_PARAMS]={
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&kquad_example.length), NULL, 0.0, 0},
    {"K1", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&kquad_example.k1), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&kquad_example.tilt), NULL, 0.0, 0},
    {"BORE", "M", IS_DOUBLE, 1, (long)((char *)&kquad_example.bore), NULL, 0.0, 0},
    {"B", "T", IS_DOUBLE, 1, (long)((char *)&kquad_example.B), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&kquad_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&kquad_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&kquad_example.dz), NULL, 0.0, 0},
    {"FSE", "M", IS_DOUBLE, 1, (long)((char *)&kquad_example.fse), NULL, 0.0, 0},
    {"N_KICKS", "", IS_LONG, 1, (long)((char *)&kquad_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&kquad_example.synch_rad), NULL, 0.0, 0},
    {"SYSTEMATIC_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&kquad_example.systematic_multipoles), NULL, 0.0, 0},
    {"RANDOM_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&kquad_example.random_multipoles), NULL, 0.0, 0},
    {"INTEGRATION_ORDER", "", IS_LONG, 0, (long)((char *)&kquad_example.integration_order), NULL, 0.0, 4},
    };

MAGNIFY magnify_example;
/* magnifier physical parameters */
PARAMETER magnify_param[N_MAGNIFY_PARAMS] = {
    {"MX", "", IS_DOUBLE, 1, (long)((char *)&magnify_example.mx), NULL, 1.0, 0},
    {"MXP", "", IS_DOUBLE, 1, (long)((char *)&magnify_example.mxp), NULL, 1.0, 0},
    {"MY", "", IS_DOUBLE, 1, (long)((char *)&magnify_example.my), NULL, 1.0, 0},
    {"MYP", "", IS_DOUBLE, 1, (long)((char *)&magnify_example.myp), NULL, 1.0, 0},
    {"MS", "", IS_DOUBLE, 1, (long)((char *)&magnify_example.ms), NULL, 1.0, 0},
    {"MDP", "", IS_DOUBLE, 1, (long)((char *)&magnify_example.mdp), NULL, 1.0, 0}
    } ;
    
SAMPLE sample_example;
/* sample physical parameters */
PARAMETER sample_param[N_SAMPLE_PARAMS] = {
    {"FRACTION", "", IS_DOUBLE, 0, (long)((char *)&sample_example.fraction), NULL, 1.0, 0},
    {"INTERVAL", "", IS_LONG, 0, (long)((char *)&sample_example.interval), NULL, 1.0, 1},
    } ;
    
HVCOR hvcor_example;
/* horizontal/vertical corrector physical parameters */
PARAMETER hvcor_param[N_HVCOR_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&hvcor_example.length), NULL, 0.0, 0},
    {"HKICK", "RAD", IS_DOUBLE, 1, (long)((char *)&hvcor_example.xkick), NULL, 0.0, 0},
    {"VKICK", "RAD", IS_DOUBLE, 1, (long)((char *)&hvcor_example.ykick), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&hvcor_example.tilt), NULL, 0.0, 0},
    {"B2", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&hvcor_example.b2), NULL, 0.0, 0},
    {"HCALIBRATION", "", IS_DOUBLE, 1, (long)((char *)&hvcor_example.xcalibration), NULL, 1.0, 0},
    {"VCALIBRATION", "", IS_DOUBLE, 1, (long)((char *)&hvcor_example.ycalibration), NULL, 1.0, 0},
    {"EDGE_EFFECTS", "", IS_LONG, 0, (long)((char *)&hvcor_example.edge_effects), NULL, 0.0, 0},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&hvcor_example.order), NULL, 0.0, 0},
    {"STEERING", "", IS_LONG, 0, (long)((char *)&hvcor_example.steering), NULL, 0.0, 1},
    };

SCATTER scatter_example;
/* scatter physical parameters */
PARAMETER scatter_param[N_SCATTER_PARAMS] = {
    {"X", "M", IS_DOUBLE, 0, (long)((char*)&scatter_example.x), NULL, 0.0, 0},
    {"XP", "M", IS_DOUBLE, 0, (long)((char*)&scatter_example.xp), NULL, 0.0, 0},
    {"Y", "M", IS_DOUBLE, 0, (long)((char*)&scatter_example.y), NULL, 0.0, 0},
    {"YP", "M", IS_DOUBLE, 0, (long)((char*)&scatter_example.yp), NULL, 0.0, 0},
    {"DP", "M", IS_DOUBLE, 0, (long)((char*)&scatter_example.dp), NULL, 0.0, 0},
    } ;
    
NIBEND nibend_example;
/* integrated bending magnet physical parameters */
PARAMETER nibend_param[N_NIBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.length), NULL, 0.0, 0},
    {"ANGLE", "RAD", IS_DOUBLE, 1, (long)((char *)&nibend_example.angle), NULL, 0.0, 0},
    {"E1", "RAD", IS_DOUBLE, 1, (long)((char *)&nibend_example.e1), NULL, 0.0, 0},
    {"E2", "RAD", IS_DOUBLE, 1, (long)((char *)&nibend_example.e2), NULL, 0.0, 0},
    {"TILT", "", IS_DOUBLE, 1, (long)((char *)&nibend_example.tilt), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.dz), NULL, 0.0, 0},
    {"FINT", "", IS_DOUBLE, 1, (long)((char *)&nibend_example.fint), NULL, DEFAULT_FINT, 0},
    {"HGAP", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.hgap), NULL, 0.0, 0},
    {"FP1", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.fp1), NULL, 10.0, 0},
    {"FP2", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.fp2), NULL, 1.0, 0},
    {"FSE", "", IS_DOUBLE, 1, (long)((char *)&nibend_example.fse), NULL, 0.0, 0},
    {"ETILT", "", IS_DOUBLE, 1, (long)((char *)&nibend_example.etilt), NULL, 0.0, 0},
    {"ACCURACY", "", IS_DOUBLE, 1, (long)((char *)&nibend_example.accuracy), NULL, DEFAULT_ACCURACY, 0},
    {"MODEL", "", IS_STRING, 0, (long)((char *)&nibend_example.model), DEFAULT_NIBEND_TYPE, 0.0, 0},
    {"METHOD", "", IS_STRING, 0, (long)((char *)&nibend_example.method), DEFAULT_INTEG_METHOD, 0.0, 0},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&nibend_example.synch_rad), NULL, 0.0, 0},
    };

KPOLY kpoly_example;
/* kick-polynomial physical parameters */
PARAMETER kpoly_param[N_KPOLY_PARAMS] = {
    {"COEFFICIENT", "M$A-ORDER$N", IS_DOUBLE, 0, (long)((char *)&kpoly_example.coefficient), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&kpoly_example.tilt), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&kpoly_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&kpoly_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 0, (long)((char *)&kpoly_example.dz), NULL, 0.0, 0},
    {"FACTOR", "", IS_DOUBLE, 0, (long)((char *)&kpoly_example.factor), NULL, 1.0, 0},
    {"ORDER", "", IS_LONG, 0, (long)((char *)&kpoly_example.order), NULL, 0.0, 0},
    {"PLANE", "", IS_STRING, 0, (long)((char *)&kpoly_example.plane), "x", 0.0, 0}
    };

NISEPT nisept_example;
/* integrated septum magnet physical parameters */
PARAMETER nisept_param[N_NISEPT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&nisept_example.length), NULL, 0.0, 0},
    {"ANGLE", "RAD", IS_DOUBLE, 1, (long)((char *)&nisept_example.angle), NULL, 0.0, 0},
    {"E1", "RAD", IS_DOUBLE, 1, (long)((char *)&nisept_example.e1), NULL, 0.0, 0},
    {"B1", "M", IS_DOUBLE, 1, (long)((char *)&nisept_example.b1), NULL, 0.0, 0},
    {"Q1REF", "M", IS_DOUBLE, 1, (long)((char *)&nisept_example.q1_ref), NULL, 0.0, 0},
    {"FLEN", "M", IS_DOUBLE, 1, (long)((char *)&nisept_example.flen), NULL, 0.0, 0},
    {"ACCURACY", "", IS_DOUBLE, 1, (long)((char *)&nisept_example.accuracy), NULL, DEFAULT_ACCURACY, 0},
    {"METHOD", "", IS_STRING, 0, (long)((char *)&nisept_example.method), DEFAULT_INTEG_METHOD, 0.0, 0},
    {"MODEL", "", IS_STRING, 0, (long)((char *)&nisept_example.model), DEFAULT_NIBEND_TYPE, 0.0, 0},
    };

RAMPRF ramprf_example;
/* ramped rf cavity physical parameters */
PARAMETER ramprf_param[N_RAMPRF_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&ramprf_example.length), NULL, 0.0, 0},
    {"VOLT", "V", IS_DOUBLE, 0, (long)((char *)&ramprf_example.volt), NULL, 0.0, 0},
    {"PHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&ramprf_example.phase), NULL, 0.0, 0},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&ramprf_example.freq), NULL, 500.0e6, 0},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&ramprf_example.phase_reference), NULL, 0.0, 0},
    {"VOLT_WAVEFORM", "", IS_STRING, 0, (long)((char *)&ramprf_example.vwaveform), NULL, 0.0, 0},
    {"PHASE_WAVEFORM", "", IS_STRING, 0, (long)((char *)&ramprf_example.pwaveform), NULL, 0.0, 0},
    {"FREQ_WAVEFORM", "", IS_STRING, 0, (long)((char *)&ramprf_example.fwaveform), NULL, 0.0, 0},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&ramprf_example.fiducial), NULL, 0.0, 0},
    };

/* momentum ramp physical parameters */
PARAMETER rampp_param[N_RAMPP_PARAMS] = {
    {"WAVEFORM", "", IS_STRING, 0, 0, NULL, 0.0, 0}
    };

STRAY stray_example;
/* stray field physical parameters */
PARAMETER stray_param[N_STRAY_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&stray_example.length), NULL, 0.0, 0},
    {"LBX", "T", IS_DOUBLE, 1, (long)((char *)&stray_example.lBx), NULL, 0.0, 0},
    {"LBY", "T", IS_DOUBLE, 1, (long)((char *)&stray_example.lBy), NULL, 0.0, 0},
    {"GBX", "T", IS_DOUBLE, 1, (long)((char *)&stray_example.gBx), NULL, 0.0, 0},
    {"GBY", "T", IS_DOUBLE, 1, (long)((char *)&stray_example.gBy), NULL, 0.0, 0},
    {"GBZ", "T", IS_DOUBLE, 1, (long)((char *)&stray_example.gBz), NULL, 0.0, 0},
    {"ORDER", "", IS_LONG, 0, (long)((char *)&stray_example.order), NULL, 0.0, 0},
    };

CSBEND csbend_example;
/* canonically-integrated sector bending magnet physical parameters */
PARAMETER csbend_param[N_CSBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&csbend_example.length), NULL, 0.0, 0},
    {"ANGLE", "RAD", IS_DOUBLE, 1, (long)((char *)&csbend_example.angle), NULL, 0.0, 0},
    {"K1", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&csbend_example.k1), NULL, 0.0, 0},
    {"K2", "1/M$a3$n", IS_DOUBLE, 1, (long)((char *)&csbend_example.k2), NULL, 0.0, 0},
    {"K3", "1/M$a3$n", IS_DOUBLE, 0, (long)((char *)&csbend_example.k3), NULL, 0.0, 0},
    {"K4", "1/M$a4$n", IS_DOUBLE, 0, (long)((char *)&csbend_example.k4), NULL, 0.0, 0},
    {"E1", "RAD", IS_DOUBLE, 1, (long)((char *)&csbend_example.e1), NULL, 0.0, 0},
    {"E2", "RAD", IS_DOUBLE, 1, (long)((char *)&csbend_example.e2), NULL, 0.0, 0},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&csbend_example.tilt), NULL, 0.0, 0},
    {"H1", "1/M", IS_DOUBLE, 1, (long)((char *)&csbend_example.h1), NULL, 0.0, 0},
    {"H2", "1/M", IS_DOUBLE, 1, (long)((char *)&csbend_example.h2), NULL, 0.0, 0},
    {"HGAP", "M", IS_DOUBLE, 1, (long)((char *)&csbend_example.hgap), NULL, 0.0, 0},
    {"FINT", "", IS_DOUBLE, 1, (long)((char *)&csbend_example.fint), NULL, DEFAULT_FINT, 0},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&csbend_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&csbend_example.dy), NULL, 0.0, 0},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&csbend_example.dz), NULL, 0.0, 0},
    {"FSE", "", IS_DOUBLE, 1, (long)((char *)&csbend_example.fse), NULL, 0.0, 0},
    {"ETILT", "", IS_DOUBLE, 1, (long)((char *)&csbend_example.etilt), NULL, 0.0, 0},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&csbend_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS},
    {"NONLINEAR", "", IS_LONG, 0, (long)((char *)&csbend_example.nonlinear), NULL, 0.0, 1},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&csbend_example.synch_rad), NULL, 0.0, 0},
    {"EDGE1_EFFECTS", "", IS_LONG, 0, (long)((char *)&csbend_example.edge1_effects), NULL, 0.0, 1},
    {"EDGE2_EFFECTS", "", IS_LONG, 0, (long)((char *)&csbend_example.edge2_effects), NULL, 0.0, 1},
    {"INTEGRATION_ORDER", "", IS_LONG, 0, (long)((char *)&csbend_example.integration_order), NULL, 0.0, 2},
    {"EDGE1_KICK_LIMIT", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.edge1_kick_limit), NULL, -1, 0},
    {"EDGE2_KICK_LIMIT", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.edge2_kick_limit), NULL, -1, 0},
    {"KICK_LIMIT_SCALING", "", IS_LONG, 0, (long)((char *)&csbend_example.kick_limit_scaling), NULL, 0, 0},
    };

TUBEND tubend_example;
/* special bending magnet for top-up with entry through the side! */
PARAMETER tubend_param[N_TUBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&tubend_example.length), NULL, 0.0, 0},
    {"ANGLE", "RAD", IS_DOUBLE, 1, (long)((char *)&tubend_example.angle), NULL, 0.0, 0},
    {"FSE", "", IS_DOUBLE, 1, (long)((char *)&tubend_example.fse), NULL, 0.0, 0},
    {"OFFSET", "", IS_DOUBLE, 1, (long)((char *)&tubend_example.offset), NULL, 0.0, 0},
    {"MAGNET_WIDTH", "", IS_DOUBLE, 1, (long)((char *)&tubend_example.magnet_width), NULL, 0.0, 0},
    {"MAGNET_ANGLE", "", IS_DOUBLE, 1, (long)((char *)&tubend_example.magnet_angle), NULL, 0.0, 0},
    };

TWMTA twmta_example;
/* names for traveling-wave muffin-tin linac parameters
 */
PARAMETER twmta_param[N_TWMTA_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&twmta_example.length), NULL, 0.0, 0},
    {"FREQUENCY", "HZ", IS_DOUBLE, 0, (long)((char *)&twmta_example.frequency), NULL, DEFAULT_FREQUENCY, 0},
    {"PHASE", "RAD", IS_DOUBLE, 0, (long)((char *)&twmta_example.phase), NULL, 0.0, 0},
    {"EZ", "V/M", IS_DOUBLE, 0, (long)((char *)&twmta_example.Ez), NULL, 0.0, 0},
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&twmta_example.accuracy), NULL, DEFAULT_ACCURACY, 0},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twmta_example.x_max), NULL, 0.0, 0},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twmta_example.y_max), NULL, 0.0, 0},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&twmta_example.dx), NULL, 0.0, 0},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&twmta_example.dy), NULL, 0.0, 0},
    {"KX", "1/M", IS_DOUBLE, 0, (long)((char *)&twmta_example.kx), NULL, 0.0, 0},
    {"BETA_WAVE", "", IS_DOUBLE, 0, (long)((char *)&twmta_example.beta_wave), NULL, DEFAULT_BETA_WAVE, 0},
    {"BSOL", "", IS_DOUBLE, 0, (long)((char *)&twmta_example.Bsol), NULL, 0.0, 0},
    {"ALPHA", "1/M", IS_DOUBLE, 0, (long)((char *)&twmta_example.alpha), NULL, 0.0, 0},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&twmta_example.phase_reference), NULL, 0.0, 0},
    {"N_STEPS", "", IS_LONG, 0, (long)((char *)&twmta_example.n_steps), NULL, 0.0, 100},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&twmta_example.method), DEFAULT_INTEG_METHOD, 0.0, 0},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&twmta_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0}
    } ;

MATTER matter_example;
/* matter physical parameters */
PARAMETER matter_param[N_MATTER_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&matter_example.length), NULL, 0.0, 0},
    {"XO", "M", IS_DOUBLE, 1, (long)((char *)&matter_example.Xo), NULL, 0.0, 0},
    {"ELASTIC", "", IS_LONG, 0, (long)((char *)&matter_example.elastic), NULL, 0.0, 0},
    };

RFMODE rfmode_example;
/* RFMODE physical parameters */
PARAMETER rfmode_param[N_RFMODE_PARAMS] = {
    {"RA", "Ohm", IS_DOUBLE, 0, (long)((char *)&rfmode_example.Ra), NULL, 0.0, 0},
    {"RS", "Ohm", IS_DOUBLE, 0, (long)((char *)&rfmode_example.Rs), NULL, 0.0, 0},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&rfmode_example.Q), NULL, 0.0, 1},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&rfmode_example.freq), NULL, 0.0, 0},
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&rfmode_example.charge), NULL, 0.0, 0},
    {"INITIAL_V", "V", IS_DOUBLE, 0, (long)((char *)&rfmode_example.initial_V), NULL, 0.0, 0},
    {"INITIAL_PHASE", "DEGREES", IS_DOUBLE, 0, (long)((char *)&rfmode_example.initial_phase), NULL, 0.0, 0},
    {"BETA", "", IS_DOUBLE, 0, (long)((char *)&rfmode_example.beta), NULL, 0.0, 0},
    {"BIN_SIZE", "S", IS_DOUBLE, 0, (long)((char *)&rfmode_example.bin_size), NULL, 0.0, 0},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&rfmode_example.n_bins), NULL, 0.0, 20},
    {"PRELOAD", "", IS_LONG, 0, (long)((char *)&rfmode_example.preload), NULL, 0.0, 0},
    {"PRELOAD_FACTOR", "", IS_DOUBLE, 0, (long)((char *)&rfmode_example.preload_factor), NULL, 1.0, 0},
    {"RIGID_UNTIL_PASS", "", IS_LONG, 0, (long)((char *)&rfmode_example.rigid_until_pass), NULL, 0.0, 0},
    {"SAMPLE_INTERVAL", "", IS_LONG, 0, (long)((char *)&rfmode_example.sample_interval), NULL, 0.0, 1},
    {"RECORD", "", IS_STRING, 0, (long)((char *)&rfmode_example.record), NULL, 0.0, 0},
    {"SINGLE_PASS", "", IS_LONG, 0, (long)((char *)&rfmode_example.single_pass), NULL, 0.0, 0},
    {"PASS_INTERVAL", "", IS_LONG, 0, (long)((char *)&rfmode_example.pass_interval), NULL, 0.0, 1},
    };

TRFMODE trfmode_example;
/* TRFMODE physical parameters */
PARAMETER trfmode_param[N_TRFMODE_PARAMS] = {
    {"RA", "Ohm", IS_DOUBLE, 0, (long)((char *)&trfmode_example.Ra), NULL, 0.0, 0},
    {"RS", "Ohm", IS_DOUBLE, 0, (long)((char *)&trfmode_example.Rs), NULL, 0.0, 0},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&trfmode_example.Q), NULL, 0.0, 1},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&trfmode_example.freq), NULL, 0.0, 0},
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&trfmode_example.charge), NULL, 0.0, 0},
    {"BETA", "", IS_DOUBLE, 0, (long)((char *)&trfmode_example.beta), NULL, 0.0, 0},
    {"BIN_SIZE", "S", IS_DOUBLE, 0, (long)((char *)&trfmode_example.bin_size), NULL, 0.0, 0},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&trfmode_example.n_bins), NULL, 0.0, 20},
    {"PLANE", "", IS_STRING, 0, (long)((char *)&trfmode_example.plane), "both", 0.0, 0},
    {"SINGLE_PASS", "", IS_LONG, 0, (long)((char *)&trfmode_example.single_pass), NULL, 0.0, 0},
    };

ZLONGIT zlongit_example;
/* ZLONGIT physical parameters */
PARAMETER zlongit_param[N_ZLONGIT_PARAMS] = {
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&zlongit_example.charge), NULL, 0.0, 0},
    {"BROAD_BAND", "", IS_LONG, 0, (long)((char *)&zlongit_example.broad_band), NULL, 0.0, 0},
    {"RA", "Ohm", IS_DOUBLE, 0, (long)((char *)&zlongit_example.Ra), NULL, 0.0, 0},
    {"RS", "Ohm", IS_DOUBLE, 0, (long)((char *)&zlongit_example.Rs), NULL, 0.0, 0},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&zlongit_example.Q), NULL, 0.0, 1},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&zlongit_example.freq), NULL, 0.0, 0},
    {"ZREAL", "", IS_STRING, 0, (long)((char *)&zlongit_example.Zreal), NULL, 0.0, 0},
    {"ZIMAG", "", IS_STRING, 0, (long)((char *)&zlongit_example.Zimag), NULL, 0.0, 0},
    {"BIN_SIZE", "S", IS_DOUBLE, 0, (long)((char *)&zlongit_example.bin_size), NULL, 0.0, 0},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&zlongit_example.n_bins), NULL, 0.0, 128},
    {"WAKES", "", IS_STRING, 0, (long)((char *)&zlongit_example.wakes), NULL, 0.0, 0},
    {"WAKE_INTERVAL", "", IS_LONG, 0, (long)((char *)&zlongit_example.wake_interval), NULL, 0.0, 1},
    {"AREA_WEIGHT", "", IS_LONG, 0, (long)((char *)&zlongit_example.area_weight), NULL, 0.0, 0},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&zlongit_example.interpolate), NULL, 0.0, 0},
    {"SMOOTH_PASSES", "", IS_LONG, 0, (long)((char *)&zlongit_example.smooth_passes), NULL, 0.0, 0},
    };

SREFFECTS SReffects_example;
/* SREFFECTS physical parameters */
PARAMETER sreffects_param[N_SREFFECTS_PARAMS] = {
    {"JX", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.Jx), NULL, 1.0, 0},
    {"JY", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.Jy), NULL, 1.0, 0},
    {"JDELTA", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.Jdelta), NULL, 2.0, 0},
    {"EXREF", "m", IS_DOUBLE, 0, (long)((char *)&SReffects_example.exRef), NULL, 0.0, 0},
    {"EYREF", "m", IS_DOUBLE, 0, (long)((char *)&SReffects_example.eyRef), NULL, 0.0, 0},
    {"SDELTAREF", "m", IS_DOUBLE, 0, (long)((char *)&SReffects_example.SdeltaRef), NULL, 0.0, 0},
    {"DDELTAREF", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.DdeltaRef), NULL, 0.0, 0},
    {"PREF", "m$be$nc", IS_DOUBLE, 0, (long)((char *)&SReffects_example.pRef), NULL, 0.0, 0},
    {"COUPLING", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.coupling), NULL, 0.0, 0},
    {"FRACTION", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.fraction), NULL, 1.0, 0},
    };

MODRF modrf_example;
PARAMETER modrf_param[N_MODRF_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&modrf_example.length), NULL, 0.0, 0},
    {"VOLT", "V", IS_DOUBLE, 0, (long)((char *)&modrf_example.volt), NULL, 0.0, 0},
    {"PHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&modrf_example.phase), NULL, 0.0, 0},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&modrf_example.freq), NULL, 500.0e6, 0},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&modrf_example.Q), NULL, 0.0, 0},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&modrf_example.phase_reference), NULL, 0.0, 0},
    {"AMMAG", "", IS_DOUBLE, 0, (long)((char *)&modrf_example.amMag), NULL, 0.0, 0},
    {"AMPHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&modrf_example.amPhase), NULL, 0.0, 0},
    {"AMFREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&modrf_example.amFreq), NULL, 0.0, 0},
    {"PMMAG", "DEG", IS_DOUBLE, 0, (long)((char *)&modrf_example.pmMag), NULL, 0.0, 0},
    {"PMPHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&modrf_example.pmPhase), NULL, 0.0, 0},
    {"PMFREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&modrf_example.pmFreq), NULL, 0.0, 0},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&modrf_example.fiducial), NULL, 0.0, 0},
    };    

BMAPXY bmapxy_example;
PARAMETER bmapxy_param[N_BMAPXY_PARAMS] = {
{"L", "M", IS_DOUBLE, 0, (long)((char *)&bmapxy_example.length), NULL, 0.0, 0},
{"STRENGTH", NULL, IS_DOUBLE, 0, (long)((char *)&bmapxy_example.strength), NULL, 0.0, 0},
{"ACCURACY", NULL, IS_DOUBLE, 0, (long)((char *)&bmapxy_example.accuracy), NULL, 0.0, 0},
{"METHOD", NULL, IS_STRING, 0, (long)((char*)&bmapxy_example.method), NULL, 0.0, 0},
{"FILENAME", NULL, IS_STRING, 0, (long)((char*)&bmapxy_example.filename), NULL, 0.0, 0},
};  

ZTRANSVERSE ztransverse_example;
PARAMETER ztransverse_param[N_ZTRANSVERSE_PARAMS] = {
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.charge), NULL, 0.0, 0},
    {"BROAD_BAND", "", IS_LONG, 0, (long)((char *)&ztransverse_example.broad_band), NULL, 0.0, 0},
    {"RS", "Ohm", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.Rs), NULL, 0.0, 0},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.Q), NULL, 0.0, 1},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.freq), NULL, 0.0, 0},
    {"INPUTFILE", "", IS_STRING, 0, (long)((char*)&ztransverse_example.inputFile), NULL, 0.0, 0},
    {"FREQCOLUMN", "", IS_STRING, 0, (long)((char*)&ztransverse_example.freqColumn), NULL, 0.0, 0},
    {"ZXREAL", "", IS_STRING, 0, (long)((char *)&ztransverse_example.ZxReal), NULL, 0.0, 0},
    {"ZXIMAG", "", IS_STRING, 0, (long)((char *)&ztransverse_example.ZxImag), NULL, 0.0, 0},
    {"ZYREAL", "", IS_STRING, 0, (long)((char *)&ztransverse_example.ZyReal), NULL, 0.0, 0},
    {"ZYIMAG", "", IS_STRING, 0, (long)((char *)&ztransverse_example.ZyImag), NULL, 0.0, 0},
    {"BIN_SIZE", "S", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.bin_size), NULL, 0.0, 0},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&ztransverse_example.interpolate), NULL, 0.0, 0},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&ztransverse_example.n_bins), NULL, 0.0, 128},
};

IBSCATTER ibs_example;
PARAMETER ibscatter_param[N_IBSCATTER_PARAMS] = {
  {"COUPLING", "", IS_DOUBLE, 0, (long)((char *)&ibs_example.coupling), NULL, 1.0, 0},
  {"FRACTION", "", IS_DOUBLE, 0, (long)((char *)&ibs_example.fraction), NULL, 1.0, 0},
  {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&ibs_example.charge), NULL, 0.0, 0},
};

WAKE wake_example;
/* WAKE physical parameters */
PARAMETER wake_param[N_WAKE_PARAMS] = {
    {"INPUTFILE", "", IS_STRING, 0, (long)((char *)&wake_example.inputFile), NULL, 0.0, 0},
    {"TCOLUMN", "", IS_STRING, 0, (long)((char *)&wake_example.tColumn), NULL, 0.0, 0},
    {"WCOLUMN", "", IS_STRING, 0, (long)((char *)&wake_example.WColumn), NULL, 0.0, 0},
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&wake_example.charge), NULL, 0.0, 0},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&wake_example.n_bins), NULL, 0.0, 128},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&wake_example.interpolate), NULL, 0.0, 0},
    {"SMOOTH_PASSES", "", IS_LONG, 0, (long)((char *)&wake_example.smooth_passes), NULL, 0.0, 0},
    };

TRWAKE trwake_example;
/* TRWAKE physical parameters */
PARAMETER trwake_param[N_TRWAKE_PARAMS] = {
    {"INPUTFILE", "", IS_STRING, 0, (long)((char *)&trwake_example.inputFile), NULL, 0.0, 0},
    {"TCOLUMN", "", IS_STRING, 0, (long)((char *)&trwake_example.tColumn), NULL, 0.0, 0},
    {"WXCOLUMN", "", IS_STRING, 0, (long)((char *)&trwake_example.WxColumn), NULL, 0.0, 0},
    {"WYCOLUMN", "", IS_STRING, 0, (long)((char *)&trwake_example.WyColumn), NULL, 0.0, 0},
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&trwake_example.charge), NULL, 0.0, 0},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&trwake_example.n_bins), NULL, 0.0, 128},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&trwake_example.interpolate), NULL, 0.0, 0},
    {"SMOOTH_PASSES", "", IS_LONG, 0, (long)((char *)&trwake_example.smooth_passes), NULL, 0.0, 0},
    };

/* array of parameter structures */

#define MAT_LEN     HAS_MATRIX|HAS_LENGTH
#define MAT_LEN_NCAT HAS_MATRIX|HAS_LENGTH|DONT_CONCAT

ELEMENT_DESCRIPTION entity_description[N_TYPES] = {
    {                0,           0,                  0,    NULL           },
    {    N_QUAD_PARAMS,     MAT_LEN|IS_MAGNET|MATRIX_TRACKING,       sizeof(QUAD),    quad_param     },
    {    N_BEND_PARAMS,     MAT_LEN|IS_MAGNET|MATRIX_TRACKING,       sizeof(BEND),    bend_param     },
    {    N_BEND_PARAMS,     MAT_LEN|IS_MAGNET|MATRIX_TRACKING,       sizeof(BEND),    bend_param     }, 
    {   N_DRIFT_PARAMS,     MAT_LEN|MATRIX_TRACKING,      sizeof(DRIFT),    drift_param    }, 
    {    N_SEXT_PARAMS,     MAT_LEN|IS_MAGNET|MATRIX_TRACKING,       sizeof(SEXT),    sext_param     },
    {                0,           0,                  0,    NULL           },
    {    N_MULT_PARAMS,  MAT_LEN_NCAT|IS_MAGNET,       sizeof(MULT),    mult_param     }, 
    {    N_SOLE_PARAMS,     MAT_LEN|IS_MAGNET|MATRIX_TRACKING,       sizeof(SOLE),    sole_param     }, 
    {    N_HCOR_PARAMS,     MAT_LEN|IS_MAGNET|MATRIX_TRACKING,       sizeof(HCOR),    hcor_param     }, 
    {    N_VCOR_PARAMS,     MAT_LEN|IS_MAGNET|MATRIX_TRACKING,       sizeof(VCOR),    vcor_param     }, 
    {    N_RFCA_PARAMS,     MAT_LEN_NCAT,       sizeof(RFCA),    rfca_param     }, 
    {                0,           0,                  0,    NULL           },
    {    N_HMON_PARAMS,     MAT_LEN_NCAT|MATRIX_TRACKING,       sizeof(HMON),    hmon_param     }, 
    {    N_VMON_PARAMS,     MAT_LEN_NCAT|MATRIX_TRACKING,       sizeof(VMON),    vmon_param     }, 
    {    N_MONI_PARAMS,     MAT_LEN_NCAT|MATRIX_TRACKING,       sizeof(MONI),    moni_param     }, 
    {    N_RCOL_PARAMS,  MAT_LEN_NCAT,       sizeof(RCOL),    rcol_param     }, 
    {    N_ECOL_PARAMS,  MAT_LEN_NCAT,       sizeof(ECOL),    ecol_param     }, 
    {    N_MARK_PARAMS,           0,       sizeof(MARK),    mark_param     }, 
    {    N_MATR_PARAMS,  MAT_LEN|MATRIX_TRACKING,       sizeof(MATR),    matr_param     }, 
    {    N_ALPH_PARAMS,  HAS_MATRIX|IS_MAGNET,       sizeof(ALPH),    alph_param     }, 
    {    N_RFDF_PARAMS,  MAT_LEN_NCAT,       sizeof(RFDF),    rfdf_param     }, 
    {    N_RFTM_PARAMS,  MAT_LEN_NCAT,    sizeof(TM_MODE),    rftm_param     }, 
    {    N_RMDF_PARAMS,  MAT_LEN_NCAT,       sizeof(RMDF),    rmdf_param     }, 
    {    N_TMCF_PARAMS,  MAT_LEN_NCAT,  sizeof(TMCF_MODE),    tmcf_param     }, 
    {    N_CEPL_PARAMS,  MAT_LEN_NCAT,  sizeof(CE_PLATES),    cepl_param     }, 
    {   N_WATCH_PARAMS,           0,      sizeof(WATCH),    watch_param    }, 
    {    N_TWPL_PARAMS,  MAT_LEN_NCAT,  sizeof(TW_PLATES),    twpl_param     }, 
    {  N_MALIGN_PARAMS,  HAS_MATRIX|DONT_CONCAT,
                                         sizeof(MALIGN),    malign_param   },
    {    N_TWLA_PARAMS,  MAT_LEN_NCAT,   sizeof(TW_LINAC),    twla_param     },
    {  N_PEPPOT_PARAMS,  MAT_LEN_NCAT,     sizeof(PEPPOT),    peppot_param   },
    {  N_ENERGY_PARAMS,           0,     sizeof(ENERGY),    energy_param   },
    {  N_MAXAMP_PARAMS,           0,     sizeof(MAXAMP),    maxamp_param   },
    {  N_ROTATE_PARAMS,  HAS_MATRIX|MATRIX_TRACKING,     sizeof(ROTATE),    rotate_param   },
    { N_TRCOUNT_PARAMS,           0,    sizeof(TRCOUNT),    trcount_param  },
    {  N_RECIRC_PARAMS,           0,     sizeof(RECIRC),    recirc_param   },
    {  N_QFRING_PARAMS,     MAT_LEN|MATRIX_TRACKING,     sizeof(QFRING),    qfring_param   },
    { N_SCRAPER_PARAMS,  MAT_LEN_NCAT,    sizeof(SCRAPER),    scraper_param  },
    {  N_CENTER_PARAMS,           0,     sizeof(CENTER),    center_param   },
    {  N_KICKER_PARAMS,  MAT_LEN_NCAT|IS_MAGNET,     sizeof(KICKER),    kicker_param   },
    {   N_KSEXT_PARAMS, MAT_LEN_NCAT|IS_MAGNET,      
                                          sizeof(KSEXT),    ksext_param    },
    {  N_KSBEND_PARAMS, MAT_LEN_NCAT|IS_MAGNET,
                                         sizeof(KSBEND),    ksbend_param   },
    {   N_KQUAD_PARAMS, MAT_LEN_NCAT|IS_MAGNET, 
                                          sizeof(KQUAD),    kquad_param    },
    { N_MAGNIFY_PARAMS, HAS_MATRIX|MATRIX_TRACKING,     sizeof(MAGNIFY),    magnify_param  },
    {  N_SAMPLE_PARAMS,          0,      sizeof(SAMPLE),    sample_param   },
    {   N_HVCOR_PARAMS,    MAT_LEN|IS_MAGNET|MATRIX_TRACKING,       sizeof(HVCOR),    hvcor_param    }, 
    { N_SCATTER_PARAMS,          0,     sizeof(SCATTER),    scatter_param  },
    {  N_NIBEND_PARAMS, MAT_LEN_NCAT|IS_MAGNET,
                                         sizeof(NIBEND),    nibend_param   },
    {   N_KPOLY_PARAMS,          0,       sizeof(KPOLY),    kpoly_param    }, 
    {  N_NISEPT_PARAMS, MAT_LEN_NCAT|IS_MAGNET,
                                         sizeof(NISEPT),    nisept_param   },
    {  N_RAMPRF_PARAMS, MAT_LEN_NCAT,    sizeof(RAMPRF),    ramprf_param   },
    {   N_RAMPP_PARAMS,          0,       sizeof(RAMPP),    rampp_param    },
    {   N_STRAY_PARAMS,    MAT_LEN|MATRIX_TRACKING,       sizeof(STRAY),    stray_param    },
    {  N_CSBEND_PARAMS, MAT_LEN_NCAT|IS_MAGNET,
                                         sizeof(CSBEND),    csbend_param   },
    {   N_TWMTA_PARAMS, MAT_LEN_NCAT,     sizeof(TWMTA),    twmta_param    },
    {  N_MATTER_PARAMS,    MAT_LEN,      sizeof(MATTER),   matter_param    },
    {  N_RFMODE_PARAMS,          0,      sizeof(RFMODE),   rfmode_param    },
    { N_TRFMODE_PARAMS,          0,     sizeof(TRFMODE),  trfmode_param    },
    { N_ZLONGIT_PARAMS,          0,     sizeof(ZLONGIT),  zlongit_param    },
    { N_SREFFECTS_PARAMS,        0,   sizeof(SREFFECTS),  sreffects_param  },
    { N_MODRF_PARAMS, MAT_LEN_NCAT,       sizeof(MODRF),    modrf_param     }, 
    { N_BMAPXY_PARAMS,     HAS_LENGTH,   sizeof(BMAPXY),  bmapxy_param      },
    { N_ZTRANSVERSE_PARAMS,          0,     sizeof(ZTRANSVERSE),  ztransverse_param    },
    { N_IBSCATTER_PARAMS,        0,   sizeof(IBSCATTER),  ibscatter_param  },
    { N_FMULT_PARAMS,  MAT_LEN_NCAT|IS_MAGNET,       sizeof(FMULT),    fmult_param     }, 
    { N_WAKE_PARAMS, 0, sizeof(WAKE), wake_param},
    { N_TRWAKE_PARAMS, 0, sizeof(TRWAKE), trwake_param},
    { N_TUBEND_PARAMS, 0, sizeof(TUBEND), tubend_param},
} ;
 

void compute_offsets()
{
    long i, j;
    for (i=0; i<N_TYPES; i++) {
        for (j=entity_description[i].n_params-1; j>=0; j--)
            entity_description[i].parameter[j].offset -= entity_description[i].parameter[0].offset;
        }
    }

