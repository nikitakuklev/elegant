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
    "MARK", "MATR", "ALPH", "RFDF", "RFTMEZ0", "RMDF", "TMCF", "CEPL", "WATCH",
    "TWPL", "MALIGN", "TWLA", "PEPPOT", "ENERGY", "MAXAMP", "ROTATE",
    "TRCOUNT", "RECIRC", "QUFRINGE", "SCRAPER", "CENTER", "BUMPER",
    "KSEXT", "KSBEND", "KQUAD", "MAGNIFY", "SAMPLE", "KICKER", "SCATTER",
    "NIBEND", "KPOLY", "NISEPT", "RAMPRF", "RAMPP", "STRAY", "CSBEND",
    "TWMTA", "MATTER", "RFMODE", "TRFMODE", "ZLONGIT", "SREFFECTS",
    "MODRF", "BMAPXY", "ZTRANSVERSE", "IBSCATTER", "FMULT",
    "WAKE", "TRWAKE", "TUBEND", "CHARGE", "PFILTER", "HISTOGRAM",
    "CSRCSBEND", "CSRDRIFT", "RFCW", "REMCOR", "MAPSOLENOID",
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
the magnet into halves.  XS<n> and DP<n> allow momentum filtration at the midpoint.",
    "A simple traveling-wave (beta=1) deflecting RF cavity.",
    "Not implemented.",
    "A linearly-ramped electric field deflector, using an approximate analytical solution FOR LOW ENERGY PARTICLES.",
    "A numerically-integrated accelerating TM RF cavity with spatially-constant fields.",
    "A numerically-integrated linearly-ramped electric field deflector.",
    "A beam property/motion monitor--allowed modes are centroid, parameter, coordinate, and fft.",
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
matrix implementation.  Recommend using CSBEND instead.",
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
defined relative to the initial beamline direction.  ** Not correct if there are tilts\n\
in the beamline. **",
    "A canonical kick sector dipole magnet.",
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
    "An element to establish the total charge of a beam.  Active on first pass only.  If given, overrides all charge specifications on other elements.",
    "An element for energy and momentum filtration.",
    "Request for histograms of particle coordinates to be output to SDDS file.",
    "Like CSBEND, but incorporates a simulation of Coherent Synchrotron radiation.",
    "A follow-on element for CSRCSBEND that applies the CSR wake over a drift.",
    "A combination of RFCA, WAKE, and TRWAKE.",
    "An element to remove correlations from the tracked beam to simulate certain types of correction.",
    "A numerically-integrated solenoid specified as a map of (Bz, Br) vs (z, r)."
    } ;

QUAD quad_example;
/* quadrupole physical parameters */
PARAMETER quad_param[N_QUAD_PARAMS]={
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&quad_example.length), NULL, 0.0, 0, "length"},
    {"K1", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&quad_example.k1), NULL, 0.0, 0, "geometric strength"},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&quad_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"FFRINGE", "", IS_DOUBLE, 1, (long)((char *)&quad_example.ffringe), NULL, 0.0, 0, "fraction of length occupied by linear fringe region"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&quad_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&quad_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&quad_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "M", IS_DOUBLE, 1, (long)((char *)&quad_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&quad_example.order), NULL, 0.0, 0, "matrix order"}
    };

BEND bend_example;
/* bending magnet physical parameters */
PARAMETER bend_param[N_BEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&bend_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, 1, (long)((char *)&bend_example.angle), NULL, 0.0, 0, "bend angle"},
    {"K1", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&bend_example.k1), NULL, 0.0, 0, "geometric focusing strength"},
    {"E1", "RAD", IS_DOUBLE, 1, (long)((char *)&bend_example.e1), NULL, 0.0, 0, "entrance edge angle"},
    {"E2", "RAD", IS_DOUBLE, 1, (long)((char *)&bend_example.e2), NULL, 0.0, 0, "exit edge angle"},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&bend_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"K2", "1/M$a3$n", IS_DOUBLE, 1, (long)((char *)&bend_example.k2), NULL, 0.0, 0, "geometric sextupole strength"},
    {"H1", "1/M", IS_DOUBLE, 1, (long)((char *)&bend_example.h1), NULL, 0.0, 0, "entrace pole-face curvature"},
    {"H2", "1/M", IS_DOUBLE, 1, (long)((char *)&bend_example.h2), NULL, 0.0, 0, "exit pole-face curvature"},
    {"HGAP", "M", IS_DOUBLE, 1, (long)((char *)&bend_example.hgap), NULL, 0.0, 0, "half-gap between poles"},
    {"FINT", "", IS_DOUBLE, 1, (long)((char *)&bend_example.fint), NULL, DEFAULT_FINT, 0, "edge-field integral"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&bend_example.dx), NULL, 0.0, 0, "misaligment of entrance"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&bend_example.dy), NULL, 0.0, 0, "misalignment of entrace"},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&bend_example.dz), NULL, 0.0, 0, "misalignment of entrance"},
    {"FSE", "", IS_DOUBLE, 1, (long)((char *)&bend_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"ETILT", "RAD", IS_DOUBLE, 1, (long)((char *)&bend_example.etilt), NULL, 0.0, 0, "error rotation about incoming longitudinal axis"},
    {"EDGE1_EFFECTS", "", IS_LONG, 1, (long)((char *)&bend_example.edge1_effects), NULL, 0.0, 1, "include entrace edge effects?"},
    {"EDGE2_EFFECTS", "", IS_LONG, 1, (long)((char *)&bend_example.edge2_effects), NULL, 0.0, 1, "include exit edge effects?"},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&bend_example.order), NULL, 0.0, 0, "matrix order"},
    {"EDGE_ORDER", "", IS_LONG, 1, (long)((char *)&bend_example.edge_order), NULL, 0.0, 0, "edge matrix order"},
    {"TRANSPORT", "", IS_LONG, 1, (long)((char *)&bend_example.TRANSPORT), NULL, 0.0, 0, "use (incorrect) TRANSPORT equations for T436 of edge?"},
    {"USE_BN", "", IS_LONG, 0, (long)((char *)&bend_example.use_bn), NULL, 0.0, 0, "use B1 and B2 instead of K1 and K2 values?"},
    {"B1", "1/M", IS_DOUBLE, 0, (long)((char *)&bend_example.b1), NULL, 0.0, 0, "K1 = B1*rho, where rho is bend radius"},
    {"B2", "1/M$a2$n", IS_DOUBLE, 0, (long)((char *)&bend_example.b2), NULL, 0.0, 0, "K2 = B2*rho"},
    };

DRIFT drift_example;
/* drift length physical parameters */
PARAMETER drift_param[N_DRIFT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&drift_example.length), NULL, 0.0, 0, "length"},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&drift_example.order), NULL, 0.0, 0, "matrix order"}
    };

SEXT sext_example;
/* sextupole physical parameters */
PARAMETER sext_param[N_SEXT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&sext_example.length), NULL, 0.0, 0, "length"},
    {"K2", "1/M$a3$n", IS_DOUBLE, 1, (long)((char *)&sext_example.k2), NULL, 0.0, 0, "geometric strength"},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&sext_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&sext_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&sext_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&sext_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "M", IS_DOUBLE, 1, (long)((char *)&sext_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&sext_example.order), NULL, 0.0, 0, "matrix order"}
    };
   
MULT mult_example;
/* multipole physical parameters */
PARAMETER mult_param[N_MULT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.length), NULL, 0.0, 0, "length"},
    {"KNL", "M$a(1-ORDER)$n", IS_DOUBLE, 0, (long)((char *)&mult_example.KnL), NULL, 0.0, 0, "integrated geometric strength"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&mult_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"BORE", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.bore), NULL, 0.0, 0, "bore radius"},
    {"BNL", "T M", IS_DOUBLE, 0, (long)((char *)&mult_example.BnL), NULL, 0.0, 0, "integrated field at pole tip, used if BORE nonzero"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 0, (long)((char *)&mult_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FACTOR", "", IS_DOUBLE, 0, (long)((char *)&mult_example.factor), NULL, 1.0, 0, "factor by which to multiply strength"},
    {"ORDER", "", IS_LONG, 0, (long)((char *)&mult_example.order), NULL, 0.0, 1, "multipole order"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&mult_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks"},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&mult_example.synch_rad), NULL, 0.0, 0, "include classical synchrotron radiation?"},
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
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&fmult_example.synch_rad), NULL, 0.0, 0, "include classical synchrotron radiation?"},
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&fmult_example.filename), NULL, 0.0, 0, "name of file containing multipole data"},
    };

SOLE sole_example;
/* solenoid physical parameters */
PARAMETER sole_param[N_SOLE_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&sole_example.length), NULL, 0.0, 0, "length"},
    {"KS", "RAD/M", IS_DOUBLE, 1, (long)((char *)&sole_example.ks), NULL, 0.0, 0, "geometric strength"},
    {"B", "T", IS_DOUBLE, 1, (long)((char *)&sole_example.B), NULL, 0.0, 0, "field strength (used if KS is zero)"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&sole_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&sole_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&sole_example.dz), NULL, 0.0, 0, "misalignment"},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&sole_example.order), NULL, 0.0, 0, "matrix order"},
    };
   
HCOR hcor_example;
/* horizontal corrector physical parameters */
PARAMETER hcor_param[N_HCOR_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&hcor_example.length), NULL, 0.0, 0, "length"},
    {"KICK", "RAD", IS_DOUBLE, 1, (long)((char *)&hcor_example.kick), NULL, 0.0, 0, "kick strength"},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&hcor_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"B2", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&hcor_example.b2), NULL, 0.0, 0, "normalized sextupole strength (kick = KICK*(1+B2*x^2))"},
    {"CALIBRATION", "", IS_DOUBLE, 1, (long)((char *)&hcor_example.calibration), NULL, 1.0, 0, "strength multiplier"},
    {"EDGE_EFFECTS", "", IS_LONG, 0, (long)((char *)&hcor_example.edge_effects), NULL, 0.0, 0, "include edge effects?"},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&hcor_example.order), NULL, 0.0, 0, "matrix order"},
    {"STEERING", "", IS_LONG, 0, (long)((char *)&hcor_example.steering), NULL, 0.0, 1, "use for steering?"},
    };

VCOR vcor_example;
/* vertical corrector physical parameters */
PARAMETER vcor_param[N_VCOR_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&vcor_example.length), NULL, 0.0, 0, "length"},
    {"KICK", "RAD", IS_DOUBLE, 1, (long)((char *)&vcor_example.kick), NULL, 0.0, 0, "kick strength"},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&vcor_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"B2", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&vcor_example.b2), NULL, 0.0, 0, "normalized sextupole strength (kick = KICK*(1+B2*y^2))"},
    {"CALIBRATION", "", IS_DOUBLE, 1, (long)((char *)&vcor_example.calibration), NULL, 1.0, 0, "strength multiplier"},
    {"EDGE_EFFECTS", "", IS_LONG, 0, (long)((char *)&vcor_example.edge_effects), NULL, 0.0, 0, "include edge effects?"},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&vcor_example.order), NULL, 0.0, 0, "matrix order"},
    {"STEERING", "", IS_LONG, 0, (long)((char *)&vcor_example.steering), NULL, 0.0, 1, "use for steering?"},
    };

RFCA rfca_example;
/* rf cavity physical parameters */
PARAMETER rfca_param[N_RFCA_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&rfca_example.length), NULL, 0.0, 0, "length"},
    {"VOLT", "V", IS_DOUBLE, 0, (long)((char *)&rfca_example.volt), NULL, 0.0, 0, "peak voltage"},
    {"PHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&rfca_example.phase), NULL, 0.0, 0, "phase"},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&rfca_example.freq), NULL, 500.0e6, 0, "frequency"},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&rfca_example.Q), NULL, 0.0, 0, "cavity Q"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&rfca_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"CHANGE_P0", "", IS_LONG, 0, (long)((char *)&rfca_example.change_p0), NULL, 0.0, 0, "does cavity change central momentum?"}, 
    {"CHANGE_T", "", IS_LONG, 0, (long)((char *)&rfca_example.change_t), NULL, 0.0, 0, "not recommended"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&rfca_example.fiducial), NULL, 0.0, 0, "mode for determining fiducial arrival time (light, tmean, first, pmaximum)"},
    {"END1_FOCUS", "", IS_LONG, 0, (long)((char *)&rfca_example.end1Focus), NULL, 0.0, 0, "include focusing at entrance?"},
    {"END2_FOCUS", "", IS_LONG, 0, (long)((char *)&rfca_example.end2Focus), NULL, 0.0, 0, "include focusing at exit?"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&rfca_example.nKicks), NULL, 0.0, 1, "number of kicks to use.  Set to zero for matrix method."},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&rfca_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&rfca_example.dy), NULL, 0.0, 0, "misalignment"},
    {"T_REFERENCE", "S", IS_DOUBLE, 1, (long)((char *)&rfca_example.tReference), NULL, -1.0, 0, "arrival time of reference particle"},
    };
   
HMON hmon_example;
/* name for horizontal monitor physical parameters */
PARAMETER hmon_param[N_HMON_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&hmon_example.length), NULL, 0.0, 0, "length"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&hmon_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&hmon_example.dy), NULL, 0.0, 0, "misalignment"},
    {"WEIGHT", "", IS_DOUBLE, 0, (long)((char *)&hmon_example.weight), NULL, 1.0, 0, "weight in correction"},
    {"TILT", "", IS_DOUBLE, 0, (long)((char *)&hmon_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"CALIBRATION", "", IS_DOUBLE, 0, (long)((char *)&hmon_example.calibration), NULL, 1.0, 0, "calibration factor for readout"},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&hmon_example.order), NULL, 1.0, 0, "matrix order"},
    {"READOUT", "", IS_STRING, 0, (long)((char *)&hmon_example.readout), NULL, 0.0, 1, "rpn expression for readout (actual position supplied in variable x)"}
    } ;
   
VMON vmon_example;
/* name for vertical monitor physical parameters */
PARAMETER vmon_param[N_VMON_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&vmon_example.length), NULL, 0.0, 0, "length"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&vmon_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&vmon_example.dy), NULL, 0.0, 0, "misalignment"},
    {"WEIGHT", "", IS_DOUBLE, 0, (long)((char *)&vmon_example.weight), NULL, 1.0, 0, "weight in correction"},
    {"TILT", "", IS_DOUBLE, 0, (long)((char *)&vmon_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"CALIBRATION", "", IS_DOUBLE, 0, (long)((char *)&vmon_example.calibration), NULL, 1.0, 0, "calibration factor for readout"},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&vmon_example.order), NULL, 1.0, 0, "matrix order"},
    {"READOUT", "", IS_STRING, 0, (long)((char *)&vmon_example.readout), NULL, 0.0, 1, "rpn expression for readout (actual position supplied in variable y)"}
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
    {"ORDER", "", IS_LONG, 1, (long)((char *)&moni_example.order), NULL, 1.0, 0, "matrix order"},
    {"XREADOUT", "", IS_STRING, 0, (long)((char *)&moni_example.x_readout), NULL, 0.0, 1, "rpn expression for x readout (actual position supplied in variables x, y"},
    {"YREADOUT", "", IS_STRING, 0, (long)((char *)&moni_example.y_readout), NULL, 0.0, 1, "rpn expression for y readout (actual position supplied in variables x, y"}
    } ;

RCOL rcol_example;   
/* name for rectangular collimator physical parameters */
PARAMETER rcol_param[N_RCOL_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&rcol_example.length), NULL, 0.0, 0, "length"},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&rcol_example.x_max), NULL, 0.0, 0, "half-width in x"},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&rcol_example.y_max), NULL, 0.0, 0, "half-width in y"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&rcol_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&rcol_example.dy), NULL, 0.0, 0, "misalignment"},
    } ;
   
ECOL ecol_example;
/* name for elliptical collimator physical parameters */
PARAMETER ecol_param[N_ECOL_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&ecol_example.length), NULL, 0.0, 0, "length"},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&ecol_example.x_max), NULL, 0.0, 0, "half-axis in x"},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&ecol_example.y_max), NULL, 0.0, 0, "half-axis in y"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&ecol_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&ecol_example.dy), NULL, 0.0, 0, "misalignment"},
    } ;

MARK mark_example;
/* name for marker parameters */
PARAMETER mark_param[N_MARK_PARAMS] = {
    {"FITPOINT", "", IS_LONG, 0, (long)((char *)&mark_example.fitpoint), NULL, 0.0, 0, "supply Twiss parameters, moments, floor coordinates for optimization?"},
    } ;

MATR matr_example;
/* name for matrix parameters */
PARAMETER matr_param[N_MATR_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&matr_example.length), NULL, 0.0, 0, "length"},
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&matr_example.filename), "", 0.0, 1, "input file"},
    {"ORDER", "", IS_LONG, 0, (long)((char *)&matr_example.order), NULL, 0.0, 1, "matrix order"},
    } ;

ALPH alph_example;
/* names for alpha magnet parameters */
PARAMETER alph_param[N_ALPH_PARAMS] = {
    {"XMAX", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.xmax), NULL, 0.0, 0, "size of alpha"},
    {"XS1", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.xs1), NULL, 0.0, 0, "inner scraper position"},
    {"XS2", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.xs2), NULL, 0.0, 0, "outer scraper position"},
    {"DP1", "", IS_DOUBLE, 1, (long)((char *)&alph_example.dp1), NULL, -1, 0, "inner scraper momentum deviation"},
    {"DP2", "", IS_DOUBLE, 1, (long)((char *)&alph_example.dp2), NULL, 1, 0, "outer scraper momentum deviation"},
    {"XPUCK", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.xPuck), NULL, -1, 0, "position of scraper puck"},
    {"WIDTHPUCK", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.widthPuck), NULL, 0, 0, "size of scraper puck"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&alph_example.dz), NULL, 0.0, 0, "misalignment"},
    {"TILT", "", IS_DOUBLE, 0, (long)((char *)&alph_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"PART", "", IS_LONG, 0, (long)((char *)&alph_example.part), NULL, 0.0, 0, "0=full, 1=first half, 2=second half"},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&alph_example.order), NULL, 0.0, 0, "matrix order [1,3]"}
    } ;

RFDF rfdf_example;
/* names for rf deflector parameters */
PARAMETER rfdf_param[N_RFDF_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&rfdf_example.length), NULL, 0.0, 0, "length"},
    {"PHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&rfdf_example.phase), NULL, 0.0, 0, "phase"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&rfdf_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"FREQUENCY", "HZ", IS_DOUBLE, 0, (long)((char *)&rfdf_example.frequency), NULL, DEFAULT_FREQUENCY, 0, "frequency"},
    {"VOLTAGE", "V", IS_DOUBLE, 0, (long)((char *)&rfdf_example.voltage), NULL, 0.0, 0, "voltage"},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&rfdf_example.time_offset), NULL, 0.0, 0, "time offset (adds to phase)"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&rfdf_example.n_kicks), NULL, 0.0, 1, "number of kicks"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&rfdf_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    } ;

RFTMEZ0 rftmez0_example;
/* names for tm mode cavity from Ez(z,r=0)
 */
PARAMETER rftmez0_param[N_RFTMEZ0_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&rftmez0_example.length), NULL, 0.0, 0, "length"},
    {"FREQUENCY", "HZ", IS_DOUBLE, 0, (long)((char *)&rftmez0_example.frequency), NULL, DEFAULT_FREQUENCY, 0, "frequency"},
    {"PHASE", "RAD", IS_DOUBLE, 0, (long)((char *)&rftmez0_example.phase), NULL, 0.0, 0, "phase"},
    {"EZ_PEAK", "V", IS_DOUBLE, 0, (long)((char *)&rftmez0_example.Ez_peak), NULL, 0.0, 0, "Peak on-axis longitudinal electric field"},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&rftmez0_example.time_offset), NULL, 0.0, 0, "time offset (adds to phase)"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&rftmez0_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link to other time-dependent elements)"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&rftmez0_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&rftmez0_example.dy), NULL, 0.0, 0, "misalignment"},
    {"N_STEPS", "", IS_LONG, 0, (long)((char *)&rftmez0_example.n_steps), NULL, 0.0, 100, "number of steps (for nonadaptive integration)"},
    {"RADIAL_ORDER", "", IS_LONG, 0, (long)((char *)&rftmez0_example.radial_order), NULL, 0.0, 1, "highest order in off-axis expansion"},
    {"CHANGE_P0", "", IS_LONG, 0, (long)((char *)&rftmez0_example.change_p0), NULL, 0.0, 0, "does element change central momentum?"},
    {"INPUTFILE", "", IS_STRING, 0, (long)((char *)&rftmez0_example.inputFile), NULL, 0.0, 0, "file containing Ez vs z at r=0"},
    {"ZCOLUMN", "", IS_STRING, 0, (long)((char *)&rftmez0_example.zColumn), NULL, 0.0, 0, "column containing z values"},
    {"EZCOLUMN", "", IS_STRING, 0, (long)((char *)&rftmez0_example.EzColumn), NULL, 0.0, 0, "column containing Ez values"},
    {"SOLENOID_FILE", "", IS_STRING, 0, (long)((char *)&rftmez0_example.solenoidFile), NULL, 0.0, 0, "file containing map of Bz and Br vs z and r.  Each page contains values for a single r."},
    {"SOLENOID_ZCOLUMN", "", IS_STRING, 0, (long)((char *)&rftmez0_example.solenoid_zColumn), NULL, 0.0, 0, "column containing z values for solenoid map."},
    {"SOLENOID_RCOLUMN", "", IS_STRING, 0, (long)((char *)&rftmez0_example.solenoid_rColumn), NULL, 0.0, 0, "column containing r values for solenoid map."},
    {"SOLENOID_BZCOLUMN", "", IS_STRING, 0, (long)((char *)&rftmez0_example.solenoidBzColumn), NULL, 0.0, 0, "column containing Bz values for solenoid map."},
    {"SOLENOID_BRCOLUMN", "", IS_STRING, 0, (long)((char *)&rftmez0_example.solenoidBrColumn), NULL, 0.0, 0, "column containing Br values for solenoid map."},
    {"SOLENOID_FACTOR", "", IS_DOUBLE, 0, (long)((char *)&rftmez0_example.solenoidFactor), NULL, 1.0, 0, "factor by which to multiply solenoid fields."},        
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&rftmez0_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&rftmez0_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, non-adaptive runge-kutta, modified midpoint)"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&rftmez0_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0, "{t|p},{median|min|max|ave|first|light} (e.g., \"t,median\")"},
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
    {"INTERVAL", "", IS_LONG, 0, (long)((char *)&watch_example.interval), NULL, 0.0, 1, "interval for data output (in turns)"},
    {"START_PASS", "", IS_LONG, 0, (long)((char*)&watch_example.start_pass), NULL, 0.0, 0, "pass on which to start"},
    {"FILENAME", "", IS_STRING, 0, (long)((char *)&watch_example.filename), "", 0.0, 0, "output filename"},
    {"LABEL", "", IS_STRING, 0, (long)((char *)&watch_example.label), "", 0.0, 0, "output label"},
    {"MODE", "", IS_STRING, 0, (long)((char *)&watch_example.mode), "coordinates", 0.0, 0, "coordinate, parameter, or centroid"},
    {"X_DATA", "", IS_LONG, 0, (long)((char*)&watch_example.xData), NULL, 0.0, 1, "include x data in coordinate mode?"},
    {"Y_DATA", "", IS_LONG, 0, (long)((char*)&watch_example.yData), NULL, 0.0, 1, "include y data in coordinate mode?"},
    {"LONGIT_DATA", "", IS_LONG, 0, (long)((char*)&watch_example.longitData), NULL, 0.0, 1, "include longitudinal data in coordinate mode?"},
    {"EXCLUDE_SLOPES", "", IS_LONG, 0, (long)((char*)&watch_example.excludeSlopes), NULL, 0.0, 0, "exclude slopes in coordinate mode?"},
    {"FLUSH_INTERVAL", "", IS_LONG, 0, (long)((char *)&watch_example.flushInterval), NULL, 0.0, 0, "file flushing interval (parameter or centroid mode)"},    
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
    {"DXP", "", IS_DOUBLE, 1, (long)((char *)&malign_example.dxp), NULL, 0.0, 0, "delta x'"},
    {"DYP", "", IS_DOUBLE, 1, (long)((char *)&malign_example.dyp), NULL, 0.0, 0, "delta y'"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&malign_example.dx), NULL, 0.0, 0, "delta x"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&malign_example.dy), NULL, 0.0, 0, "delta y"},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&malign_example.dz), NULL, 0.0, 0, "delta z"},
    {"DT", "S", IS_DOUBLE, 1, (long)((char *)&malign_example.dt), NULL, 0.0, 0, "delta t"},
    {"DP", "", IS_DOUBLE, 1, (long)((char *)&malign_example.dp), NULL, 0.0, 0, "delta p/pCentral"},
    {"DE", "", IS_DOUBLE, 1, (long)((char *)&malign_example.de), NULL, 0.0, 0, "delta gamma/gammaCentral"},
    {"ON_PASS", "", IS_LONG, 0, (long)((char *)&malign_example.on_pass), NULL, 0.0, -1, "pass on which to apply"},
    } ;

TW_LINAC twla_example;
/* names for traveling-wave linac parameters
 */
PARAMETER twla_param[N_TWLA_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&twla_example.length), NULL, 0.0, 0, "length"},
    {"FREQUENCY", "HZ", IS_DOUBLE, 0, (long)((char *)&twla_example.frequency), NULL, DEFAULT_FREQUENCY, 0, "frequency"},
    {"PHASE", "RAD", IS_DOUBLE, 0, (long)((char *)&twla_example.phase), NULL, 0.0, 0, "phase"},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&twla_example.time_offset), NULL, 0.0, 0, "time offset (adds to phase)"},
    {"EZ", "V/M", IS_DOUBLE, 0, (long)((char *)&twla_example.Ez), NULL, 0.0, 0, "electric field"},
    {"B_SOLENOID", "T", IS_DOUBLE, 0, (long)((char *)&twla_example.B_solenoid), NULL, 0.0, 0, "solenoid field"},
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&twla_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twla_example.x_max), NULL, 0.0, 0, "x half-aperture"},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twla_example.y_max), NULL, 0.0, 0, "y half-aperture"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&twla_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&twla_example.dy), NULL, 0.0, 0, "misalignment"},
    {"BETA_WAVE", "", IS_DOUBLE, 0, (long)((char *)&twla_example.beta_wave), NULL, DEFAULT_BETA_WAVE, 0, "(phase velocity)/c"},
    {"ALPHA", "1/M", IS_DOUBLE, 0, (long)((char *)&twla_example.alpha), NULL, 0.0, 0, "field attenuation factor"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&twla_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"N_STEPS", "", IS_LONG, 0, (long)((char *)&twla_example.n_steps), NULL, 0.0, 100, "number of steps (for nonadaptive integration)"},
    {"FOCUSSING", "", IS_LONG, 0, (long)((char *)&twla_example.focussing), NULL, 0.0, 1, "include focusing effects?"},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&twla_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, non-adaptive runge-kutta, modified midpoint)"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&twla_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0, "{t|p},{median|min|max|ave|first|light} (e.g., \"t,median\")"},
    {"CHANGE_P0", "", IS_LONG, 0, (long)((char *)&twla_example.change_p0), NULL, 0.0, 0, "does element change central momentum?"},
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
    {"ELLIPTICAL", "", IS_LONG, 0, (long)((char *)&maxamp_example.elliptical), NULL, 0.0, 0, "is aperture elliptical?"}
    } ;

ROTATE rotate_example;
/* names for beam rotation */
PARAMETER rotate_param[N_ROTATE_PARAMS] = {
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&rotate_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    } ;

/* names for transmission count */
PARAMETER trcount_param[N_TRCOUNT_PARAMS] = {
    {"DUMMY", "", IS_LONG, 0, 0, NULL, 0.0, 0, ""},
    } ;

/* names for recirculation point */
PARAMETER recirc_param[N_RECIRC_PARAMS] = {
    {"I_RECIRC_ELEMENT", "", IS_LONG, 0, 0, NULL, 0.0, 0, ""},
    } ;

QFRING qfring_example;
/* quadrupole fringe-field physical parameters */
PARAMETER qfring_param[N_QFRING_PARAMS]={
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&qfring_example.length), NULL, 0.0, 0, "length"},
    {"K1", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&qfring_example.k1), NULL, 0.0, 0, "peak geometric strength"},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&qfring_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&qfring_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&qfring_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&qfring_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "M", IS_DOUBLE, 1, (long)((char *)&qfring_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"DIRECTION", "", IS_LONG, 1, (long)((char *)&qfring_example.direction), NULL, 0.0, 0, "1=entrance, -1=exit"},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&qfring_example.order), NULL, 0.0, 0, "matrix order"}
    };

SCRAPER scraper_example;
/* scraper physical parameters */
PARAMETER scraper_param[N_SCRAPER_PARAMS]={
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&scraper_example.length), NULL, 0.0, 0, "length"},
    {"POSITION", "M", IS_DOUBLE, 0, (long)((char *)&scraper_example.position), NULL, 0.0, 0, "position of edge"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&scraper_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&scraper_example.dy), NULL, 0.0, 0, "misalignment"},
    {"XO", "M", IS_DOUBLE, 0, (long)((char *)&scraper_example.Xo), NULL, 0.0, 0, "radiation length"},
    {"INSERT_FROM", "", IS_STRING, 0, (long)((char *)&scraper_example.insert_from), NULL, 0.0, 0, "direction from which inserted (+x, -x, +y, -y"},
    {"ELASTIC", "", IS_LONG,  0, (long)((char *)&scraper_example.elastic), NULL, 0.0, 0, "elastic scattering?"},
    {"DIRECTION", "", IS_LONG,  0, (long)((char *)&scraper_example.direction), NULL, 0.0, -1, "obsolete"},
    };

CENTER center_example;
/* beam centering physical parameters */
PARAMETER center_param[N_CENTER_PARAMS]={
    {"X" , "", IS_LONG, 0, (long)((char *)&center_example.x), NULL, 0.0, 1, "center x coordinates?"},
    {"XP", "", IS_LONG, 0, (long)((char *)&center_example.xp), NULL, 0.0, 1, "center x' coordinates?"},
    {"Y" , "", IS_LONG, 0, (long)((char *)&center_example.y), NULL, 0.0, 1, "center y coordinates?" },
    {"YP", "", IS_LONG, 0, (long)((char *)&center_example.yp), NULL, 0.0, 1, "center y' coordinates?"},
    {"ONCE_ONLY", "", IS_LONG, 0, (long)((char *)&center_example.onceOnly), NULL, 0.0, 0, "compute centering offsets for first beam only, apply to all?"},
    };

KICKER kicker_example;
/* kicker physical parameters */
PARAMETER kicker_param[N_KICKER_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&kicker_example.length), NULL, 0.0, 0, "length"},
    {"ANGLE", "RAD", IS_DOUBLE, 0, (long)((char *)&kicker_example.angle), NULL, 0.0, 0, "kick angle"},
    {"TILT", "RAD", IS_DOUBLE, 0, (long)((char *)&kicker_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"TIME_OFFSET", "S", IS_DOUBLE, 0, (long)((char *)&kicker_example.time_offset), NULL, 0.0, 0, "time offset of waveform"},
    {"PERIODIC", "", IS_LONG, 0, (long)((char *)&kicker_example.periodic), NULL, 0.0, 0, "is waveform periodic?"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&kicker_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"FIRE_ON_PASS", "", IS_LONG, 0, (long)((char *)&kicker_example.fire_on_pass), NULL, 0.0, 0, "pass number to fire on"},
    {"WAVEFORM", "", IS_STRING, 0, (long)((char *)&kicker_example.waveform), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving kick factor vs time"},
    {"SPATIAL_DEPENDENCE", "", IS_STRING, 0, (long)((char *)&kicker_example.spatial_dependence), NULL, 0.0, 0, "rpn function of x and y giving spatial dependence of kick"},
    } ;

KSEXT ksext_example;
/* kick sextupole physical parameters */
PARAMETER ksext_param[N_KSEXT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&ksext_example.length), NULL, 0.0, 0, "length"},
    {"K2", "1/M$a3$n", IS_DOUBLE, 1, (long)((char *)&ksext_example.k2), NULL, 0.0, 0, "geometric strength"},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&ksext_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"BORE", "M", IS_DOUBLE, 1, (long)((char *)&ksext_example.bore), NULL, 0.0, 0, "bore radius"},
    {"B", "T", IS_DOUBLE, 1, (long)((char *)&ksext_example.B), NULL, 0.0, 0, "field at pole tip (used if bore nonzero)"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&ksext_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&ksext_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&ksext_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "M", IS_DOUBLE, 1, (long)((char *)&ksext_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&ksext_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks"},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&ksext_example.synch_rad), NULL, 0.0, 0, "include classical synchrotron radiation?"},
    {"SYSTEMATIC_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ksext_example.systematic_multipoles), NULL, 0.0, 0, "input file for systematic multipoles"},
    {"RANDOM_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&ksext_example.random_multipoles), NULL, 0.0, 0, "input file for random multipoles"},
    {"INTEGRATION_ORDER", "", IS_LONG, 0, (long)((char *)&ksext_example.integration_order), NULL, 0.0, 4, "integration order (2 or 4)"},
    };

KSBEND ksbend_example;
/* symplectic sector bending magnet physical parameters */
PARAMETER ksbend_param[N_KSBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, 1, (long)((char *)&ksbend_example.angle), NULL, 0.0, 0, "bend angle"},
    {"K1", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&ksbend_example.k1), NULL, 0.0, 0, "geometric quadrupole strength"},
    {"K2", "1/M$a3$n", IS_DOUBLE, 1, (long)((char *)&ksbend_example.k2), NULL, 0.0, 0, "geometric sextupole strength"},
    {"K3", "1/M$a3$n", IS_DOUBLE, 0, (long)((char *)&ksbend_example.k3), NULL, 0.0, 0, "geometric octupole strength"},
    {"K4", "1/M$a4$n", IS_DOUBLE, 0, (long)((char *)&ksbend_example.k4), NULL, 0.0, 0, "geometric decapole strength"},
    {"E1", "RAD", IS_DOUBLE, 1, (long)((char *)&ksbend_example.e1), NULL, 0.0, 0, "entrance edge angle"},
    {"E2", "RAD", IS_DOUBLE, 1, (long)((char *)&ksbend_example.e2), NULL, 0.0, 0, "exit edge angle"},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&ksbend_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"H1", "1/M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.h1), NULL, 0.0, 0, "entrance pole-face curvature"},
    {"H2", "1/M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.h2), NULL, 0.0, 0, "exit pole-face curvature"},
    {"HGAP", "M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.hgap), NULL, 0.0, 0, "half-gap between poles"},
    {"FINT", "", IS_DOUBLE, 1, (long)((char *)&ksbend_example.fint), NULL, DEFAULT_FINT, 0, "edge-field integral"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&ksbend_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, 1, (long)((char *)&ksbend_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"ETILT", "", IS_DOUBLE, 1, (long)((char *)&ksbend_example.etilt), NULL, 0.0, 0, "error rotation about incoming longitudinal axis"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&ksbend_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks"},
    {"NONLINEAR", "", IS_LONG, 0, (long)((char *)&ksbend_example.nonlinear), NULL, 0.0, 1, "include nonlinear field components?"},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&ksbend_example.synch_rad), NULL, 0.0, 0, "include classical synchrotron radiation?"},
    {"EDGE1_EFFECTS", "", IS_LONG, 1, (long)((char *)&ksbend_example.edge1_effects), NULL, 0.0, 1, "include entrace edge effects?"},
    {"EDGE2_EFFECTS", "", IS_LONG, 1, (long)((char *)&ksbend_example.edge2_effects), NULL, 0.0, 1, "include exit edge effects?"},
    {"EDGE_ORDER", "", IS_LONG, 1, (long)((char *)&ksbend_example.edge_order), NULL, 0.0, 1, "edge matrix order"},
    {"PARAXIAL", "", IS_LONG, 0, (long)((char *)&ksbend_example.paraxial), NULL, 0.0, 0, "use paraxial approximation?"},
    {"TRANSPORT", "", IS_LONG, 1, (long)((char *)&ksbend_example.TRANSPORT), NULL, 0.0, 0, "use (incorrect) TRANSPORT equations for T436 of edge?"},
    {"METHOD", "", IS_STRING, 0, (long)((char *)&ksbend_example.method), "modified-midpoint", 0.0, 0, "integration method (modified-midpoint, leap-frog"}
    };

KQUAD kquad_example;
/* kick quadrupole physical parameters */
PARAMETER kquad_param[N_KQUAD_PARAMS]={
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&kquad_example.length), NULL, 0.0, 0, "length"},
    {"K1", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&kquad_example.k1), NULL, 0.0, 0, "geometric strength"},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&kquad_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"BORE", "M", IS_DOUBLE, 1, (long)((char *)&kquad_example.bore), NULL, 0.0, 0, "bore radius"},
    {"B", "T", IS_DOUBLE, 1, (long)((char *)&kquad_example.B), NULL, 0.0, 0, "pole tip field (used if bore nonzero)"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&kquad_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&kquad_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&kquad_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "M", IS_DOUBLE, 1, (long)((char *)&kquad_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"N_KICKS", "", IS_LONG, 1, (long)((char *)&kquad_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks"},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&kquad_example.synch_rad), NULL, 0.0, 0, "include classical synchrotron radiation?"},
    {"SYSTEMATIC_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&kquad_example.systematic_multipoles), NULL, 0.0, 0, "input file for systematic multipoles"},
    {"RANDOM_MULTIPOLES", "", IS_STRING, 0, (long)((char *)&kquad_example.random_multipoles), NULL, 0.0, 0, "input file for random multipoles"},
    {"INTEGRATION_ORDER", "", IS_LONG, 0, (long)((char *)&kquad_example.integration_order), NULL, 0.0, 4, "integration order (2 or 4)"},
    };

MAGNIFY magnify_example;
/* magnifier physical parameters */
PARAMETER magnify_param[N_MAGNIFY_PARAMS] = {
    {"MX", "", IS_DOUBLE, 1, (long)((char *)&magnify_example.mx), NULL, 1.0, 0, "factor for x coordinates"},
    {"MXP", "", IS_DOUBLE, 1, (long)((char *)&magnify_example.mxp), NULL, 1.0, 0, "factor for x' coordinates"},
    {"MY", "", IS_DOUBLE, 1, (long)((char *)&magnify_example.my), NULL, 1.0, 0, "factor for y coordinates"},
    {"MYP", "", IS_DOUBLE, 1, (long)((char *)&magnify_example.myp), NULL, 1.0, 0, "factor for y' coordinates"},
    {"MS", "", IS_DOUBLE, 1, (long)((char *)&magnify_example.ms), NULL, 1.0, 0, "factor for s coordinates"},
    {"MDP", "", IS_DOUBLE, 1, (long)((char *)&magnify_example.mdp), NULL, 1.0, 0, "factor for (p-pCentral)/pCentral"},
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
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&hvcor_example.length), NULL, 0.0, 0, "length"},
    {"HKICK", "RAD", IS_DOUBLE, 1, (long)((char *)&hvcor_example.xkick), NULL, 0.0, 0, "x kick angle"},
    {"VKICK", "RAD", IS_DOUBLE, 1, (long)((char *)&hvcor_example.ykick), NULL, 0.0, 0, "y kick angle"},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&hvcor_example.tilt), NULL, 0.0, 0, "rotation about longitudinal axis"},
    {"B2", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&hvcor_example.b2), NULL, 0.0, 0, "normalized sextupole strength (e.g., kick = KICK*(1+B2*x^2))"},
    {"HCALIBRATION", "", IS_DOUBLE, 1, (long)((char *)&hvcor_example.xcalibration), NULL, 1.0, 0, "factor applied to obtain x kick"},
    {"VCALIBRATION", "", IS_DOUBLE, 1, (long)((char *)&hvcor_example.ycalibration), NULL, 1.0, 0, "factor applied to obtain y kick"},
    {"EDGE_EFFECTS", "", IS_LONG, 0, (long)((char *)&hvcor_example.edge_effects), NULL, 0.0, 0, "include edge effects?"},
    {"ORDER", "", IS_LONG, 1, (long)((char *)&hvcor_example.order), NULL, 0.0, 0, "matrix order"},
    {"STEERING", "", IS_LONG, 0, (long)((char *)&hvcor_example.steering), NULL, 0.0, 1, "use for steering?"},
    };

SCATTER scatter_example;
/* scatter physical parameters */
PARAMETER scatter_param[N_SCATTER_PARAMS] = {
    {"X", "M", IS_DOUBLE, 0, (long)((char*)&scatter_example.x), NULL, 0.0, 0, "rms scattering level for x"},
    {"XP", "M", IS_DOUBLE, 0, (long)((char*)&scatter_example.xp), NULL, 0.0, 0, "rms scattering level for x'"},
    {"Y", "M", IS_DOUBLE, 0, (long)((char*)&scatter_example.y), NULL, 0.0, 0, "rms scattering level for y"},
    {"YP", "M", IS_DOUBLE, 0, (long)((char*)&scatter_example.yp), NULL, 0.0, 0, "rms scattering level for y'"},
    {"DP", "M", IS_DOUBLE, 0, (long)((char*)&scatter_example.dp), NULL, 0.0, 0, "rms scattering level for (p-pCentral)/pCentral"},
    } ;
    
NIBEND nibend_example;
/* integrated bending magnet physical parameters */
PARAMETER nibend_param[N_NIBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, 1, (long)((char *)&nibend_example.angle), NULL, 0.0, 0, "bending angle"},
    {"E1", "RAD", IS_DOUBLE, 1, (long)((char *)&nibend_example.e1), NULL, 0.0, 0, "entrance edge angle"},
    {"E2", "RAD", IS_DOUBLE, 1, (long)((char *)&nibend_example.e2), NULL, 0.0, 0, "exit edge angle"},
    {"TILT", "", IS_DOUBLE, 1, (long)((char *)&nibend_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FINT", "", IS_DOUBLE, 1, (long)((char *)&nibend_example.fint), NULL, DEFAULT_FINT, 0, "edge-field integral"},
    {"HGAP", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.hgap), NULL, 0.0, 0, "half-gap between poles"},
    {"FP1", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.fp1), NULL, 10.0, 0, "fringe parameter (tanh model)"},
    {"FP2", "M", IS_DOUBLE, 1, (long)((char *)&nibend_example.fp2), NULL, 1.0, 0, "not used"},
    {"FSE", "", IS_DOUBLE, 1, (long)((char *)&nibend_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"ETILT", "", IS_DOUBLE, 1, (long)((char *)&nibend_example.etilt), NULL, 0.0, 0, "error rotation about incoming longitudinal axis"},
    {"ACCURACY", "", IS_DOUBLE, 1, (long)((char *)&nibend_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"MODEL", "", IS_STRING, 0, (long)((char *)&nibend_example.model), DEFAULT_NIBEND_TYPE, 0.0, 0, "fringe model (hard-edge, linear, cubic-spline, tanh, quintic"},
    {"METHOD", "", IS_STRING, 0, (long)((char *)&nibend_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, modified-midpoint, two-pass modified-midpoint, leap-frog, non-adaptive runge-kutta"},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&nibend_example.synch_rad), NULL, 0.0, 0, "include classical synchrotron radiation?"},
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

NISEPT nisept_example;
/* integrated septum magnet physical parameters */
PARAMETER nisept_param[N_NISEPT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&nisept_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, 1, (long)((char *)&nisept_example.angle), NULL, 0.0, 0, "bend angle"},
    {"E1", "RAD", IS_DOUBLE, 1, (long)((char *)&nisept_example.e1), NULL, 0.0, 0, "entrance edge angle"},
    {"B1", "1/M", IS_DOUBLE, 1, (long)((char *)&nisept_example.b1), NULL, 0.0, 0, "normalized gradient (K1=B1*L/ANGLE)"},
    {"Q1REF", "M", IS_DOUBLE, 1, (long)((char *)&nisept_example.q1_ref), NULL, 0.0, 0, "distance from septum at which bending radius is L/ANGLE"},
    {"FLEN", "M", IS_DOUBLE, 1, (long)((char *)&nisept_example.flen), NULL, 0.0, 0, "fringe field length"},
    {"ACCURACY", "", IS_DOUBLE, 1, (long)((char *)&nisept_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"METHOD", "", IS_STRING, 0, (long)((char *)&nisept_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, modified-midpoint, two-pass modified-midpoint, leap-frog, non-adaptive runge-kutta"},
    {"MODEL", "", IS_STRING, 0, (long)((char *)&nisept_example.model), DEFAULT_NIBEND_TYPE, 0.0, 0, "fringe model (hard-edge, linear, cubic-spline, tanh, quintic"},
    };

RAMPRF ramprf_example;
/* ramped rf cavity physical parameters */
PARAMETER ramprf_param[N_RAMPRF_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&ramprf_example.length), NULL, 0.0, 0, "length"},
    {"VOLT", "V", IS_DOUBLE, 0, (long)((char *)&ramprf_example.volt), NULL, 0.0, 0, "nominal voltage"},
    {"PHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&ramprf_example.phase), NULL, 0.0, 0, "nominal phase"},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&ramprf_example.freq), NULL, 500.0e6, 0, "nominal frequency"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&ramprf_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"VOLT_WAVEFORM", "", IS_STRING, 0, (long)((char *)&ramprf_example.vwaveform), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving voltage waveform factor vs time"},
    {"PHASE_WAVEFORM", "", IS_STRING, 0, (long)((char *)&ramprf_example.pwaveform), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving phase offset vs time (requires FREQ_WAVEFORM)"},
    {"FREQ_WAVEFORM", "", IS_STRING, 0, (long)((char *)&ramprf_example.fwaveform), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving frequency factor vs time (requires PHASE_WAVEFORM)"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&ramprf_example.fiducial), NULL, 0.0, 0, "mode for determining fiducial arrival time (light, tmean, first, pmaximum)"},
    };

/* momentum ramp physical parameters */
PARAMETER rampp_param[N_RAMPP_PARAMS] = {
    {"WAVEFORM", "", IS_STRING, 0, 0, NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving momentum factor vs time"}
    };

STRAY stray_example;
/* stray field physical parameters */
PARAMETER stray_param[N_STRAY_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&stray_example.length), NULL, 0.0, 0, "length"},
    {"LBX", "T", IS_DOUBLE, 1, (long)((char *)&stray_example.lBx), NULL, 0.0, 0, "local Bx"},
    {"LBY", "T", IS_DOUBLE, 1, (long)((char *)&stray_example.lBy), NULL, 0.0, 0, "local By"},
    {"GBX", "T", IS_DOUBLE, 1, (long)((char *)&stray_example.gBx), NULL, 0.0, 0, "global Bx"},
    {"GBY", "T", IS_DOUBLE, 1, (long)((char *)&stray_example.gBy), NULL, 0.0, 0, "global By"},
    {"GBZ", "T", IS_DOUBLE, 1, (long)((char *)&stray_example.gBz), NULL, 0.0, 0, "global Bz"},
    {"ORDER", "", IS_LONG, 0, (long)((char *)&stray_example.order), NULL, 0.0, 0, "matrix order"},
    };

CSBEND csbend_example;
/* canonically-integrated sector bending magnet physical parameters */
PARAMETER csbend_param[N_CSBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&csbend_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, 1, (long)((char *)&csbend_example.angle), NULL, 0.0, 0, "bend angle"},
    {"K1", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&csbend_example.k1), NULL, 0.0, 0, "geometric quadrupole strength"},
    {"K2", "1/M$a3$n", IS_DOUBLE, 1, (long)((char *)&csbend_example.k2), NULL, 0.0, 0, "geometric sextupole strength"},
    {"K3", "1/M$a3$n", IS_DOUBLE, 0, (long)((char *)&csbend_example.k3), NULL, 0.0, 0, "geometric octupole strength"},
    {"K4", "1/M$a4$n", IS_DOUBLE, 0, (long)((char *)&csbend_example.k4), NULL, 0.0, 0, "geometric decapole strength"},
    {"E1", "RAD", IS_DOUBLE, 1, (long)((char *)&csbend_example.e1), NULL, 0.0, 0, "entrance edge angle"},
    {"E2", "RAD", IS_DOUBLE, 1, (long)((char *)&csbend_example.e2), NULL, 0.0, 0, "exit edge angle"},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&csbend_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"H1", "1/M", IS_DOUBLE, 1, (long)((char *)&csbend_example.h1), NULL, 0.0, 0, "entrance pole-face curvature"},
    {"H2", "1/M", IS_DOUBLE, 1, (long)((char *)&csbend_example.h2), NULL, 0.0, 0, "exit pole-face curvature"},
    {"HGAP", "M", IS_DOUBLE, 1, (long)((char *)&csbend_example.hgap), NULL, 0.0, 0, "half-gap between poles"},
    {"FINT", "", IS_DOUBLE, 1, (long)((char *)&csbend_example.fint), NULL, DEFAULT_FINT, 0, "edge-field integral"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&csbend_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&csbend_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&csbend_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, 1, (long)((char *)&csbend_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"ETILT", "", IS_DOUBLE, 1, (long)((char *)&csbend_example.etilt), NULL, 0.0, 0, "error rotation about incoming longitudinal axis"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&csbend_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks"},
    {"NONLINEAR", "", IS_LONG, 0, (long)((char *)&csbend_example.nonlinear), NULL, 0.0, 1, "include nonlinear field components?"},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&csbend_example.synch_rad), NULL, 0.0, 0, "include classical synchrotron radiation?"},
    {"EDGE1_EFFECTS", "", IS_LONG, 0, (long)((char *)&csbend_example.edge1_effects), NULL, 0.0, 1, "include entrace edge effects?"},
    {"EDGE2_EFFECTS", "", IS_LONG, 0, (long)((char *)&csbend_example.edge2_effects), NULL, 0.0, 1, "include exit edge effects?"},
    {"INTEGRATION_ORDER", "", IS_LONG, 0, (long)((char *)&csbend_example.integration_order), NULL, 0.0, 2, "integration order (2 or 4)"},
    {"EDGE1_KICK_LIMIT", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.edge1_kick_limit), NULL, -1., 0, "maximum kick entrance edge can deliver"},
    {"EDGE2_KICK_LIMIT", "", IS_DOUBLE, 0, (long)((char *)&csbend_example.edge2_kick_limit), NULL, -1., 0, "maximum kick exit edge can deliver"},
    {"KICK_LIMIT_SCALING", "", IS_LONG, 0, (long)((char *)&csbend_example.kick_limit_scaling), NULL, 0, 0, "scale maximum edge kick with FSE?"},
    {"USE_BN", "", IS_LONG, 0, (long)((char *)&csbend_example.use_bn), NULL, 0.0, 0, "use B<n> instead of K<n>?"},
    {"B1", "1/M", IS_DOUBLE, 0, (long)((char *)&csbend_example.b1), NULL, 0.0, 0, "K1 = B1*rho, where rho is bend radius"},
    {"B2", "1/M$a2$n", IS_DOUBLE, 0, (long)((char *)&csbend_example.b2), NULL, 0.0, 0, "K2 = B2*rho"},
    {"B3", "1/M$a3$n", IS_DOUBLE, 0, (long)((char *)&csbend_example.b3), NULL, 0.0, 0, "K3 = B3*rho"},
    {"B4", "1/M$a4$n", IS_DOUBLE, 0, (long)((char *)&csbend_example.b4), NULL, 0.0, 0, "K4 = B4*rho"},
    {"ISR", "", IS_LONG, 0, (long)((char *)&csbend_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (scattering)?"},
    };


CSRCSBEND csrcsbend_example;
/* canonically-integrated sector bending magnet with CSR physical parameters */
PARAMETER csrcsbend_param[N_CSRCSBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.angle), NULL, 0.0, 0, "bend angle"},
    {"K1", "1/M$a2$n", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.k1), NULL, 0.0, 0, "geometric quadrupole strength"},
    {"K2", "1/M$a3$n", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.k2), NULL, 0.0, 0, "geometric sextupole strength"},
    {"K3", "1/M$a3$n", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.k3), NULL, 0.0, 0, "geometric octupole strength"},
    {"K4", "1/M$a4$n", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.k4), NULL, 0.0, 0, "geometric decapole strength"},
    {"E1", "RAD", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.e1), NULL, 0.0, 0, "entrance edge angle"},
    {"E2", "RAD", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.e2), NULL, 0.0, 0, "exit edge angle"},
    {"TILT", "RAD", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.tilt), NULL, 0.0, 0, "rotation about incoming longitudinal axis"},
    {"H1", "1/M", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.h1), NULL, 0.0, 0, "entrance pole-face curvature"},
    {"H2", "1/M", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.h2), NULL, 0.0, 0, "exit pole-face curvature"},
    {"HGAP", "M", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.hgap), NULL, 0.0, 0, "half-gap between poles"},
    {"FINT", "", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.fint), NULL, DEFAULT_FINT, 0, "edge-field integral"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.dy), NULL, 0.0, 0, "misalignment"},
    {"DZ", "M", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.dz), NULL, 0.0, 0, "misalignment"},
    {"FSE", "", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"ETILT", "", IS_DOUBLE, 1, (long)((char *)&csrcsbend_example.etilt), NULL, 0.0, 0, "error rotation about incoming longitudinal axis"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.n_kicks), NULL, 0.0, DEFAULT_N_KICKS, "number of kicks"},
    {"NONLINEAR", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.nonlinear), NULL, 0.0, 1, "include nonlinear field components?"},
    {"SYNCH_RAD", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.synch_rad), NULL, 0.0, 0, "include classical synchrotron radiation?"},
    {"EDGE1_EFFECTS", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.edge1_effects), NULL, 0.0, 1, "include entrace edge effects?"},
    {"EDGE2_EFFECTS", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.edge2_effects), NULL, 0.0, 1, "include exit edge effects?"},
    {"INTEGRATION_ORDER", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.integration_order), NULL, 0.0, 2, "integration order (2 or 4)"},
    {"BINS", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.bins), NULL, 0.0, 0, "number of bins for CSR wake"},
    {"BIN_RANGE_FACTOR", "", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.binRangeFactor), NULL, 1.2, 0, "Factor by which to increase the range of histogram compared to total bunch length.  Large value eliminates binning problems in CSRDRIFTs."},
    {"SG_HALFWIDTH", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.SGHalfWidth), NULL, 0.0, 0, "Savitzky-Golay filter half-width for smoothing current histogram"},
    {"SG_ORDER", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing current histogram"},
    {"SGDERIV_HALFWIDTH", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.SGDerivHalfWidth), NULL, 0.0, 0, "Savitzky-Golay filter half-width for taking derivative of current histogram"},
    {"SGDERIV_ORDER", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.SGDerivOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for taking derivative of current histogram"},
    {"OUTPUT_FILE", "", IS_STRING, 0, (long)((char *)&csrcsbend_example.histogramFile), NULL, 0.0, 0, "output file for CSR wakes"},
    {"OUTPUT_INTERVAL", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.outputInterval), NULL, 0.0, 1, "interval (in kicks) of output to OUTPUT_FILE"},
    {"OUTPUT_LAST_WAKE_ONLY", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.outputLastWakeOnly), NULL, 0.0, 0, "output final wake only?"},
    {"STEADY_STATE", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.steadyState), NULL, 0.0, 0, "use steady-state wake equations?"},
    {"USE_BN", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.use_bn), NULL, 0.0, 0, "use B<n> instead of K<n>?"},
    {"B1", "1/M", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.b1), NULL, 0.0, 0, "K1 = B1*rho, where rho is bend radius"},
    {"B2", "1/M$a2$n", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.b2), NULL, 0.0, 0, "K2 = B2*rho"},
    {"B3", "1/M$a3$n", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.b3), NULL, 0.0, 0, "K3 = B3*rho"},
    {"B4", "1/M$a4$n", IS_DOUBLE, 0, (long)((char *)&csrcsbend_example.b4), NULL, 0.0, 0, "K4 = B4*rho"},
    {"ISR", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.isr), NULL, 0.0, 0, "include incoherent synchrotron radiation (scattering)?"},
    {"CSR", "", IS_LONG, 0, (long)((char *)&csrcsbend_example.csr), NULL, 0.0, 1, "enable CSR computations?"},
    };

TUBEND tubend_example;
/* special bending magnet for top-up with entry through the side! */
PARAMETER tubend_param[N_TUBEND_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&tubend_example.length), NULL, 0.0, 0, "arc length"},
    {"ANGLE", "RAD", IS_DOUBLE, 1, (long)((char *)&tubend_example.angle), NULL, 0.0, 0, "bend angle"},
    {"FSE", "", IS_DOUBLE, 1, (long)((char *)&tubend_example.fse), NULL, 0.0, 0, "fractional strength error"},
    {"OFFSET", "", IS_DOUBLE, 1, (long)((char *)&tubend_example.offset), NULL, 0.0, 0, "horizontal offset of magnet center from arc center"},
    {"MAGNET_WIDTH", "", IS_DOUBLE, 1, (long)((char *)&tubend_example.magnet_width), NULL, 0.0, 0, "horizontal width of the magnet pole"},
    {"MAGNET_ANGLE", "", IS_DOUBLE, 1, (long)((char *)&tubend_example.magnet_angle), NULL, 0.0, 0, "angle that the magnet was designed for"},
    };

TWMTA twmta_example;
/* names for traveling-wave muffin-tin linac parameters
 */
PARAMETER twmta_param[N_TWMTA_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&twmta_example.length), NULL, 0.0, 0, "length"},
    {"FREQUENCY", "HZ", IS_DOUBLE, 0, (long)((char *)&twmta_example.frequency), NULL, DEFAULT_FREQUENCY, 0, "frequency"},
    {"PHASE", "RAD", IS_DOUBLE, 0, (long)((char *)&twmta_example.phase), NULL, 0.0, 0, "phase"},
    {"EZ", "V/M", IS_DOUBLE, 0, (long)((char *)&twmta_example.Ez), NULL, 0.0, 0, "electric field"},
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&twmta_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"X_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twmta_example.x_max), NULL, 0.0, 0, "x half-aperture"},
    {"Y_MAX", "M", IS_DOUBLE, 0, (long)((char *)&twmta_example.y_max), NULL, 0.0, 0, "y half-aperture"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&twmta_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&twmta_example.dy), NULL, 0.0, 0, "misalignment"},
    {"KX", "1/M", IS_DOUBLE, 0, (long)((char *)&twmta_example.kx), NULL, 0.0, 0, "horizontal wave number"},
    {"BETA_WAVE", "", IS_DOUBLE, 0, (long)((char *)&twmta_example.beta_wave), NULL, DEFAULT_BETA_WAVE, 0, "(phase velocity)/c"},
    {"BSOL", "", IS_DOUBLE, 0, (long)((char *)&twmta_example.Bsol), NULL, 0.0, 0, "solenoid field"},
    {"ALPHA", "1/M", IS_DOUBLE, 0, (long)((char *)&twmta_example.alpha), NULL, 0.0, 0, "field attenuation factor"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&twmta_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"N_STEPS", "", IS_LONG, 0, (long)((char *)&twmta_example.n_steps), NULL, 0.0, 100, "number of kicks"},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&twmta_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, non-adaptive runge-kutta, modified midpoint)"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&twmta_example.fiducial), DEFAULT_FIDUCIAL_MODE, 0.0, 0, "{t|p},{median|min|max|ave|first|light} (e.g., \"t,median\")"}
    } ;

MATTER matter_example;
/* matter physical parameters */
PARAMETER matter_param[N_MATTER_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&matter_example.length), NULL, 0.0, 0, "length"},
    {"XO", "M", IS_DOUBLE, 1, (long)((char *)&matter_example.Xo), NULL, 0.0, 0, "radiation length"},
    {"ELASTIC", "", IS_LONG, 0, (long)((char *)&matter_example.elastic), NULL, 0.0, 0, "elastic scattering?"},
    };

RFMODE rfmode_example;
/* RFMODE physical parameters */
PARAMETER rfmode_param[N_RFMODE_PARAMS] = {
    {"RA", "Ohm", IS_DOUBLE, 0, (long)((char *)&rfmode_example.Ra), NULL, 0.0, 0, "shunt impedance"},
    {"RS", "Ohm", IS_DOUBLE, 0, (long)((char *)&rfmode_example.Rs), NULL, 0.0, 0, "shunt impedance (Ra=2*Rs)"},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&rfmode_example.Q), NULL, 0.0, 1, "cavity Q"},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&rfmode_example.freq), NULL, 0.0, 0, "frequency"},
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&rfmode_example.charge), NULL, 0.0, 0, "beam charge (or use CHARGE element)"},
    {"INITIAL_V", "V", IS_DOUBLE, 0, (long)((char *)&rfmode_example.initial_V), NULL, 0.0, 0, "initial voltage"},
    {"INITIAL_PHASE", "RAD", IS_DOUBLE, 0, (long)((char *)&rfmode_example.initial_phase), NULL, 0.0, 0, "initial phase"},
    {"INITIAL_T", "S", IS_DOUBLE, 0, (long)((char *)&rfmode_example.initial_t), NULL, 0.0, 0, "time at which INITIAL_V and INITIAL_PHASE held"},
    {"BETA", "", IS_DOUBLE, 0, (long)((char *)&rfmode_example.beta), NULL, 0.0, 0, "normalized load impedance"},
    {"BIN_SIZE", "S", IS_DOUBLE, 0, (long)((char *)&rfmode_example.bin_size), NULL, 0.0, 0, "bin size for current histogram (use 0 for autosize)"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&rfmode_example.n_bins), NULL, 0.0, 20, "number of bins for current histogram"},
    {"PRELOAD", "", IS_LONG, 0, (long)((char *)&rfmode_example.preload), NULL, 0.0, 0, "preload cavity with steady-state field"},
    {"PRELOAD_FACTOR", "", IS_DOUBLE, 0, (long)((char *)&rfmode_example.preload_factor), NULL, 1.0, 0, "multiply preloaded field by this value"},
    {"RIGID_UNTIL_PASS", "", IS_LONG, 0, (long)((char *)&rfmode_example.rigid_until_pass), NULL, 0.0, 0, "don't affect the beam until this pass"},
    {"SAMPLE_INTERVAL", "", IS_LONG, 0, (long)((char *)&rfmode_example.sample_interval), NULL, 0.0, 1, "passes between output to RECORD file"},
    {"RECORD", "", IS_STRING, 0, (long)((char *)&rfmode_example.record), NULL, 0.0, 0, "output file for cavity fields"},
    {"SINGLE_PASS", "", IS_LONG, 0, (long)((char *)&rfmode_example.single_pass), NULL, 0.0, 0, "if nonzero, don't accumulate field from pass to pass"},
    {"PASS_INTERVAL", "", IS_LONG, 0, (long)((char *)&rfmode_example.pass_interval), NULL, 0.0, 1, "interval in passes at which to apply PASS_INTERVAL times the field (may increase speed)"},
    };

TRFMODE trfmode_example;
/* TRFMODE physical parameters */
PARAMETER trfmode_param[N_TRFMODE_PARAMS] = {
    {"RA", "Ohm", IS_DOUBLE, 0, (long)((char *)&trfmode_example.Ra), NULL, 0.0, 0, "shunt impedance"},
    {"RS", "Ohm", IS_DOUBLE, 0, (long)((char *)&trfmode_example.Rs), NULL, 0.0, 0, "shunt impedance (Ra=2*Rs)"},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&trfmode_example.Q), NULL, 0.0, 1, "cavity Q"},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&trfmode_example.freq), NULL, 0.0, 0, "frequency"},
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&trfmode_example.charge), NULL, 0.0, 0, "beam charge (or use CHARGE element)"},
    {"BETA", "", IS_DOUBLE, 0, (long)((char *)&trfmode_example.beta), NULL, 0.0, 0, "normalized load impedance"},
    {"BIN_SIZE", "S", IS_DOUBLE, 0, (long)((char *)&trfmode_example.bin_size), NULL, 0.0, 0, "bin size for current histogram (use 0 for autosize)"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&trfmode_example.n_bins), NULL, 0.0, 20, "number of bins for current histogram"},
    {"PLANE", "", IS_STRING, 0, (long)((char *)&trfmode_example.plane), "both", 0.0, 0, "x, y, or both"},
    {"SINGLE_PASS", "", IS_LONG, 0, (long)((char *)&trfmode_example.single_pass), NULL, 0.0, 0, "if nonzero, don't accumulate field from pass to pass"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&trfmode_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&trfmode_example.dy), NULL, 0.0, 0, "misalignment"},
    };

ZLONGIT zlongit_example;
/* ZLONGIT physical parameters */
PARAMETER zlongit_param[N_ZLONGIT_PARAMS] = {
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&zlongit_example.charge), NULL, 0.0, 0, "beam charge (or use CHARGE element)"},
    {"BROAD_BAND", "", IS_LONG, 0, (long)((char *)&zlongit_example.broad_band), NULL, 0.0, 0, "broad-band impedance?"},
    {"RA", "Ohm", IS_DOUBLE, 0, (long)((char *)&zlongit_example.Ra), NULL, 0.0, 0, "shunt impedance"},
    {"RS", "Ohm", IS_DOUBLE, 0, (long)((char *)&zlongit_example.Rs), NULL, 0.0, 0, "shunt impedance (Ra=2*Rs)"},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&zlongit_example.Q), NULL, 0.0, 1, "cavity Q"},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&zlongit_example.freq), NULL, 0.0, 0, "frequency (BROAD_BAND=1)"},
    {"ZREAL", "", IS_STRING, 0, (long)((char *)&zlongit_example.Zreal), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving real part of impedance vs f (BROAD_BAND=0)"},
    {"ZIMAG", "", IS_STRING, 0, (long)((char *)&zlongit_example.Zimag), NULL, 0.0, 0, "<filename>=<x>+<y> form specification of input file giving imaginary part of impedance vs f (BROAD_BAND=0)"},
    {"BIN_SIZE", "S", IS_DOUBLE, 0, (long)((char *)&zlongit_example.bin_size), NULL, 0.0, 0, "bin size for current histogram (use 0 for autosize)"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&zlongit_example.n_bins), NULL, 0.0, 128, "number of bins for current histogram"},
    {"WAKES", "", IS_STRING, 0, (long)((char *)&zlongit_example.wakes), NULL, 0.0, 0, "filename for output of wake"},
    {"WAKE_INTERVAL", "", IS_LONG, 0, (long)((char *)&zlongit_example.wake_interval), NULL, 0.0, 1, "interval in passes at which to output wake"},
    {"AREA_WEIGHT", "", IS_LONG, 0, (long)((char *)&zlongit_example.area_weight), NULL, 0.0, 0, "use area-weighting in assigning charge to histogram?"},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&zlongit_example.interpolate), NULL, 0.0, 0, "interpolate wake?"},
    {"SMOOTHING", "", IS_LONG, 0, (long)((char *)&zlongit_example.smoothing), NULL, 0.0, 0, "smooth current histogram?"},
    {"SG_ORDER", "", IS_LONG, 0, (long)((char *)&zlongit_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing"},
    {"SG_HALFWIDTH", "", IS_LONG, 0, (long)((char *)&zlongit_example.SGHalfWidth), NULL, 0.0, 4, "Savitzky-Golay filter halfwidth for smoothing"},
    };

SREFFECTS SReffects_example;
/* SREFFECTS physical parameters */
PARAMETER sreffects_param[N_SREFFECTS_PARAMS] = {
    {"JX", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.Jx), NULL, 1.0, 0, "x damping partition number"},
    {"JY", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.Jy), NULL, 1.0, 0, "y damping partition number"},
    {"JDELTA", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.Jdelta), NULL, 2.0, 0, "momentum damping partition number"},
    {"EXREF", "m", IS_DOUBLE, 0, (long)((char *)&SReffects_example.exRef), NULL, 0.0, 0, "reference equilibrium x emittance"},
    {"EYREF", "m", IS_DOUBLE, 0, (long)((char *)&SReffects_example.eyRef), NULL, 0.0, 0, "reference equilibrium y emittance"},
    {"SDELTAREF", "m", IS_DOUBLE, 0, (long)((char *)&SReffects_example.SdeltaRef), NULL, 0.0, 0, "reference equilibrium fractional momentum spread"},
    {"DDELTAREF", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.DdeltaRef), NULL, 0.0, 0, "reference fractional momentum loss (per turn)"},
    {"PREF", "m$be$nc", IS_DOUBLE, 0, (long)((char *)&SReffects_example.pRef), NULL, 0.0, 0, "reference momentum (to which other reference values pertain)"},
    {"COUPLING", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.coupling), NULL, 0.0, 0, "x-y coupling"},
    {"FRACTION", "", IS_DOUBLE, 0, (long)((char *)&SReffects_example.fraction), NULL, 1.0, 0, "fraction of implied SR effect to simulate with each instance"},
    {"DAMPING", "", IS_LONG, 0, (long)((char *)&SReffects_example.damp), NULL, 0, 1, "include damping?"},
    {"QEXCITATION", "", IS_LONG, 0, (long)((char *)&SReffects_example.qExcite), NULL, 0, 1, "include quantum excitation?"},
    {"LOSSES", "", IS_LONG, 0, (long)((char *)&SReffects_example.loss), NULL, 0, 1, "include average losses?"},
    };

MODRF modrf_example;
PARAMETER modrf_param[N_MODRF_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&modrf_example.length), NULL, 0.0, 0, "length"},
    {"VOLT", "V", IS_DOUBLE, 0, (long)((char *)&modrf_example.volt), NULL, 0.0, 0, "nominal voltage"},
    {"PHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&modrf_example.phase), NULL, 0.0, 0, "nominal phase"},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&modrf_example.freq), NULL, 500.0e6, 0, "nominal frequency"},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&modrf_example.Q), NULL, 0.0, 0, "cavity Q"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&modrf_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"AMMAG", "", IS_DOUBLE, 0, (long)((char *)&modrf_example.amMag), NULL, 0.0, 0, "magnitude of amplitude modulation"},
    {"AMPHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&modrf_example.amPhase), NULL, 0.0, 0, "phase of amplitude modulation"},
    {"AMFREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&modrf_example.amFreq), NULL, 0.0, 0, "frequency of amplitude modulation"},
    {"PMMAG", "DEG", IS_DOUBLE, 0, (long)((char *)&modrf_example.pmMag), NULL, 0.0, 0, "magnitude of phase modulation"},
    {"PMPHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&modrf_example.pmPhase), NULL, 0.0, 0, "phase of phase modulation"},
    {"PMFREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&modrf_example.pmFreq), NULL, 0.0, 0, "frequency of phase modulation"},
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&modrf_example.fiducial), NULL, 0.0, 0, "mode for determining fiducial arrival time (light, tmean, first, pmaximum)"},
    };    

BMAPXY bmapxy_example;
PARAMETER bmapxy_param[N_BMAPXY_PARAMS] = {
{"L", "M", IS_DOUBLE, 0, (long)((char *)&bmapxy_example.length), NULL, 0.0, 0, "length"},
{"STRENGTH", NULL, IS_DOUBLE, 0, (long)((char *)&bmapxy_example.strength), NULL, 0.0, 0, "factor by which to multiply field"},
{"ACCURACY", NULL, IS_DOUBLE, 0, (long)((char *)&bmapxy_example.accuracy), NULL, 0.0, 0, "integration accuracy"},
{"METHOD", NULL, IS_STRING, 0, (long)((char*)&bmapxy_example.method), NULL, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, modified-midpoint, two-pass modified-midpoint, leap-frog, non-adaptive runge-kutta"},
{"FILENAME", NULL, IS_STRING, 0, (long)((char*)&bmapxy_example.filename), NULL, 0.0, 0, "name of file containing columns (x, y, Fx, Fy) giving normalized field (Fx, Fy) vs (x, y)"},
};  

ZTRANSVERSE ztransverse_example;
PARAMETER ztransverse_param[N_ZTRANSVERSE_PARAMS] = {
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.charge), NULL, 0.0, 0, "beam charge (or use CHARGE element)"},
    {"BROAD_BAND", "", IS_LONG, 0, (long)((char *)&ztransverse_example.broad_band), NULL, 0.0, 0, "broad-band impedance?"},
    {"RS", "Ohm", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.Rs), NULL, 0.0, 0, "shunt impedance (Ra=2*Rs)"},
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
    {"SMOOTHING", "", IS_LONG, 0, (long)((char *)&ztransverse_example.smoothing), NULL, 0.0, 0, "smooth current histogram?"},
    {"SG_ORDER", "", IS_LONG, 0, (long)((char *)&ztransverse_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing"},
    {"SG_HALFWIDTH", "", IS_LONG, 0, (long)((char *)&ztransverse_example.SGHalfWidth), NULL, 0.0, 4, "Savitzky-Golay filter halfwidth for smoothing"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&ztransverse_example.dy), NULL, 0.0, 0, "misalignment"},
};

IBSCATTER ibs_example;
PARAMETER ibscatter_param[N_IBSCATTER_PARAMS] = {
  {"COUPLING", "", IS_DOUBLE, 0, (long)((char *)&ibs_example.coupling), NULL, 1.0, 0, "x-y coupling"},
  {"FRACTION", "", IS_DOUBLE, 0, (long)((char *)&ibs_example.fraction), NULL, 1.0, 0, "fraction of full effect from each instance"},
  {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&ibs_example.charge), NULL, 0.0, 0, "beam charge (or use CHARGE element)"},
};

WAKE wake_example;
/* WAKE physical parameters */
PARAMETER wake_param[N_WAKE_PARAMS] = {
    {"INPUTFILE", "", IS_STRING, 0, (long)((char *)&wake_example.inputFile), NULL, 0.0, 0, "name of file giving Green function"},
    {"TCOLUMN", "", IS_STRING, 0, (long)((char *)&wake_example.tColumn), NULL, 0.0, 0, "column in INPUTFILE containing time data"},
    {"WCOLUMN", "", IS_STRING, 0, (long)((char *)&wake_example.WColumn), NULL, 0.0, 0, "column in INPUTFILE containing Green function"},
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&wake_example.charge), NULL, 0.0, 0, "beam charge (or use CHARGE element)"},
    {"FACTOR", "C", IS_DOUBLE, 0, (long)((char *)&wake_example.factor), NULL, 1.0, 0, "factor to multiply wake by"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&wake_example.n_bins), NULL, 0.0, 128, "number of bins for current histogram"},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&wake_example.interpolate), NULL, 0.0, 0, "interpolate wake?"},
    {"SMOOTHING", "", IS_LONG, 0, (long)((char *)&wake_example.smoothing), NULL, 0.0, 0, "smooth current histogram?"},
    {"SG_HALFWIDTH", "", IS_LONG, 0, (long)((char *)&wake_example.SGHalfWidth), NULL, 0.0, 4, "Savitzky-Golay filter half-width for smoothing"},
    {"SG_ORDER", "", IS_LONG, 0, (long)((char *)&wake_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing"},
    {"CHANGE_P0", "", IS_LONG, 0, (long)((char *)&wake_example.change_p0), NULL, 0.0, 0, "change central momentum?"},
    };

TRWAKE trwake_example;
/* TRWAKE physical parameters */
PARAMETER trwake_param[N_TRWAKE_PARAMS] = {
    {"INPUTFILE", "", IS_STRING, 0, (long)((char *)&trwake_example.inputFile), NULL, 0.0, 0, "name of file giving Green functions"},
    {"TCOLUMN", "", IS_STRING, 0, (long)((char *)&trwake_example.tColumn), NULL, 0.0, 0, "column in INPUTFILE containing time data"},
    {"WXCOLUMN", "", IS_STRING, 0, (long)((char *)&trwake_example.WxColumn), NULL, 0.0, 0, "column in INPUTFILE containing x Green function"},
    {"WYCOLUMN", "", IS_STRING, 0, (long)((char *)&trwake_example.WyColumn), NULL, 0.0, 0, "column in INPUTFILE containing y Green function"},
    {"CHARGE", "C", IS_DOUBLE, 0, (long)((char *)&trwake_example.charge), NULL, 0.0, 0, "beam charge (or use CHARGE element)"},
    {"FACTOR", "C", IS_DOUBLE, 0, (long)((char *)&trwake_example.factor), NULL, 1.0, 0, "factor to multiply wake by"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&trwake_example.n_bins), NULL, 0.0, 128, "number of bins for current histogram"},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&trwake_example.interpolate), NULL, 0.0, 0, "interpolate wake?"},
    {"SMOOTHING", "", IS_LONG, 0, (long)((char *)&trwake_example.smoothing), NULL, 0.0, 0, "smooth current histogram?"},
    {"SG_HALFWIDTH", "", IS_LONG, 0, (long)((char *)&trwake_example.SGHalfWidth), NULL, 0.0, 4, "Savitzky-Golay filter half-width for smoothing"},
    {"SG_ORDER", "", IS_LONG, 0, (long)((char *)&trwake_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&trwake_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&trwake_example.dy), NULL, 0.0, 0, "misalignmnet"},
    };

CHARGE charge_example;
/* CHARGE physical parameters */
PARAMETER charge_param[N_CHARGE_PARAMS] = {
    {"TOTAL", "C", IS_DOUBLE, 0, (long)((char *)&charge_example.charge), NULL, 0.0, 0, "total charge in beam"},
    {"PER_PARTICLE", "C", IS_DOUBLE, 0, (long)((char *)&charge_example.chargePerParticle), NULL, 0.0, 0, "charge per macroparticle"},
};

PFILTER pfilter_example;
/* PFILTER physical parameters */
PARAMETER pfilter_param[N_PFILTER_PARAMS] = {
    {"DELTALIMIT", "", IS_DOUBLE, 0, (long)((char *)&pfilter_example.deltaLimit), NULL, -1.0, 0, "maximum fractional momentum deviation"},
    {"LOWERFRACTION", "", IS_DOUBLE, 0, (long)((char *)&pfilter_example.lowerFraction), NULL, 0.0, 0, "fraction of lowest-momentum particles to remove"},
    {"UPPERFRACTION", "", IS_DOUBLE, 0, (long)((char *)&pfilter_example.upperFraction), NULL, 0.0, 0, "fraction of highest-momentum particles to remove"},
    {"FIXPLIMITS", "", IS_LONG, 0, (long)((char *)&pfilter_example.fixPLimits), NULL, 0.0, 0, "fix the limits in p from LOWERFRACTION and UPPERFRACTION applied to first beam"},
    {"BEAMCENTERED", "", IS_LONG, 0, (long)((char *)&pfilter_example.beamCentered), NULL, 0.0, 0, "if nonzero, center for DELTALIMIT is average beam momentum"},
};

HISTOGRAM histogram_example;
/* HISTOGRAM physical parameters */
PARAMETER histogram_param[N_HISTOGRAM_PARAMS] = {
  {"FILENAME", "", IS_STRING, 0, (long)((char *)&histogram_example.filename), "", 0.0, 0, "filename for histogram output"},
  {"INTERVAL", "", IS_LONG, 0, (long)((char *)&histogram_example.interval), NULL, 0.0, 1, "interval in passes between output"},
  {"START_PASS", "", IS_LONG, 0, (long)((char*)&histogram_example.startPass), NULL, 0.0, 0, "starting pass for output"},
  {"BINS", "", IS_LONG, 0, (long)((char*)&histogram_example.bins), NULL, 0.0, 50, "number of bins"},
  {"FIXED_BIN_SIZE", "", IS_LONG, 0, (long)((char*)&histogram_example.fixedBinSize), NULL, 0.0, 0, "if nonzero, bin size is fixed at given value"},
  {"X_DATA", "", IS_LONG, 0, (long)((char*)&histogram_example.xData), NULL, 0.0, 1, "histogram x and x'?"},
  {"Y_DATA", "", IS_LONG, 0, (long)((char*)&histogram_example.yData), NULL, 0.0, 1, "histogram y and y'?"},
  {"LONGIT_DATA", "", IS_LONG, 0, (long)((char*)&histogram_example.longitData), NULL, 0.0, 1, "histogram t and p?"},
  {"BIN_SIZE_FACTOR", "", IS_DOUBLE, 0, (long)((char*)&histogram_example.binSizeFactor), NULL, 1.0, 0, "multiply computed bin size by this factor before histogramming"},
};
  
CSRDRIFT csrdrift_example;
/* CSR drift length physical parameters */
PARAMETER csrdrift_param[N_CSRDRIFT_PARAMS] = {
    {"L", "M", IS_DOUBLE, 1, (long)((char *)&csrdrift_example.length), NULL, 0.0, 0, "length"},
    {"ATTENUATION_LENGTH", "M", IS_DOUBLE, 0, (long)((char *)&csrdrift_example.attenuationLength), NULL, 0.0, 0, "exponential attenuation length for wake"},
    {"DZ", "", IS_DOUBLE, 0, (long)((char *)&csrdrift_example.dz), NULL, 0.0, 0, "interval between kicks"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&csrdrift_example.nKicks), NULL, 0.0, 1, "number of kicks (if DZ is zero)"},
    {"SPREAD", "", IS_LONG, 0, (long)((char *)&csrdrift_example.spread), NULL, 0.0, 0, "use spreading function?"},
    {"USE_OVERTAKING_LENGTH", "", IS_LONG, 0, (long)((char *)&csrdrift_example.useOvertakingLength), NULL, 0.0, 0, "use overtaking length for ATTENUATION_LENGTH?"},
    {"USE_SALDIN54", "", IS_LONG, 0, (long)((char*)&csrdrift_example.useSaldin54), NULL, 0, 0, "Use Saldin et al eq. 54 (NIM A 398 (1997) 373-394 for decay vs z?"},
    {"SALDIN54POINTS", "", IS_LONG, 0, (long)((char*)&csrdrift_example.nSaldin54Points), NULL, 0.0, 20, "Number of values of position inside bunch to average for Saldin eq 54."},
    {"CSR", "", IS_LONG, 0, (long)((char *)&csrdrift_example.csr), NULL, 0.0, 1, "do CSR calcuations"},
    {"SPREAD_MODE", "", IS_STRING, 0, (long)((char *)&csrdrift_example.spreadMode), "full", 0.0, 0, "full, simple, or radiation-only"},
    {"WAVELENGTH_MODE", "", IS_STRING, 0, (long)((char *)&csrdrift_example.wavelengthMode), "sigmaz", 0.0, 0, "sigmaz or peak-to-peak"},
    {"BUNCHLENGTH_MODE", "", IS_STRING, 0, (long)((char *)&csrdrift_example.bunchlengthMode), "68-percentile", 0.0, 0, "rms, 68-percentile, or 90-percentile"},
    {"SALDIN54_OUTPUT", "", IS_STRING, 0, (long)((char*)&csrdrift_example.Saldin54Output), NULL, 0.0, 0, "Filename for output of CSR intensity vs. z as computed using Saldin eq 54."},
    {"OL_MULTIPLIER", "", IS_DOUBLE, 0, (long)((char *)&csrdrift_example.overtakingLengthMultiplier), NULL, 1.0, 0, "factor by which to multiply the overtaking length to get the attenuation length"},
    };

RFCW rfcw_example;
/* rf cavity with wakes physical parameters */
PARAMETER rfcw_param[N_RFCW_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&rfcw_example.length), NULL, 0.0, 0, "length"},
    {"CELL_LENGTH", "M", IS_DOUBLE, 0, (long)((char *)&rfcw_example.cellLength), NULL, 0.0, 0, "cell length (used to scale wakes)"},
    {"VOLT", "V", IS_DOUBLE, 0, (long)((char *)&rfcw_example.volt), NULL, 0.0, 0, "voltage"},
    {"PHASE", "DEG", IS_DOUBLE, 0, (long)((char *)&rfcw_example.phase), NULL, 0.0, 0, "phase"},
    {"FREQ", "Hz", IS_DOUBLE, 0, (long)((char *)&rfcw_example.freq), NULL, 500.0e6, 0, "frequency"},
    {"Q", "", IS_DOUBLE, 0, (long)((char *)&rfcw_example.Q), NULL, 0.0, 0, "cavity Q"},
    {"PHASE_REFERENCE", "", IS_LONG, 0, (long)((char *)&rfcw_example.phase_reference), NULL, 0.0, 0, "phase reference number (to link with other time-dependent elements)"},
    {"CHANGE_P0", "", IS_LONG, 0, (long)((char *)&rfcw_example.change_p0), NULL, 0.0, 0, "does element change central momentum?"}, 
    {"CHANGE_T", "", IS_LONG, 0, (long)((char *)&rfcw_example.change_t), NULL, 0.0, 0, "not recommended"}, 
    {"FIDUCIAL", "", IS_STRING, 0, (long)((char *)&rfcw_example.fiducial), NULL, 0.0, 0, "mode for determining fiducial arrival time (light, tmean, first, pmaximum)"},
    {"END1_FOCUS", "", IS_LONG, 0, (long)((char *)&rfcw_example.end1Focus), NULL, 0.0, 0, "include focusing at entrance?"},
    {"END2_FOCUS", "", IS_LONG, 0, (long)((char *)&rfcw_example.end2Focus), NULL, 0.0, 0, "include focusing at exit?"},
    {"N_KICKS", "", IS_LONG, 0, (long)((char *)&rfcw_example.nKicks), NULL, 0.0, 1, "number of kicks to use.  Set to zero for matrix method."},
    {"WAKEFILE", "", IS_STRING, 0, (long)((char *)&rfcw_example.wakeFile), NULL, 0.0, 0, "name of file containing Green functions"},
    {"ZWAKEFILE", "", IS_STRING, 0, (long)((char *)&rfcw_example.zWakeFile), NULL, 0.0, 0, "if WAKEFILE=NULL, optional name of file containing longitudinal Green function"},
    {"TRWAKEFILE", "", IS_STRING, 0, (long)((char *)&rfcw_example.trWakeFile), NULL, 0.0, 0, "if WAKEFILE=NULL, optional name of file containing transverse Green functions"},
    {"TCOLUMN", "", IS_STRING, 0, (long)((char *)&rfcw_example.tColumn), NULL, 0.0, 0, "column containing time data"},
    {"WXCOLUMN", "", IS_STRING, 0, (long)((char *)&rfcw_example.WxColumn), NULL, 0.0, 0, "column containing x Green function"},
    {"WYCOLUMN", "", IS_STRING, 0, (long)((char *)&rfcw_example.WyColumn), NULL, 0.0, 0, "column containing y Green function"},
    {"WZCOLUMN", "", IS_STRING, 0, (long)((char *)&rfcw_example.WzColumn), NULL, 0.0, 0, "column containing longitudinal Green function"},
    {"N_BINS", "", IS_LONG, 0, (long)((char *)&rfcw_example.n_bins), NULL, 0.0, 0, "number of bins for current histogram"},
    {"INTERPOLATE", "", IS_LONG, 0, (long)((char *)&rfcw_example.interpolate), NULL, 0.0, 0, "interpolate wake?"},
    {"SMOOTHING", "", IS_LONG, 0, (long)((char *)&rfcw_example.smoothing), NULL, 0.0, 0, "smooth current histogram?"},
    {"SG_HALFWIDTH", "", IS_LONG, 0, (long)((char *)&rfcw_example.SGHalfWidth), NULL, 0.0, 4, "Savitzky-Golay filter half-width for smoothing"},
    {"SG_ORDER", "", IS_LONG, 0, (long)((char *)&rfcw_example.SGOrder), NULL, 0.0, 1, "Savitzky-Golay filter order for smoothing"},
    {"DX", "M", IS_DOUBLE, 1, (long)((char *)&rfcw_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 1, (long)((char *)&rfcw_example.dy), NULL, 0.0, 0, "misalignment"},
    };
   
REMCOR remcor_example;
/* beam centering physical parameters */
PARAMETER remcor_param[N_REMCOR_PARAMS]={
    {"X" , "", IS_LONG, 0, (long)((char *)&remcor_example.x), NULL, 0.0, 1, "remove correlations in x?"},
    {"XP", "", IS_LONG, 0, (long)((char *)&remcor_example.xp), NULL, 0.0, 1, "remove correlations in x'?"},
    {"Y" , "", IS_LONG, 0, (long)((char *)&remcor_example.y), NULL, 0.0, 1, "remove correlations in y?"},
    {"YP", "", IS_LONG, 0, (long)((char *)&remcor_example.yp), NULL, 0.0, 1, "remove correlations in y'?"},
    {"WITH", "", IS_LONG, 0, (long)((char *)&remcor_example.with), NULL, 0.0, 6, "coordinate to remove correlations with (1,2,3,4,5,6)=(x,x',y,y',s,dP/Po)"},
    {"ONCE_ONLY", "", IS_LONG, 0, (long)((char *)&remcor_example.onceOnly), NULL, 0.0, 0, "compute correction only for first beam, apply to all?"},
    };

MAP_SOLENOID mapSol_example;

PARAMETER mapSolenoid_param[N_MAPSOLENOID_PARAMS] = {
    {"L", "M", IS_DOUBLE, 0, (long)((char *)&mapSol_example.length), NULL, 0.0, 0, "length"},
    {"DX", "M", IS_DOUBLE, 0, (long)((char *)&mapSol_example.dx), NULL, 0.0, 0, "misalignment"},
    {"DY", "M", IS_DOUBLE, 0, (long)((char *)&mapSol_example.dy), NULL, 0.0, 0, "misalignment"},
    {"N_STEPS", "", IS_LONG, 0, (long)((char *)&mapSol_example.n_steps), NULL, 0.0, 100, "number of steps (for nonadaptive integration)"},
    {"INPUTFILE", "", IS_STRING, 0, (long)((char *)&mapSol_example.inputFile), NULL, 0.0, 0, "SDDS file containing (Br, Bz) vs (r, z).  Each page should have values for a fixed r."},
    {"RCOLUMN", "", IS_STRING, 0, (long)((char *)&mapSol_example.rColumn), NULL, 0.0, 0, "column containing r values"},
    {"ZCOLUMN", "", IS_STRING, 0, (long)((char *)&mapSol_example.zColumn), NULL, 0.0, 0, "column containing z values"},
    {"BRCOLUMN", "", IS_STRING, 0, (long)((char *)&mapSol_example.BrColumn), NULL, 0.0, 0, "column containing Br values"},
    {"BZCOLUMN", "", IS_STRING, 0, (long)((char *)&mapSol_example.BzColumn), NULL, 0.0, 0, "column containing Bz values"},
    {"FACTOR", "", IS_DOUBLE, 0, (long)((char *)&mapSol_example.factor), NULL, DEFAULT_ACCURACY, 0, "factor by which to multiply fields"},
    {"ACCURACY", "", IS_DOUBLE, 0, (long)((char *)&mapSol_example.accuracy), NULL, DEFAULT_ACCURACY, 0, "integration accuracy"},
    {"METHOD", " ", IS_STRING, 0, (long)((char *)&mapSol_example.method), DEFAULT_INTEG_METHOD, 0.0, 0, "integration method (runge-kutta, bulirsch-stoer, non-adaptive runge-kutta, modified midpoint)"},
    } ;

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
    {    N_RFCA_PARAMS,     MAT_LEN_NCAT|HAS_RF_MATRIX,       sizeof(RFCA),    rfca_param     }, 
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
    {    N_RFTMEZ0_PARAMS,  MAT_LEN_NCAT,    sizeof(RFTMEZ0),    rftmez0_param     }, 
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
    {  N_RAMPRF_PARAMS, MAT_LEN_NCAT|HAS_RF_MATRIX,    sizeof(RAMPRF),    ramprf_param   },
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
    { N_MODRF_PARAMS, MAT_LEN_NCAT|HAS_RF_MATRIX,       sizeof(MODRF),    modrf_param     }, 
    { N_BMAPXY_PARAMS,     HAS_LENGTH,   sizeof(BMAPXY),  bmapxy_param      },
    { N_ZTRANSVERSE_PARAMS,          0,     sizeof(ZTRANSVERSE),  ztransverse_param    },
    { N_IBSCATTER_PARAMS,        0,   sizeof(IBSCATTER),  ibscatter_param  },
    { N_FMULT_PARAMS,  MAT_LEN_NCAT|IS_MAGNET,       sizeof(FMULT),    fmult_param     }, 
    { N_WAKE_PARAMS, 0, sizeof(WAKE), wake_param},
    { N_TRWAKE_PARAMS, 0, sizeof(TRWAKE), trwake_param},
    { N_TUBEND_PARAMS, 0, sizeof(TUBEND), tubend_param},
    { N_CHARGE_PARAMS, 0, sizeof(CHARGE), charge_param},
    { N_PFILTER_PARAMS, 0, sizeof(PFILTER), pfilter_param},
    { N_HISTOGRAM_PARAMS, 0, sizeof(HISTOGRAM), histogram_param},
    {  N_CSRCSBEND_PARAMS, MAT_LEN_NCAT|IS_MAGNET,
                                         sizeof(CSRCSBEND),    csrcsbend_param   },
    {  N_CSRDRIFT_PARAMS, MAT_LEN_NCAT,
                                         sizeof(CSRDRIFT),    csrdrift_param   },
    {    N_RFCW_PARAMS,     MAT_LEN_NCAT|HAS_RF_MATRIX,       sizeof(RFCW),    rfcw_param     }, 
    {  N_REMCOR_PARAMS,           0,     sizeof(REMCOR),    remcor_param   },
    { N_MAPSOLENOID_PARAMS,  MAT_LEN_NCAT,    sizeof(MAP_SOLENOID),    mapSolenoid_param    }, 
} ;
 

void compute_offsets()
{
    long i, j;
    for (i=0; i<N_TYPES; i++) {
        for (j=entity_description[i].n_params-1; j>=0; j--)
            entity_description[i].parameter[j].offset -= entity_description[i].parameter[0].offset;
        }
    }

