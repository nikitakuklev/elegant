;;; elegant.el --- Major mode for editing Elegant files in Emacs

;; Copyright (C) 2021 Oscar BLANCO

;; Author     : orblancog <orblancog@gmail.com>
;; Maintainer : orblancog
;; Created    : 01 Nov 2021
;; Keywords   : languages
;; Homepage   : https://github.com/orblancog/
;; Version    : 0.1

;; This file is not part of GNU Emacs

;; This program is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.

;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <https://www.gnu.org/licenses/>.

;;; Commentary:
;; ** FEATURES in version 0.1 **
;;   * Highlights commands, parameters and special operators in Elegant
;;   * If the file extension is '.ele' or '.lte' then the buffer is automatically highlighted,
;;     but, any buffer can be highlighted by typing the following sequence :
;;        M-x elegant-mode RETURN
;;     where, M+x means to type the META character (i.e. ALT in most systems) and
;;     the x key at the same time, and then type elegant-mode and the RETURN key.
;;   * If the limit of 80 characters is desired, they could be highlighted differently
;;     by uncomment the line
;;       ;(require 'whitespace)
;;       ;(setq whitespace-line-column 80) ;; limit line length
;;       ;(setq whitespace-style '(face lines-tail))
;;       ;(add-hook 'elegant-mode-hook 'whitespace-mode)
;;     in this (elegant.el) file.
;; ** HOW TO INSTALL **
;;  1) Download the file elegant.el from https://github.com/orblancog/mad-x_syntax
;;  2)
;;     - Emacs version >= 25. Open emacs Use the command
;;          M+X package-install-file RETURN elegant.el RETURN
;;       And go to 3)
;;     - Emacs version previous to 25
;;         1. Copy this file (elegant.el) to the highlighting definition folder,
;;            the exact location depends on the emacs version
;;              a)  ~/.emacs.d/lisp/  ---> (Emacs v24.5.X or v25)
;;              b)  ~/.emacs.d/       ---> (Emacs v23.X.X)
;;              c)  ~/                ---> (Emacs v21.X.X)
;;            e.g.
;;               $ cp elegant.el ~/.emacs.d/lisp/
;;         2. Edit or create your .emacs file,
;;            the exact location depends on the emacs version,
;;            typically in ~/ or in ~/.emacs.d/
;;              Add the following block where the load-path must match point 1.
;;
;;              ;;;; START OF BLOCK TO COPY AND UNCOMMENT
;;              (global-font-lock-mode t);; Enable syntax highlighting
;;              (setq font-lock-maximum-decoration t)
;;              (add-to-list 'load-path "~/.emacs.d");; <---!!! edit according to point 1!!!
;;              (autoload 'elegant-mode "elegant" "ELEGANT-mode" t)
;;              (setq auto-mode-alist (append '(("\\.\\(\\(elegant\\)\\|\\(seq\\)\\)$" . elegant-mode))
;;                auto-mode-alist))
;;              ;;;; END OF BLOCK TO COPY AND UNCOMMENT
;;
;;  3) Close and reopen emacs to reload the environment variables.
;; **INFO**
;;   * For mad instructions, visit
;;       http://elegant.web.cern.ch/elegant
;;   * Other syntax highlightings could be found inside the
;;     elegant sources.  Check the 'syntax' folder in the elegant dir !
;;        for example: VIM, ATOM, KATE+KWRITE
;;   * Write me to the email address above about any bug including an example.

;;; History:
;; v 0.1 2021nov01   First release at SOLEIL. File is also available in the
;;                   Elegant sources "syntax" folder.

;;; Code:
(defgroup elegant nil
  "Major mode to edit Elegant files in emacs."
  :group 'languages)

(defvar elegant-mode-hook nil
  "Hook for elegant-model initialize nil.")

;;;; add  80 characters line
;;;; (global-whitespace-mode +1)
;(require 'whitespace)
;(setq whitespace-line-column 80) ;; limit line length
;(setq whitespace-style '(face lines-tail))
;(add-hook 'elegant-mode-hook 'whitespace-mode)

(defconst elegant-font-lock-keywords-face-all
  ;; elegant-font-lock-keywords-programflow
  `((,(regexp-opt '(;; font-lock-keyword-face
		    ;; elegant-font-lock-keywords-programflow
		   "STRING"
		   "double"
		   "long"
		   )
		  'words)
     .  font-lock-keyword-face))
  "Highlighting expressions for MAD-X mode (keywords-all).")

(defconst elegant-font-lock-builtin-face-all
  `((,(regexp-opt '(;; font-lock-builtin-face
		   ;; elegant-font-lock-keywords-commands
		   "&alter_elements"
		   "&amplification_factors"
		   "&analyze_map"
		   "&aperture_data"
		   "&aperture_input"
		   "&awe_beam"
		   "&bunched_beam"
		   "&bunched_beam_moments"
		   "&change_end"
		   "&change_particle"
		   "&change_start"
		   "&chaos_map"
		   "&chromaticity"
		   "&closed_orbit"
		   "&correct"
		   "&correct_tunes"
		   "&correction_matrix_output"
		   "&coupled_twiss_output"
		   "&divide_elements"
		   "&elastic_scattering"
		   "&end"
		   "&error"
		   "&error_control"
		   "&error_element"
		   "&find_aperture"
		   "&fit_traces"
		   "&floor_coordinates"
		   "&frequency_map"
		   "&global_settings"
		   "&ignore_elements"
		   "&include_commands"
		   "&inelastic_scattering"
		   "&insert_elements"
		   "&insert_sceffects"
		   "&ion_effects"
		   "&linear_chromatic_tracking_setup"
		   "&link_control"
		   "&link_elements"
		   "&load_parameters"
		   "&matrix_output"
		   "&modulate_elements"
		   "&moments_output"
		   "&momentum_aperture"
		   "&obstruction_data"
		   "&optimization_constraint"
		   "&optimization_covariable"
		   "&optimization_setup"
		   "&optimization_term"
		   "&optimization_variable"
		   "&optimize"
		   "&parallel_optimization_setup"
		   "&print_dictionary"
		   "&ramp_elements"
		   "&replace_elements"
		   "&rf_setup"
		   "&rpn_expression"
		   "&rpn_load"
		   "&run_control"
		   "&run_setup"
		   "&sasefel"
		   "&save_lattice"
		   "&sdds_beam"
		   "&semaphores"
		   "&set_reference_particle_output"
		   "&setup_linear_chromatic_tracking"
		   "&slice_analysis"
		   "&steering_element"
		   "&stop"
		   "&subprocess"
		   "&touschek_scatter"
		   "&trace"
		   "&track"
		   "&transmute_elements"
		   "&tune_footprint"
		   "&tune_shift_with_amplitude"
		   "&twiss_analysis"
		   "&twiss_output"
		   "&vary_element"
		   )
		  'words)
     . font-lock-builtin-face))
  "Highlighting expressions for MAD-X mode (builtin-all).")

(defconst elegant-font-lock-type-face-all
  ;; elegant-font-lock-keywords-elements
  `((,(regexp-opt '(; font-lock-type-face
		    ;; elegant-font-lock-keywords-elements
		   "alph"
		   "apcontour"
		   "beambeam"
		   "bggexp"
		   "bmapxy"
		   "bmxyz"
		   "boffaxe"
		   "branch"
		   "brat"
		   "bumper"
		   "ccbend"
		   "center"
		   "cepl"
		   "charge"
		   "ckicker"
		   "clean"
		   "corgpipe"
		   "cpickup"
		   "csbend"
		   "csrcsbend"
		   "csrdrift"
		   "cwiggler"
		   "drif"
		   "dscatter"
		   "ecol"
		   "edrift"
		   "ehkick"
		   "ekicker"
		   "else"
		   "ematrix"
		   "emittance"
		   "energy"
		   "evkick"
		   "floor"
		   "fmult"
		   "frfmode"
		   "ftable"
		   "ftrfmode"
		   "gfwiggler"
		   "gkickmap"
		   "histogram"
		   "hkick"
		   "hkpoly"
		   "hmon"
		   "ibscatter"
		   "ilmatrix"
		   "ioneffects"
		   "kicker"
		   "koct"
		   "kpoly"
		   "kquad"
		   "kquse"
		   "ksbend"
		   "ksext"
		   "lgbend"
		   "line"
		   "lmirror"
		   "lrwake"
		   "lscdrift"
		   "lsrmdltr"
		   "lthinlens"
		   "magnify"
		   "malign"
		   "mapsolenoid"
		   "mark"
		   "matr"
		   "matter"
		   "maxamp"
		   "mbumper"
		   "mhistogram"
		   "modrf"
		   "moni"
		   "mrfdf"
		   "mult"
		   "nibend"
		   "nisept"
		   "octu"
		   "peppot"
		   "pfilter"
		   "polyseries"
		   "quad"
		   "qufringe"
		   "rampp"
		   "ramprf"
		   "rben"
		   "rcol"
		   "recirc"
		   "reflect"
		   "remcor"
		   "rfca"
		   "rfcw"
		   "rfdf"
		   "rfmode"
		   "rftm110"
		   "rftmez0"
		   "rimult"
		   "rmdf"
		   "rotate"
		   "sample"
		   "sben"
		   "scatter"
		   "scmult"
		   "scraper"
		   "script"
		   "sext"
		   "shrfdf"
		   "slice"
		   "sole"
		   "speedbump"
		   "sreffects"
		   "stray"
		   "taperapc"
		   "taperape"
		   "taperapr"
		   "tfbdriver"
		   "tfbpickup"
		   "tmcf"
		   "trcount"
		   "trfmode"
		   "trwake"
		   "tscatter"
		   "tubend"
		   "twiss"
		   "twla"
		   "twmta"
		   "twpl"
		   "ukickmap"
		   "vkick"
		   "vmon"
		   "wake"
		   "watch"
		   "wiggler"
		   "zlongit"
		   "ztransverse"
		   )
		  'words)
     . font-lock-type-face))
  "Highlighting expressions for MAD-X mode (type-all).")

;; (defconst elegant-font-lock-warning-face-all
;;   ;; elegant-font-lock-keywords-errordef
;;   `((,(regexp-opt '(; font-lock-warning-face
;; 		   ;; elegant-font-lock-keywords-errordef
;; 		   "ERROR"
;; 		   "EALIGN"
;; 		   "EFCOMP"
;; 		   "SETERR"
;; 		   )
;; 		  'words)
;;      . font-lock-warning-face))
;;   "Highlighting expressions for MAD-X mode (warning-all).")

;; (defconst elegant-font-lock-special_operators
;;   ;; elegant-font-lock-special_operators
;;   `((,(regexp-opt '(; font-lock-warning-face
;; 		   ;; elegant-font-lock-special_operators
;; 		   ":="
;; 		   "->"
;; 		   )
;; 		  t)
;;      . font-lock-warning-face))
;;   "Highlighting expressions for MAD-X mode (special-operators).")

(defconst elegant-font-lock-constant-face-all
  ;; elegant-font-lock-keywords-constants
  `((,(regexp-opt '(; font-lock-constant-face
		   ;; elegant-font-lock-keywords-constants
		   "POSITRON"
		   "ELECTRON"
		   )
		  'words)
     . font-lock-constant-face))
  "Highlighting expressions for MAD-X mode (constant-all).")

(defconst elegant-font-lock-special_constants
  ;; elegant-font-lock-special_constants
  `((,(regexp-opt '(; font-lock-constant-face
		   ;; elegant-font-lock-special_constants
		   "NULL"
		   )
		  t)
     . font-lock-constant-face))
  "Highlighting expressions for MAD-X mode (special_constants).")

(defconst elegant-font-lock-doc-face-all
  ;; elegant-font-lock-keywords-parameters
  `((,(regexp-opt '(;; font-lock-doc-face
		    ;; elegant-font-lock-keywords-parameters
;; 		   "NOEXPR" ;2016.08
;; 		   "KNLL" ;2016.08
;; 		   "CNLL" ;2016.08
;; 		   "ROOTMACRO"
;; 		   "MOMENT_S"
;; 		   "MOMENT"
;; 		   "SURVEYALL"
;; 		   "DELTA_DEPENDENCY"
;; 		   "DELTAP_DEPENDENCY"
;; 		   "IGNORE_MAP_ORBIT"
;; 		   "RING_PARAMETERS"
;; 		   "CENTER_MAGNETS"
;; 		   "SUMMARY_FILE"
;; 		   "SUMMARY_TABLE"
;; 		   "MAPTABLE"
;; 		   "DELS"
;; 		   "TXI"
;; 		   "TYI"
;; 		   "TLI"
;; 		   "UPDATE"
;; 		   "NCORR"
;; 		   "SNGVAL"
;; 		   "SNGCUT"
;; 		   "CORRLIM"
;; 		   "TWORING"
;; 		   "UNITS"
;; 		   "CORZERO"
;; 		   "BEAM1TAB"
;; 		   "BEAM2TAB"
;; 		   "EXTERN"
;; 		   "NAME_COL"
;; 		   "X_COL"
;; 		   "Y_COL"
;; 		   "TWISSUM"
;; 		   "LSQ"
;; 		   "MICADO"
;; 		   "SVD"
;; 		   "THICK"
;; 		   "VECTOR"
;; 		   "expr"
;; 		   "VLENGTH"
;; 		   "SLOW"
;; 		   "PATCH_ANG"
;; 		   "PATCH_TRANS"
;; 		   "ADD_ANGLE"
;; 		   "ADD_PASS"
;; 		   "NEXT_SEQU"
;; 		   "ZERO_SUPPR"
;; 		   "N_BESSEL"
;; 		   "RIPKEN"
;; 		   "NO_CAVITY_TOTALPATH"
;; 		   "PNL"
;; 		   "PSL"
;; 		   "MAKEDIPEDGE"
;; 		   "TRUNCATE"
;; 		   "APPEND"
;; 		   "ROW1"
;; 		   "ROW1"
;; 		   "PARAM"
;; 		   "SINKICK"
;; 		   "SINPEAK"
;; 		   "SINTUNE"
;; 		   "SINPHASE"
;; 		   "L"
;; 		   "K1"
;; 		   "K1S"
;; 		   "TILT"
;; 		   "K2"
;; 		   "K2S"
;; 		   "K3"
;; 		   "K3S"
;; 		   "KS"
;; 		   "K3S"
;; 		   "KSI"
;; 		   "VOLT"
;; 		   "LAG"
;; 		   "FREQ"
;; 		   "HARMON"
;; 		   "rv1"
;; 		   "rv2"
;; 		   "rv3"
;; 		   "rv4"
;; 		   "rph1"
;; 		   "rph2"
;; 		   "lagf"
;; 		   "PARTICLE"
;; 		   "MASS"
;; 		   "CHARGE"
;; 		   "ENERGY"
;; 		   "PC"
;; 		   "GAMMA"
;; 		   "BETA"
;; 		   "BRHO"
;; 		   "EX"
;; 		   "EXN"
;; 		   "EY"
;; 		   "EYN"
;; 		   "ET"
;; 		   "SIGT"
;; 		   "SIGE"
;; 		   "KBUNCH"
;; 		   "NPART"
;; 		   "PLANE"
;; 		   "BCURRENT"
;; 		   "CURRENT"
;; 		   "BUNCHED"
;; 		   "RADIATE"
;; 		   "NORM_NO"
;; 		   "BV"
;; 		   "SEQUENCE"
;; 		   "ENDSEQUENCE"
;; 		   "REFER"
;; 		   "MAD8"
;; 		   "GNFU"
;; 		   "KICK"
;; 		   "HKICK"
;; 		   "VKICK"
;; 		   "ANGLE"
;; 		   "K0"
;; 		   "K0S"
;; 		   "RESPLIT"
;; 		   "DAMP"
;; 		   "QUANTUM"
;; 		   "RECLOSS"
;; 		   "ELEMENT_BY_ELEMENT"
;; 		   "NORM"
;; 		   "NORM_OUT"
;; 		   "RADIATION_MODEL1"
;; 		   "RADIATION_ENERGY_LOSS"
;; 		   "RADIATION_QUADr"
;; 		   "BEAM_ENVELOPE"
;; 		   "SPACE_CHARGE"
;; 		   "FX"
;; 		   "FY"
;; 		   "FT"
;; 		   "FFILE"
;; 		   "E1"
;; 		   "E2"
;; 		   "FINT"
;; 		   "FINTX"
;; 		   "HGAP"
;; 		   "H1"
;; 		   "H2"
;; 		   "FLAG"
;; 		   "RANGE"
;; 		   "PLACE"
;; 		   "PARENT"
;; 		   "KEYWORD"
;; 		   "FROM"
;; 		   "AT"
;; 		   "REFPOS"
;; 		   "LENGTH"
;; 		   "EXACT_MIS"
;; 		   "CLASS"
;; 		   "PATTERN"
;; 		   "FILE"
;; 		   "FORMAT"
;; 		   "RBARC"
;; 		   "TWISS_PRINT"
;; 		   "THREADER"
;; 		   "THIN_FOC"
;; 		   "NO_FATAL_STOP"
;; 		   "TEXT"
;; 		   "BARE"
;; 		   "SLICE"
;; 		   "THICK"
;; 		   "COMMAND"
;; 		   "TO"
;; 		   "NEWNAME"
;; 		   "BY"
;; 		   "SELECTED"
;; 		   "H"
;; 		   "LRAD"
;; 		   "KNL"
;; 		   "KSL"
;; 		   "TYPE"
;; 		   "SIGX"
;; 		   "SIGY"
;; 		   "XMA"
;; 		   "YMA"
;; 		   "BBSHAPE"
;; 		   "WIDTH"
;; 		   "BBDIR"
;; 		   "ECHO"
;; 		   "ECHOMACRO"
;; 		   "TRACE"
;; 		   "VERIFY"
;; 		   "PERIOD"
;; 		   "SECTORMAP"
;; 		   "SECTORFILE"
;; 		   "KEEPORBIT"
;; 		   "USEORBIT"
;; 		   "COUPLE"
;; 		   "FULL"
;; 		   "COLUMN"
;; 		   "MOMENTS"
;; 		   "PARAMETRIC"
;; 		   "CLEAR"
;; 		   "POS"
;; 		   "POLYNOMIAL"
;; 		   "MONOMIAL"
;; 		   "PARAMETRIC"
;; 		   "QUANTITY"
;; 		   "ROW"
;; 		   "SEED"
;; 		   "ADD"
;; 		   "INFO"
;; 		   "DEBUG"
;; 		   "VERBOSE"
;; 		   "TELL"
;; 		   "RESET"
;; 		   "WARN"
;; 		   "LABEL"
;; 		   "APERTYPE"
;; 		   "HAXIS"
;; 		   "HMIN"
;; 		   "HMAX"
;; 		   "VAXIS"
;; 		   "VAXIS1"
;; 		   "VAXIS2"
;; 		   "VAXIS3"
;; 		   "VAXIS4"
;; 		   "VMIN"
;; 		   "VMAX"
;; 		   "BARS"
;; 		   "SYMBOL"
;; 		   "NOVERSION"
;; 		   "NO_FATAL_ERROR"
;; 		   "NO_FATAL_STOP"
;; 		   "INTERPOLATE"
;; 		   "NOLINE"
;; 		   "NOTITLE"
;; 		   "MARKER_PLOT"
;; 		   "RANGE_PLOT"
;; 		   "MULTIPLE"
;; 		   "PTC"
;; 		   "PTC_TABLE"
;; 		   "TRACKFILE"
;; 		   "CAVALL"
;; 		   "MULT_AUTO_OFF"
;; 		   "MAX_MULT_ORD"
;; 		   "SPLIT"
;; 		   "RADIUS"
;; 		   "WARNING"
;; 		   "STYLE"
;; 		   "COLOUR"
;; 		   "TURNS"
;; 		   "EVERYSTEP"
;; 		   "ONETABLE"
;; 		   "TABLEALLSTEPS"
;; 		   "GCS"
;; 		   "ROOTNTUPLE"
;; 		   "EXTENSION"
;; 		   "FASTUNE"
;; 		   "MAXAPER"
;; 		   "LYAPUNOV"
;; 		   "ORBIT"
;; 		   "TOL"
;; 		   "DS"
;; 		   "DPHI"
;; 		   "DTHETA"
;; 		   "DPSI"
;; 		   "MREX"
;; 		   "MREY"
;; 		   "MSCALX"
;; 		   "MSCALY"
;; 		   "AREX"
;; 		   "AREY"
;; 		   "ORDER"
;; 		   "DKN"
;; 		   "DKS"
;; 		   "DKNR"
;; 		   "DKSR"
;; 		   "HYSTER"
;; 		   "HCOEFFN"
;; 		   "HCOEFFS"
;; 		   "BETA0"
;; 		   "RMATRIX"
;; 		   "STEP"
;; 		   "CHROM"
;; 		   "LOWER"
;; 		   "UPPER"
;; 		   "SLOPE"
;; 		   "OPT"
;; 		   "CALLS"
;; 		   "NO"
;; 		   "XDISTR"
;; 		   "YDISTR"
;; 		   "ZDISTR"
;; 		   "TOLERANCE"
;; 		   "STRATEGY"
;; 		   "REPEAT"
;; 		   "BISEC"
;; 		   "COOL"
;; 		   "BALANCE"
;; 		   "RANDOM"
;; 		   "MODEL"
;; 		   "TARGET"
;; 		   "MODE"
;; 		   "MONERROR"
;; 		   "MONON"
;; 		   "MONSCALE"
;; 		   "PLANEX"
;; 		   "COND"
;; 		   "RESOUT"
;; 		   "CLIST"
;; 		   "MLIST"
;; 		   "STATUS"
;; 		   "POST"
;; 		   "FONT"
;; 		   "LWIDTH"
;; 		   "APER_TOL"
;; 		   "APER_OFFSET"
;; 		   "HALOFILE"
;; 		   "PIPEFILE"
;; 		   "DQF"
;; 		   "BETAQFX"
;; 		   "DP"
;; 		   "DPARX"
;; 		   "DPARY"
;; 		   "COR"
;; 		   "BBEAT"
;; 		   "NCO"
;; 		   "HALO"
;; 		   "INTERVAL"
;; 		   "SPEC"
;; 		   "NOTSIMPLE"
;; 		   "TRUEPROFILE"
;; 		   "OFFSETELEM"
;; 		   "XSIZE"
;; 		   "YSIZE"
;; 		   "ASCALE"
;; 		   "LSCALE"
;; 		   "SSCALE"
;; 		   "RSCALE"
;; 		   "DETUNE"
;; 		   "DISTORT1"
;; 		   "DISTORT2"
;; 		   "START_STOP"
;; 		   "MULTIPOLE_ORDER_RANGE"
;; 		   "NOPRINT"
;; 		   "PRINT_ALL"
;; 		   "PRINT_AT_END"
;; 		   "NOSIXTRACK"
;; 		   "X0"
;; 		   "Y0"
;; 		   "Z0"
;; 		   "THETA0"
;; 		   "PHI0"
;; 		   "PSI0"
;; 		   "SUMM"
;; 		   "CENTRE"
;; 		   "SECTOR_NMUL_MAX"
;; 		   "SECTOR_nMUL"
;; 		   "NTPSA"
;; 		   "SYMPRINT"
;; 		   "TIME"
;; 		   "METHOD"
;; 		   "NST"
;; 		   "EXACT"
;; 		   "OFFSET_DELTAP"
;; 		   "ERRORS_OUT"
;; 		   "ERRORS_IN"
;; 		   "MAGNET_NAME"
;; 		   "RESPLIT"
;; 		   "THIN"
;; 		   "XBEND"
;; 		   "EVEN"
;; 		   "OVERWRITE"
;; 		   "INDEX"
;; 		   "ONEPASS"
;; 		   "DUMP"
;; 		   "DEBUGLEVEL"
;; 		   "LEVEL"
;; 		   "BBORBIT"
;; 		   "SYMPL"
;; 		   "MAXACCELERATION"
;; 		   "EXACT_MISS"
;; 		   "TOTALPATH"
;; 		   "RADIATION"
;; 		   "FRINGE"
;; 		   "ICASE"
;; 		   "CLOSED_ORBIT"
;; 		   "SLICE_MAGNETS"
;; 		   "INITIAL_MATRIX_TABLE"
;; 		   "MATRIX_MANUAL"
;; 		   "INITIAL_MAP_MANUAL"
;; 		   "INITIAL"
;; 		   "ELEMENT"
;; 		   "TRUSTRANGE"
;; 		   "ANHX"
;; 		   "ANHY"
;; 		   "GNUF"
;; 		   "HAML"
;; 		   "EIGN"
;; 		   "INITIAL_MATRIX_MANUAL"
;; 		   "ELEMENTNAME"
;; 		   "KN"
;; 		   "KS"
;; 		   "EXACTMATCH"
;; 		   "ONLYPOSITION"
;; 		   "ONLYORIENTATION"
;; 		   "AUTOPLACEDOWNSTREAM"
;; 		   "REFFRAME"
;; 		   "USE_PTCKNOBS"
		   )
		  'words)
     . font-lock-doc-face))
  "Highlighting expressions for MAD-X mode (doc-all).")

;; (defconst elegant-font-lock-function-name-face-all
;;   ;; elegant-font-lock-keywords-functions
;;   `((,(regexp-opt '(;;  font-lock-function-name-face
;; 		   ;; elegant-font-lock-keywords-functions
;; 		   "SQRT"
;; 		   "LOG"
;; 		   "LOG10"
;; 		   "EXP"
;; 		   "SIN"
;; 		   "COS"
;; 		   "TAN"
;; 		   "ASIN"
;; 		   "ACOS"
;; 		   "ATAN"
;; 		   "SINH"
;; 		   "COSH"
;; 		   "TANH"
;; 		   "SINC"
;; 		   "ABS"
;; 		   "ERF"
;; 		   "ERFC"
;; 		   "FLOOR"
;; 		   "CEIL"
;; 		   "ROUND"
;; 		   "RANF"
;; 		   "GAUSS"
;; 		   "TGAUSS"
;; 		   "FLAT5"
;; 		   "FLAT56"
;; 		   )
;; 		  'words)
;;      . font-lock-function-name-face))
;;   "Highlighting expressions for MAD-X mode (name-all)." )

;; (defconst elegant-font-lock-variable-name-face-all
;;   ;; elegant-font-lock-keywords-variables_elegant
;;   `((,(concat (regexp-opt '(;; font-lock-variable-name-face
;; 			    ;; elegant-font-lock-keywords-variables_elegant
;; 			    "mvar1"
;; 			    "mvar2"
;; 			    "mvar3"
;; 			    "mvar4"
;; 			    "CIRC"
;; 			    "FREQ0"
;; 			    "DTBYDS"
;; 			    "U0"
;; 			    "QS"
;; 			    "ARAD"
;; 			    "PDAMP"
;; 			    "N1MIN"
;; 			    "Z"
;; 			    "PHI"
;; 			    "PSI"
;; 			    "X"
;; 			    "Y"
;; 			    "BETX"
;; 			    "BETY"
;; 			    "NAME"
;; 			    "S"
;; 			    "k0l"
;; 			    "k1l"
;; 			    "k2l"
;; 			    "k3l"
;; 			    "k4l"
;; 			    "K1"
;; 			    "K2"
;; 			    "K3"
;; 			    "K4"
;; 			    "K5"
;; 			    "K6"
;; 			    "K1L"
;; 			    "K2L"
;; 			    "K3L"
;; 			    "K4L"
;; 			    "K5L"
;; 			    "K6L"       
;; 			    "KICK1"
;; 			    "KICK2"
;; 			    "KICK3"
;; 			    "KICK4"
;; 			    "KICK5"
;; 			    "KICK6"
;; 			    "MU1"
;; 			    "MU2"
;; 			    "MU3"
;; 			    "MUX"
;; 			    "MUY"
;; 			    "PX"
;; 			    "PY"
;; 			    "PT"
;; 			    "DELTAP"
;; 			    "XN"
;; 			    "PXN"
;; 			    "WX"
;; 			    "PHI"
;; 			    "THETA"
;; 			    "PHIX"
;; 			    "YN"
;; 			    "PYN"
;; 			    "WY"
;; 			    "PHIY"
;; 			    "TN"
;; 			    "PTN"
;; 			    "WT"
;; 			    "PHIT"
;; 			    "ALFX"
;; 			    "DX"
;; 			    "DPX"
;; 			    "ALFY"
;; 			    "DY"
;; 			    "DPY"
;; 			    "ENERGY"
;; 			    "DMUX"
;; 			    "DDX"
;; 			    "DDPX"
;; 			    "DMUY"
;; 			    "DDY"
;; 			    "DDPY"
;; 			    "Q1"
;; 			    "Q2"
;; 			    "DQ1"
;; 			    "DQ2"
;; 			    "DDQ1"
;; 			    "DDQ2"
;; 			    "N1"
;; 			    "N1X_M"
;; 			    "N1Y_M"
;; 			    "APER_1"
;; 			    "APER_2"
;; 			    "APER_3"
;; 			    "APER_4"
;; 			    "RTOL"
;; 			    "XTOL"
;; 			    "YTOL"
;; 			    "ON_AP"
;; 			    "ON_ELEM"
;; 			    "LENGTH"
;; 			    "ORBIT5"
;; 			    "BETXMAX"
;; 			    "DXMAX"
;; 			    "DXRMS"
;; 			    "XCOMAX"
;; 			    "XRMS"
;; 			    "BETYMAX"
;; 			    "DYMAX"
;; 			    "DYRMS"
;; 			    "YCOMAX"
;; 			    "YCORMS"
;; 			    "SYNCH_1"
;; 			    "SYNCH_2"
;; 			    "SYNCH_3"
;; 			    "SYNCH_4"
;; 			    "SYNCH_5"
;; 			    "DISTANCE"
;; 			    "LYAPUNOV"
;; 			    "LOGDIST"
;; 			    "LOGTURNS"
;; 			    "RE"
;; 					;		    "RE11";REPLACE BY RE[1-6][1-6]
;; 			    "T"
;; 					;		    "T111";; REPLACE BY T[1-6][1-6][1-6]
;; 					;		    "TM111";; REPLACE BY TM[1-6][1-6][1-6]
;; 					;		    "BETA11"; REPLACE BY BETA[1-3][1-3]
;; 					;		    "BETA11P"; REPLACE BY BETA[1-3][1-3]P
;; 			    "ALFA"
;; 					;		    "ALFA11"; REPLACE BY ALFA[1-3][1-3]
;; 					;		    "ALFA11P"; REPLACE BY ALFA[1-3][1-3]P
;; 			    "GAMMATR"
;; 			    "GAMAX"
;; 			    "GAMAY"
;; 					;		    "GAMA11"; REPLACE BY GAMA[1-3][1-3]
;; 					;		    "GAMA11P"; REPLACE BY GAMA[1-3][1-3]P
;; 					;		    "GAMMA11"; REPLACE BY GAMA[1-3][1-3]
;; 					;		    "DISP1; REPLACE BY DISP[1-4]P?[1-3]?"
;; 					;		    "DISP1P1"; REPLACE BY DISP[1-4]P[1-3]
;; 					;		    "EIGN11"; REPLACE BY EIGN[1-6][1-6]
;; 			    "R"
;; 					;		    "R11";; REPLACE BY R[1-6][1-6] AFTER REGEXP-OPT
;; 					;		    "RM11";; REPLACE BY RM[1-6][1-6] AFTER REGEXP-OPT		    
;; 			    )
;; 			  'words)
;; 	      ;; some variables already optimized
;; 	      "\\|\\<RE[1-6][1-6]\\>"
;; 	      "\\|\\<TM?[1-6][1-6][1-6]\\>"
;; 	      "\\|\\<BETA[1-3][1-3]P?\\>"
;; 	      "\\|\\<ALFA[1-3][1-3]P?\\>"
;; 	      "\\|\\<GAMA[1-3][1-3]P?\\>"
;; 	      "\\|\\<DISP[1-4]P?[1-3]?\\>"
;; 	      "\\|\\<EIGN[1-6][1-6]\\>"
;; 	      "\\|\\<RM?[1-6][1-6]\\>")
;;      . font-lock-variable-name-face))
;;   "Highlighting expressions for MAD-X mode (variable-name-all).")

(defconst elegant-font-lock-intfp-name-face-all
  ;; elegant- fonts for integers and floating point numbers
  (list
   '("\\<\\(\\([0-9]+\\.?[0-9]*\\|\\.[0-9]+\\)\\([eE][+-]?\\([0-9]+\\.?[0-9]*\\|[0-9]*\\.[0-9]+\\)\\)?\\)\\>"
     . font-lock-keyword-face))
  "Highlighting expresssions for MAD-X mode (integers and floats).")

(defconst elegant-font-lock-keywords-4
  (append
   elegant-font-lock-special_constants
   ;; elegant-font-lock-special_operators
   elegant-font-lock-keywords-face-all
   elegant-font-lock-constant-face-all
   ;; elegant-font-lock-function-name-face-all
   elegant-font-lock-type-face-all
   ;; elegant-font-lock-variable-name-face-all
   elegant-font-lock-builtin-face-all
   ;; elegant-font-lock-warning-face-all
   ;; elegant-font-lock-doc-face-all
   elegant-font-lock-intfp-name-face-all)
 " highlighting in Elegant mode.")

(defvar elegant-font-lock-keywords elegant-font-lock-keywords-4
  "Default highlighting expressions for MAD-X mode.")

(defvar elegant-mode-syntax-table
  (let ((elegant-mode-syntax-table (make-syntax-table)))
	
  ;; This is added so entity names with underscores and dots can be more easily parsed
  (modify-syntax-entry ?& "w" elegant-mode-syntax-table)
  (modify-syntax-entry ?_ "w" elegant-mode-syntax-table)
  (modify-syntax-entry ?. "w" elegant-mode-syntax-table)
	
  ;;  Comment styles are similar to C++
  (modify-syntax-entry ?/ ". 124 b" elegant-mode-syntax-table)
  (modify-syntax-entry ?* ". 23" elegant-mode-syntax-table)
  (modify-syntax-entry ?\n "> b" elegant-mode-syntax-table)
  (modify-syntax-entry ?! "< b" elegant-mode-syntax-table)
  (modify-syntax-entry ?' "|" elegant-mode-syntax-table)
     elegant-mode-syntax-table)
  "Syntax table for `elegant-mode'.")

;;;###autoload
(define-derived-mode elegant-mode fundamental-mode "elegant"
  "Major mode for editing Elegant script files."
  (setq font-lock-defaults '(elegant-font-lock-keywords nil t))
  ;; Setting up Imenu
  (setq imenu-generic-expression nil)
  (setq imenu-prev-index-position-function nil)
  (setq imenu-extract-index-name-function nil)
  ;;  (imenu-create-index-function)
  ;; Set up search
  (setq case-fold-search t)
  )
;; Highlighting .elegant and .seq buffers
;;;###autoload
(add-to-list 'auto-mode-alist '("\\.\\(ele\\|lte\\)\\'" . elegant-mode))
(provide 'elegant-mode)
;;; elegant.el ends here


