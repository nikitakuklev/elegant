set turns 400
set xchrom 1
set ychrom 1
set queue all.q
set cores 8
set inputList [list matchTemplate.ele optim1Template1.ele optim1Template2.ele optim1Template3.ele optim1Template4.ele]
set modeList [list p s p p p]

# Ring parameters
set sectors 40

# Emittance constraint 
set ex0Target 80e-12
set exEffTarget 80e-12
set exTargetQuantity ex0

# Error levels
set fseError 2e-4
set tiltError 10e-4
set dxError 10e-6
set dyError 10e-6

# DA and MMAP parameters
set rfHarmonic 1296
set xmax 0.015
set ymax 0.0015
set s_end [expr 27.6*6]
set name_pattern S*:S?

# Script for processing jobs, plus related parameters
set processScript processJob1
# Largest momentum deviation (artificially imposed)
set deltaLimit 6e-2
# Beam parameters
set current 100.0
set bunches 24
set coupling 0.01
set sigmaz 10
set energy 7e3
# Set to large value (e.g., 1.0) to disable
set clipPositive 0.007
