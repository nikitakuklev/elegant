#*************************************************************************
# Copyright (c) 2002 The University of Chicago, as Operator of Argonne
# National Laboratory.
# Copyright (c) 2002 The Regents of the University of California, as
# Operator of Los Alamos National Laboratory.
# This file is distributed subject to a Software License Agreement found
# in the file LICENSE that is included with this distribution. 
#*************************************************************************


TOP=../..
include $(TOP)/configure/CONFIG
include $(TOP)/src/elegant/Makefile.OAG
include $(TOP)/configure/RULES

alter$(OBJ): alter.h

alter.h: ../alter.nl
	nlpp -suppressSummaryVariables ../alter.nl alter.h

amplif$(OBJ): amplif.h

amplif.h: ../amplif.nl
	nlpp -suppressSummaryVariables ../amplif.nl amplif.h

analyze$(OBJ): analyze.h

analyze.h: ../analyze.nl
	nlpp -suppressSummaryVariables ../analyze.nl analyze.h

aperture_search$(OBJ): aperture_search.h

aperture_search.h: ../aperture_search.nl
	nlpp -suppressSummaryVariables ../aperture_search.nl aperture_search.h

bunched_beam$(OBJ): bunched_beam.h

bunched_beam.h: ../bunched_beam.nl
	nlpp -suppressSummaryVariables ../bunched_beam.nl bunched_beam.h

chrom$(OBJ): chrom.h

chrom.h: ../chrom.nl
	nlpp -suppressSummaryVariables ../chrom.nl chrom.h

closed_orbit$(OBJ): closed_orbit.h

closed_orbit.h: ../closed_orbit.nl
	nlpp -suppressSummaryVariables ../closed_orbit.nl closed_orbit.h

correct$(OBJ): correct.h steer_elem.h

correct.h: ../correct.nl
	nlpp -suppressSummaryVariables ../correct.nl correct.h

coupled_twiss$(OBJ): coupled_twiss.h

coupled_twiss.h: ../coupled_twiss.nl
	nlpp -suppressSummaryVariables ../coupled_twiss.nl coupled_twiss.h

divideElements$(OBJ): divideElements.h

divideElements.h: ../divideElements.nl
	nlpp -suppressSummaryVariables ../divideElements.nl divideElements.h

global_settings$(OBJ): global_settings.h

global_settings.h: ../global_settings.nl
	nlpp -suppressSummaryVariables ../global_settings.nl global_settings.h

transmuteElements$(OBJ): transmuteElements.h

transmuteElements.h: ../transmuteElements.nl
	nlpp -suppressSummaryVariables ../transmuteElements.nl transmuteElements.h

ignoreElements$(OBJ): ignoreElements.h

ignoreElements.h: ../ignoreElements.nl
	nlpp -suppressSummaryVariables ../ignoreElements.nl ignoreElements.h

elegant$(OBJ): elegant.h aperture_data.h

elegant.h: ../elegant.nl
	nlpp ../elegant.nl elegant.h

aperture_data.h: ../aperture_data.nl
	nlpp -suppressSummaryVariables ../aperture_data.nl aperture_data.h

error$(OBJ): error.h

error.h: ../error.nl
	nlpp -suppressSummaryVariables ../error.nl error.h

floor$(OBJ): floor.h

floor.h: ../floor.nl
	nlpp -suppressSummaryVariables ../floor.nl floor.h

frequencyMap$(OBJ): frequencyMap.h

frequencyMap.h: ../frequencyMap.nl
	nlpp -suppressSummaryVariables ../frequencyMap.nl frequencyMap.h

chaosMap$(OBJ): chaosMap.h

chaosMap.h: ../chaosMap.nl
	nlpp -suppressSummaryVariables ../chaosMap.nl chaosMap.h

insertSCeffects$(OBJ): insertSCeffects.h

insertSCeffects.h: ../insertSCeffects.nl
	nlpp -suppressSummaryVariables ../insertSCeffects.nl insertSCeffects.h

ionEffects$(OBJ): ionEffects.h

ionEffects.h: ../ionEffects.nl
	nlpp -suppressSummaryVariables ../ionEffects.nl ionEffects.h

insert_elements$(OBJ): insert_elements.h

insert_elements.h: ../insert_elements.nl
	nlpp -suppressSummaryVariables ../insert_elements.nl insert_elements.h

replace_elements$(OBJ): replace_elements.h

replace_elements.h: ../replace_elements.nl
	nlpp -suppressSummaryVariables ../replace_elements.nl replace_elements.h

touschekScatter$(OBJ): touschekScatter.h

touschekScatter.h: ../touschekScatter.nl
	nlpp -suppressSummaryVariables ../touschekScatter.nl touschekScatter.h

link_elements$(OBJ): link_elements.h

link_elements.h: ../link_elements.nl
	nlpp -suppressSummaryVariables ../link_elements.nl link_elements.h

load_parameters$(OBJ): load_parameters.h

load_parameters.h: ../load_parameters.nl
	nlpp -suppressSummaryVariables ../load_parameters.nl load_parameters.h

matrix_output$(OBJ): matrix_output.h

matrix_output.h: ../matrix_output.nl
	nlpp -suppressSummaryVariables ../matrix_output.nl matrix_output.h

modulate$(OBJ): modulate.h

modulate.h: ../modulate.nl
	nlpp -suppressSummaryVariables ../modulate.nl modulate.h

momentumAperture$(OBJ): momentumAperture.h

momentumAperture.h: ../momentumAperture.nl
	nlpp -suppressSummaryVariables ../momentumAperture.nl momentumAperture.h

moments$(OBJ): moments.h

moments.h: ../moments.nl
	nlpp -suppressSummaryVariables ../moments.nl moments.h

optim_covariable$(OBJ): optim_covariable.h

optim_covariable.h: ../optim_covariable.nl
	nlpp -suppressSummaryVariables ../optim_covariable.nl optim_covariable.h

optimize$(OBJ): optimize.h optim_covariable.h

optimize.h: ../optimize.nl
	nlpp -suppressSummaryVariables ../optimize.nl optimize.h

ramp$(OBJ): ramp.h

ramp.h: ../ramp.nl
	nlpp -suppressSummaryVariables ../ramp.nl ramp.h

response$(OBJ): response.h

response.h: ../response.nl
	nlpp -suppressSummaryVariables ../response.nl response.h

run_rpnexpr$(OBJ): run_rpnexpr.h

run_rpnexpr.h: ../run_rpnexpr.nl
	nlpp -suppressSummaryVariables ../run_rpnexpr.nl run_rpnexpr.h

sasefel$(OBJ): sasefel.h

sasefel.h: ../sasefel.nl
	nlpp -suppressSummaryVariables ../sasefel.nl sasefel.h

get_beamline$(OBJ): save_lattice.h

save_lattice$(OBJ): save_lattice.h

save_lattice.h: ../save_lattice.nl
	nlpp -suppressSummaryVariables ../save_lattice.nl save_lattice.h

sdds_beam$(OBJ): sdds_beam.h

sdds_beam.h: ../sdds_beam.nl
	nlpp -suppressSummaryVariables ../sdds_beam.nl sdds_beam.h

sliceAnalysis$(OBJ): sliceAnalysis.h

sliceAnalysis.h: ../sliceAnalysis.nl
	nlpp -suppressSummaryVariables ../sliceAnalysis.nl sliceAnalysis.h

inelasticScattering$(OBJ): inelasticScattering.h

inelasticScattering.h: ../inelasticScattering.nl
	nlpp -suppressSummaryVariables ../inelasticScattering.nl inelasticScattering.h

elasticScattering$(OBJ): elasticScattering.h

elasticScattering.h: ../elasticScattering.nl
	nlpp -suppressSummaryVariables ../elasticScattering.nl elasticScattering.h

steer_elem$(OBJ): steer_elem.h

steer_elem.h: ../steer_elem.nl
	nlpp -suppressSummaryVariables ../steer_elem.nl steer_elem.h

subprocess$(OBJ): subprocess.h

subprocess.h: ../subprocess.nl
	nlpp -suppressSummaryVariables ../subprocess.nl subprocess.h

trace$(OBJ): trace.h

trace.h: ../trace.nl
	nlpp -suppressSummaryVariables ../trace.nl trace.h

tune$(OBJ): tune.h

tune.h: ../tune.nl
	nlpp -suppressSummaryVariables ../tune.nl tune.h

tuneFootprint$(OBJ): tuneFootprint.h

tuneFootprint.h: ../tuneFootprint.nl
	nlpp -suppressSummaryVariables ../tuneFootprint.nl tuneFootprint.h

twiss$(OBJ): twiss.h

twiss.h: ../twiss.nl
	nlpp -suppressSummaryVariables ../twiss.nl twiss.h

vary$(OBJ): vary.h

vary.h: ../vary.nl
	nlpp -suppressSummaryVariables ../vary.nl vary.h

fitTraces$(OBJ): fitTraces.h

fitTraces.h: ../fitTraces.nl
	nlpp -suppressSummaryVariables ../fitTraces.nl fitTraces.h

obstructionData$(OBJ): obstructionData.h

obstructionData.h: ../obstructionData.nl
	nlpp -suppressSummaryVariables ../obstructionData.nl obstructionData.h

elegantLocation = $(wildcard ../../bin/$(EPICS_HOST_ARCH)/elegant elegant)
PelegantLocation = $(wildcard ../../bin/$(EPICS_HOST_ARCH)/Pelegant Pelegant)

elegant = $(words $(notdir $(elegantLocation)))
Pelegant = $(words $(notdir $(PelegantLocation)))


PelegantNewer=0

ifeq ($(EPICS_HOST_ARCH),darwin-x86)
ifeq ($(Pelegant), 1)
PelegantTime=$(shell stat -f %m $(PelegantLocation))
ifeq ($(elegant), 1)
elegantTime=$(shell stat -f %m $(elegantLocation))
PelegantNewer=$(shell rpnl "$(elegantTime) $(PelegantTime) < ? 1 : 0 $$")
else
PelegantNewer=1
endif
endif
else
ifeq ($(Pelegant), 1)
PelegantTime=$(shell stat --format=%Y $(PelegantLocation))
ifeq ($(elegant), 1)
elegantTime=$(shell stat --format=%Y $(elegantLocation))
PelegantNewer=$(shell rpnl "$(elegantTime) $(PelegantTime) < ? 1 : 0 $$")
else
PelegantNewer=1
endif
endif
endif


ifeq ($(PelegantNewer), 0)
Pelegant:
	$(RM) *.o O.$(EPICS_HOST_ARCH)/*.o O.$(EPICS_HOST_ARCH)/*.h
	$(MV) $(TOP)/src/elegant/Makefile $(TOP)/src/elegant/Makefile.TMP
	$(CP) $(TOP)/src/elegant/Makefile.Pelegant $(TOP)/src/elegant/Makefile
	$(MAKE) MPI=1 -f Makefile
	$(RM) O.$(EPICS_HOST_ARCH)/insertSCeffects.o O.$(EPICS_HOST_ARCH)/drand_oag.o O.$(EPICS_HOST_ARCH)/mad_parse.o
	$(RM) $(TOP)/src/elegant/Makefile
	$(MV) $(TOP)/src/elegant/Makefile.TMP $(TOP)/src/elegant/Makefile

elegant:
	$(RM) O.$(EPICS_HOST_ARCH)/link_date.o
	$(MV) $(TOP)/src/elegant/Makefile $(TOP)/src/elegant/Makefile.TMP
	$(CP) $(TOP)/src/elegant/Makefile.Pelegant $(TOP)/src/elegant/Makefile
	$(MAKE) NOMPI=1 -f Makefile
	$(RM) $(TOP)/src/elegant/Makefile
	$(MV) $(TOP)/src/elegant/Makefile.TMP $(TOP)/src/elegant/Makefile

all buildInstall:
	$(RM) O.$(EPICS_HOST_ARCH)/link_date.o
	$(MV) $(TOP)/src/elegant/Makefile $(TOP)/src/elegant/Makefile.TMP
	$(CP) $(TOP)/src/elegant/Makefile.Pelegant $(TOP)/src/elegant/Makefile
	$(MAKE) NOMPI=1 -f Makefile
	$(RM) $(TOP)/src/elegant/Makefile
	$(MV) $(TOP)/src/elegant/Makefile.TMP $(TOP)/src/elegant/Makefile

else
Pelegant:
	$(RM) O.$(EPICS_HOST_ARCH)/link_date.o O.$(EPICS_HOST_ARCH)/insertSCeffects.o O.$(EPICS_HOST_ARCH)/drand_oag.o O.$(EPICS_HOST_ARCH)/mad_parse.o
	$(MV) $(TOP)/src/elegant/Makefile $(TOP)/src/elegant/Makefile.TMP
	$(CP) $(TOP)/src/elegant/Makefile.Pelegant $(TOP)/src/elegant/Makefile
	$(MAKE) MPI=1 -f Makefile
	$(RM) O.$(EPICS_HOST_ARCH)/insertSCeffects.o O.$(EPICS_HOST_ARCH)/drand_oag.o O.$(EPICS_HOST_ARCH)/mad_parse.o
	$(RM) $(TOP)/src/elegant/Makefile
	$(MV) $(TOP)/src/elegant/Makefile.TMP $(TOP)/src/elegant/Makefile

elegant:
	$(RM) *.o O.$(EPICS_HOST_ARCH)/*.o
	$(MV) $(TOP)/src/elegant/Makefile $(TOP)/src/elegant/Makefile.TMP
	$(CP) $(TOP)/src/elegant/Makefile.Pelegant $(TOP)/src/elegant/Makefile
	$(MAKE) NOMPI=1 -f Makefile
	$(RM) $(TOP)/src/elegant/Makefile
	$(MV) $(TOP)/src/elegant/Makefile.TMP $(TOP)/src/elegant/Makefile

all buildInstall:
	$(RM) *.o O.$(EPICS_HOST_ARCH)/*.o O.$(EPICS_HOST_ARCH)/*.h
	$(MV) $(TOP)/src/elegant/Makefile $(TOP)/src/elegant/Makefile.TMP
	$(CP) $(TOP)/src/elegant/Makefile.Pelegant $(TOP)/src/elegant/Makefile
	$(MAKE) NOMPI=1 -f Makefile
	$(RM) $(TOP)/src/elegant/Makefile
	$(MV) $(TOP)/src/elegant/Makefile.TMP $(TOP)/src/elegant/Makefile

endif



ifdef BASE_3_15
clean:
else
clean::
endif
	$(RM) fitTraces.h vary.h twiss.h tune.h tuneFootprint.h trace.h subprocess.h steer_elem.h sliceAnalysis.h sdds_beam.h save_lattice.h sasefel.h run_rpnexpr.h response.h optimize.h optim_covariable.h matrix_output.h load_parameters.h link_elements.h frequencyMap.h floor.h error.h elegant.h ignoreElements.h transmuteElements.h divideElements.h correct.h steer_elem.h closed_orbit.h chrom.h bunched_beam.h aperture_search.h analyze.h amplif.h alter.h insertSCeffects.h insert_elements.h touschekScatter.h momentumAperture.h aperture_data.h replace_elements.h modulate.h ramp.h chaosMap.h ionEffects.h obstructionData.h global_settings.h



