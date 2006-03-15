#*************************************************************************
# Copyright (c) 2002 The University of Chicago, as Operator of Argonne
# National Laboratory.
# Copyright (c) 2002 The Regents of the University of California, as
# Operator of Los Alamos National Laboratory.
# This file is distributed subject to a Software License Agreement found
# in the file LICENSE that is included with this distribution. 
#*************************************************************************
#
# $Id: Makefile,v 1.5 2006-03-15 17:19:24 shang Exp $
#
#
#

TOP=../../..
#If epics/extensions/configure directory exists, build with it.
#Otherwise use epics/extensions/config.
ifeq (0, $(words $(notdir $(wildcard $(TOP)/configure))))
include $(TOP)/config/CONFIG_APPS
include $(TOP)/config/RULES_ARCHS
else
include $(TOP)/configure/CONFIG
include $(TOP)/src/elegant/sddsbrightness/Makefile.OAG
include $(TOP)/configure/RULES

sddsbrightness$(OBJ): ../sddsbrightness.h


ifdef WIN32
usb$(OBJ) :
	f77 /compile_only /threads /libs:dll /fast ..\\usb.f

brighte$(OBJ) :
	f77 /compile_only /threads /libs:dll /fast ..\\brighte.f

hunt$(OBJ) :
	f77 /compile_only /threads /libs:dll /fast ..\\hunt.f

urgent$(OBJ):
	f77 /compile_only /threads /libs:dll /fast ..\\urgent.f

us$(OBJ):
	f77 /compile_only /threads /libs:dll /fast ..\\us.f
else
F77 = f77
ifdef AIX
F77 = g77
endif
ifeq ($(HOST_ARCH),solaris)
F77 = g77
solaris_libs = -lcurses -lrt
endif
ifeq ($(HOST_ARCH),solaris-sparc)
F77 = g77
solaris_libs = -lcurses -lrt
endif
ifeq ($(EPICS_HOST_ARCH),solaris-sparc)
F77 = g77
solaris_libs = -lcurses -lrt
ifdef RATIONAL
F77 = $(RATIONAL) f77 -g
endif
endif
ifeq ($(EPICS_HOST_ARCH),solaris-sparc64)
F77 = g77 -m64
solaris_libs = -lcurses -lrt
endif
usb$(OBJ) : 
	$(F77) ../usb.f -c -o usb$(OBJ)	
brighte$(OBJ) :
	$(F77) ../brighte.f -c -o brighte$(OBJ)
hunt$(OBJ) :
	$(F77) ../hunt.f -c -o hunt$(OBJ)
urgent$(OBJ) :
	$(F77) ../urgent.f -c -o urgent$(OBJ)
us$(OBJ) :
	$(F77) ../us.f -c -o us$(OBJ)
endif

ifdef WIN32
sddsbrightness$(EXE) : $(OBJNAMES)
	f77 /exe /threads /libs:dll /link /WARN:3 /incremental:no /opt:ref /release usb.obj \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdbcommon.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\matlib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\fftpack.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\SDDS1.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\rpnlib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\namelist.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdbmth.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdblib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\z.lib \
	-out:sddsbrightness$(EXE) $(OBJNAMES)
sddsurgent$(EXE) : $(URGENT_OBJNAMES)
	f77 /exe /threads /libs:dll /link /WARN:3 /incremental:no /opt:ref /release urgent.obj us.obj \
	c:\\epics\\extensions\\lib\\$(T_A)\\oagphy.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdbcommon.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\matlib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\fftpack.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\SDDS1.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\rpnlib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\namelist.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdbmth.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdblib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\z.lib \
	-out:sddsbrightness$(EXE) $(OBJNAMES)
else
sddsbrightness$(EXE): $(OBJNAMES)
	$(F77) -o $@ $(OBJNAMES) -L$(INSTALL_LIB) -L$(EPICS_EXTENSIONS_LIB) -lmdbcommon -lmatlib -lfftpack -lSDDS1 -lrpnlib -lnamelist -lmdbmth -lmdblib -lz $(solaris_libs)

sddsurgent$(EXE): $(URGENT_OBJNAMES)
	$(F77) -o $@ $(URGENT_OBJNAMES) -L$(INSTALL_LIB) -L$(EPICS_EXTENSIONS_LIB) -loagphy -lmdbcommon -lmatlib -lfftpack -lSDDS1 -lrpnlib -lnamelist -lmdbmth -lmdblib -lz $(solaris_libs)
endif

endif
