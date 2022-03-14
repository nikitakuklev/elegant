#*************************************************************************
# Copyright (c) 2002 The University of Chicago, as Operator of Argonne
# National Laboratory.
# Copyright (c) 2002 The Regents of the University of California, as
# Operator of Los Alamos National Laboratory.
# This file is distributed subject to a Software License Agreement found
# in the file LICENSE that is included with this distribution. 
#*************************************************************************
#
# $Id: Makefile,v 1.14 2011-02-08 03:24:07 soliday Exp $
#
#
#

TOP=../../..
include $(TOP)/configure/CONFIG
include $(TOP)/src/elegant/sddsbrightness/Makefile.OAG
include $(TOP)/configure/RULES

sddsbrightness$(OBJ): ../sddsbrightness.h


ifeq ($(OS_CLASS),Darwin)
 F77 = gfortran -m64 -ffixed-line-length-132
 #F77 = /usr/local/bin/gfortran-4.2 -arch x86_64 -ffixed-line-length-132
 static_flags = -static-libgcc -L/opt/local/lib
endif
ifeq ($(HOST_ARCH),solaris)
 F77 = g77
 solaris_libs = -lrt
endif
ifeq ($(HOST_ARCH),solaris-sparc)
 F77 = g77
 solaris_libs = -lrt
endif

ifeq ($(EPICS_HOST_ARCH),solaris-sparc)
 F77 = g77
 solaris_libs = -lrt
 ifdef RATIONAL
  F77 = $(RATIONAL) f77 -g
 endif
endif

ifeq ($(EPICS_HOST_ARCH),solaris-sparc64)
 F77 = g77 -m64
 solaris_libs = -lrt
endif
ifeq ($(EPICS_HOST_ARCH),solaris-x86)
 F77 = f77
 solaris_libs = -lrt
endif
ifeq ($(EPICS_HOST_ARCH),solaris-x86-gnu)
 F77 = f77
 solaris_libs = -lrt
endif
ifeq ($(EPICS_HOST_ARCH),solaris-x86_64)
 F77 = f77 -m64
 solaris_libs = -lrt
endif

ifeq ($(EPICS_HOST_ARCH),linux-x86_64)
 F77 = gfortran -m64 -ffixed-line-length-132
 ifeq (cori, $(findstring cori, $(HOSTNAME)))
  F77 = ftn -m64 -ffixed-line-length-132
 endif
 ifeq (edison, $(findstring edison, $(HOSTNAME)))
  F77 = ftn -m64 -ffixed-line-length-132
 endif
endif

ifeq ($(EPICS_HOST_ARCH),linux-ppc64)
 F77 = gfortran -m64 -ffixed-line-length-132
endif
ifeq ($(EPICS_HOST_ARCH),linux-x86)
 F77 = gfortran -m32 -ffixed-line-length-132
endif
ifeq ($(EPICS_HOST_ARCH),win32-x86-mingw)
 F77 = gfortran -m32 -ffixed-line-length-132
endif
ifeq ($(EPICS_HOST_ARCH),windows-x64-mingw)
 F77 = x86_64-w64-mingw32-gfortran -m64 -ffixed-line-length-132
endif
ifeq ($(EPICS_HOST_ARCH),windows-x64)
 F77 = gfortran -m64 -ffixed-line-length-132
endif

#RATIONAL = purify
#RATIONAL = quantify
#RATIONAL = purecov
ifdef RATIONAL
 HOST_OPT = NO
 STATIC_BUILD=NO
 F77 := $(RATIONAL) -best-effort $(F77)
endif


ifeq ($(OS_CLASS),Linux)
 ifeq ($(STATIC_BUILD),YES)
  static_flags = -Wl,-Bstatic -static-libgcc
  linux_libs = -lc
 endif
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


sddsbrightness$(EXE): $(OBJNAMES)
	$(F77) -o $@ $(OBJNAMES) $(static_flags) -L$(INSTALL_LIB) -L$(EPICS_EXTENSIONS_LIB) -lmdbcommon -lmatlib -lfftpack -lSDDS1 -lrpnlib -lnamelist -lmdbmth -lmdblib -llzma -lz $(solaris_libs) $(linux_libs)

sddsurgent$(EXE): $(URGENT_OBJNAMES)
	$(F77) -o $@ $(URGENT_OBJNAMES) $(static_flags) -L$(INSTALL_LIB) -L$(EPICS_EXTENSIONS_LIB) -loagphy -lmdbcommon -lmatlib -lfftpack -lSDDS1 -lrpnlib -lnamelist -lmdbmth -lmdblib -llzma -lz $(solaris_libs) $(linux_libs)

sddsfluxcurve$(EXE): $(FLUXCURVE_OBJNAMES)
	$(F77) -o $@ $(FLUXCURVE_OBJNAMES) $(static_flags) -L$(INSTALL_LIB) -L$(EPICS_EXTENSIONS_LIB) -loagphy -lmdbcommon -lmatlib -lfftpack -lSDDS1 -lrpnlib -lnamelist -lmdbmth -lmdblib -llzma -lz $(solaris_libs) $(linux_libs)


