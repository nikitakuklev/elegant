#*************************************************************************
# Copyright (c) 2002 The University of Chicago, as Operator of Argonne
# National Laboratory.
# Copyright (c) 2002 The Regents of the University of California, as
# Operator of Los Alamos National Laboratory.
# This file is distributed subject to a Software License Agreement found
# in the file LICENSE that is included with this distribution. 
#*************************************************************************
#
# $Id: Makefile,v 1.4 2002-08-14 20:23:30 soliday Exp $
#
#  Lowest Level Directroy Makefile
# $Log: not supported by cvs2svn $
# Revision 1.3  1999/08/05 15:24:20  soliday
# Now uses Makefile.Host
#
# Revision 1.2  1997/10/20 14:57:02  borland
# Improved trace fitting and related routines.  Added output of traces
# after fitting.  Fixed some output-related bugs.
#
# Revision 1.1.1.1  1996/02/01  16:30:34  borland
# Imported files
#
#
#

TOP=../..
#EPICS=../../..
include $(TOP)/config/CONFIG_APPS

include $(TOP)/config/RULES_ARCHS
#include $(EPICS)/config/RULES_DIRS


