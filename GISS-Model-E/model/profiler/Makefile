#makefile for the component Ent
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF

FSRCS = 
#F90SRCS = $(wildcard *.F90)
F90SRCS= ProfileReport_mod.F90  TimeFormatUtilities_mod.F90  Timer_mod.F90 \
ReportColumn_mod.F90   TimerList_mod.F90            TimerPackage_mod.F90
F90OBJS = $(F90SRCS:%.F90=%.o)

ifeq ($(ESMF),YES)
$(F90OBJS) : F90FLAGS+=-DUSE_MPI
endif

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif


include $(CONFIG_DIR)/base.mk
include $(CONFIG_DIR)/rules.mk




