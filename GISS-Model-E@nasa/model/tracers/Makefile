# makefile for the component shared
.SUFFIXES:

SUPPORTED_OPTIONS = 

#F90SRCS = $(wildcard *.F90) # don't do this !
F90SRCS = \
EmissionRegion.F90        TracerBundle.F90        TracerHashMap.F90 \
EmissionScenario.F90      TracerBundleSubset.F90  TracerSource.F90 \
MonthlySurfaceSource.F90  TracerConstants.F90     TracerSurfaceSource.F90 \
OldTracer_mod.F90         Tracer.F90

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
include $(CONFIG_DIR)/rules.mk

