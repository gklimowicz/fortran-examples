# makefile for the component shared
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF

FSRCS = TRIDIAG.f


ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
#sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk

ifeq ($(OFFLINE_RUN),YES)
CPPFLAGS += -DOFFLINE_RUN
endif

