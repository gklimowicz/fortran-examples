.SUFFIXES:

SUPPORTED_OPTIONS =

FSRCS = DOMAIN_DECOMPcs.f
F90SRCS=


ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
include $(CONFIG_DIR)/rules.mk

no_default:
	echo "no default target!!!"
	echo "GISSCLIM_DIR = $(GISSCLIM_DIR)"
	echo "CONFIG_DIR = $(CONFIG_DIR)"
	exit 1

