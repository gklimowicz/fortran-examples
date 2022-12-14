.SUFFIXES:

HOMEDIR = $(wildcard ~)
MODELERC ?= $(HOMEDIR)/.modelErc
sinclude $(MODELERC)

export GISSCLIM_DIR := $(shell pwd)/../..
export CONFIG_DIR := $(GISSCLIM_DIR)/config

include $(CONFIG_DIR)/rules.mk

diags:
	sh ./compscr "$(F90) $(FFLAGS_ALL)" "$(LIBS)"
