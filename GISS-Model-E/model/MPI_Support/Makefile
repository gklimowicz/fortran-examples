.SUFFIXES:

SUPPORTED_OPTIONS =

FSRCS = dd2d_utils.f
ifndef ($(FVCUBED),YES)
 FSRCS += DomainDecompLatLon.f
endif

#ifeq ($(FVCUBED),YES)
 FSRCS +=pario_fbsa.f
#endif

F90SRCS+= MpiSupport_mod.F90 dist_grid_mod.F90 Halo_mod.F90
F90SRCS+= GatherScatter_mod.F90 SpecialIO_mod.F90 GlobalSum_mod.F90
F90SRCS+= ProcessTopology_mod.F90 Domain_mod.F90 Hidden_mod.F90


F90SRCS+= DomainDecomposition_mod.F90


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

